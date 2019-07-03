#ifndef FANOMAN_STIELTJES_INTERPOLANT_H
#define FANOMAN_STIELTJES_INTERPOLANT_H

#include "printer.h"

#include <exception>
#include <sstream>
#include <vector>
#include <string>
#include <boost/multiprecision/number.hpp> // fabs, pow


namespace fanoman {

	/** Implementation of a piecewise monotonicity preserving cubic hermite spline
	 * following
	 *
	 * James M. Hyman, SIAM J. Sci. Stat. Comput. 1983, 4, 645-654 {hyman:1983}.
	 *
	 * Equation and page numbers given here refer to the numbers in this
	 * publication.
	 *
	 * Other publications eventually referenced below:
	 *
	 * Randall L. Dougherty, Alan Edelman, James M. Hyman, Math. Comp. 52, 1989,
	 * 471-494 {dougherty:1989}.
	*/

template<class number_t>
class interpolant {
private:
	std::vector<number_t> m_x;
	std::vector<number_t> m_y;
	printer &m_printer;
	const std::string m_derivative_method;
	const std::string m_constrain_method;

	bool m_skip_monotonicity_constrain;
	bool m_remove_redundant_grid_points;

	size_t m_length;
	std::vector<number_t> m_dx;
	std::vector<number_t> m_dy;
	std::vector<number_t> m_s;
	std::vector<number_t> m_f_dot;
	std::vector<number_t> m_c2;
	std::vector<number_t> m_c3;
	std::vector<number_t> m_c4;

public:
	interpolant(
		const std::vector<number_t> &x_,
		const std::vector<number_t> &y_,
		printer &printer_,
		const std::string derivative_method_ = "fritsch-carlson",
		const std::string constrain_method_ = "",
		bool remove_redundant_grid_points_ = false) :
		m_x(x_), m_y(y_), m_printer(printer_),
		m_derivative_method(derivative_method_),
		m_constrain_method(constrain_method_),
		m_remove_redundant_grid_points(remove_redundant_grid_points_),
		m_skip_monotonicity_constrain(false) {
		init();
	}

	number_t interpolate(const number_t &x0, bool &is_extrapolated);
	std::vector<number_t> get_c2_coefficients() { return m_c2; }

private:
	void init();
	void check_input();
	void init_diffs_and_slopes();

	void compute_derivatives();
	void compute_fritsch_butland_derivatives();
	void compute_parabolic_derivatives();
	void compute_fritsch_carlson_derivatives();

	void compute_c2_coefficients();
	void compute_hyman_1983_c2_coefficients();
	void compute_hyman_1989_c2_coefficients();

	void compute_c3_c4_coefficients();

	int sign(number_t val);

};


template<class number_t>
void interpolant<number_t>::init() {

	check_input();
	init_diffs_and_slopes();
	compute_derivatives();
	compute_c2_coefficients();
	compute_c3_c4_coefficients();
}


template<class number_t>
void interpolant<number_t>::check_input() {

	// checks input vectors for consistent sizes

	bool do_crash = false;

	if (m_x.size() != m_y.size()) {
		std::ostringstream oss;
		oss << "Error:  got x and y vectors of differing length: x.size() = "
			<< m_x.size() << ", y.size() = " << m_y.size();
		m_printer.error(oss.str());
		do_crash = true;
 	}

	if (do_crash) {
		throw std::runtime_error("Error in fanoman::interpolant::check_input(): See above for details.");
	}

	if (m_remove_redundant_grid_points) {

		// remove (x_i,y_i) pairs if x_i == x_{i-1} if requested
		std::vector<number_t> x_, y_;
		x_.push_back(m_x.front());
		y_.push_back(m_y.front());
		for (size_t i = 1; i != m_x.size(); ++i) {
			if (m_x[i] != x_.back()) {
				x_.push_back(m_x[i]);
				y_.push_back(m_y[i]);
			} 
		}
		m_x = x_;
		m_y = y_;
	}

	// set m_length to the size of the possibly modified m_x vector
	m_length = m_x.size();
}


template<class number_t>
void interpolant<number_t>::init_diffs_and_slopes() {

	// calculate the slope between input points
	for (size_t i = 0; i != m_length-1; ++i) {
		m_dx.push_back(m_x[i+1] - m_x[i]);
		m_dy.push_back(m_y[i+1] - m_y[i]);
		m_s.push_back(m_dy.back() / m_dx.back());
	}
}


template<class number_t>
void interpolant<number_t>::compute_derivatives() {

	if (m_derivative_method == "fritsch-butland") {
		compute_fritsch_butland_derivatives();
		return;
	}
	if (m_derivative_method == "parabolic") {
		compute_parabolic_derivatives();
		return;
	}
	if (m_derivative_method == "fritsch-carlson") {
		m_skip_monotonicity_constrain = true;
		compute_fritsch_carlson_derivatives();
		return;
	}

	// if we are here, the requested method for the computation of approximate
	// derivatives is not defined;  throw an error and exit
	std::ostringstream oss1, oss2;
	oss1 << "Don't know this method for the computation of approximate derivatives";
	oss2 << "in fanoman::interpolant::compute_derivatives(): '"
		<< m_derivative_method << "'.";
	m_printer.error(oss1.str());
	m_printer.error(oss2.str());
	throw std::runtime_error("Error in fanoman::interpolant::compute_derivatives(): See message"
		" above for details.");
}


template<class number_t>
void interpolant<number_t>::compute_fritsch_butland_derivatives() {

	// compute the approximated derivatives following the Fritsch-Butland
	// algorithm as stated in {hyman:1983}, table 1, entry 2 on p. 649, which
	// in turn refers to 
	//
	// F. N. Fritsch, J. Butland: An improved monotone piecewise cubic
	// interpolation algorithm. Lawrence Livermore National Laboratory preprint
	// UCRL-85104, 1980.
	//
	// On the boundaries, the formulae given on p. 650 in section "C. Boundaries"
	// are employed.
	//
	// The computed approximate derivatives are stored in m_f_dot.

	// the lower boundary
	number_t f_dot_low = (m_dx[0] + m_dx[0] + m_dx[1]) * m_s[0] - m_dx[0] * m_s[1];
	f_dot_low /= (m_dx[0] + m_dx[1]);
	m_f_dot.push_back(f_dot_low);

	// values between boundaries
	for (size_t i = 1; i != m_length-1; ++i) {
		number_t s_i_min = std::min<number_t>(m_s[i-1], m_s[i]);
		number_t s_i_max = std::max<number_t>(m_s[i-1], m_s[i]);
		number_t f_dot_val = number_t(3.) * s_i_min * s_i_max;
		f_dot_val /= (s_i_max + s_i_min + s_i_min);
		m_f_dot.push_back(f_dot_val);
	}

	// the high boundary
	number_t f_dot_high = (m_dx[m_length-2] + m_dx[m_length-2] + m_dx[m_length-3])
		* m_s[m_length-2] - m_dx[m_length-2] * m_s[m_length-3];
	f_dot_high /= (m_dx[m_length-2] + m_dx[m_length-3]);
	m_f_dot.push_back(f_dot_high);

}


template<class number_t>
void interpolant<number_t>::compute_parabolic_derivatives() {

	// compute the approximated derivatives following the Fritsch-Butland
	// algorithm as stated in {hyman:1983}, table 1, entry 3 on p. 649, which
	// in turn refers to 
	//
	// The numercial solution of time dependent PDEs on an adaptive mesh, Los
	// Alamos Scientific Laboratory, Los Alamos, NM, LA-UR-80-3702, 1980.
	//
	// On the boundaries, the formulae given on p. 650 in section "C. Boundaries"
	// are employed.
	//
	// The computed approximate derivatives are stored in m_f_dot.

	// the lower boundary
	number_t f_dot_low = (m_dx[0] + m_dx[0] + m_dx[1]) * m_s[0] - m_dx[0] * m_s[1];
	f_dot_low /= (m_dx[0] + m_dx[1]);
	m_f_dot.push_back(f_dot_low);

	// values between boundaries
	for (size_t i = 1; i != m_length-1; ++i) {
		number_t f_dot_val = m_dx[i-1] * m_s[i] + m_dx[i] * m_s[i-1];
		f_dot_val /= (m_x[i+1] - m_x[i-1]);
		m_f_dot.push_back(f_dot_val);
	}

	// the high boundary
	number_t f_dot_high = (m_dx[m_length-2] + m_dx[m_length-2] + m_dx[m_length-3])
		* m_s[m_length-2] - m_dx[m_length-2] * m_s[m_length-3];
	f_dot_high /= (m_dx[m_length-2] + m_dx[m_length-3]);
	m_f_dot.push_back(f_dot_high);

}


template<class number_t>
void interpolant<number_t>::compute_fritsch_carlson_derivatives() {

	// compute the approximated derivatives following the Fritsch-Carlson
	// algorithm as given in
	//
	// F. N. Fritsch, R. E. Carlson, SIAM J. Numer. Anal. 17, 1980, 238-246
	// {fritsch:1980}.
	//
	// The computed approximate derivatives are stored in m_f_dot.

	// the lower boundary
	number_t f_dot_low = m_s[0];
	m_f_dot.push_back(f_dot_low);

	// values between boundaries
	for (size_t i = 1; i != m_length-1; ++i) {
		number_t f_dot_val;
		if ((m_s[i] * m_s[i-1]) <= number_t(0.)) {
			f_dot_val = 0.;
		} else {
			f_dot_val = number_t(3.) * (m_dx[i] + m_dx[i-1]);
			f_dot_val /= (
				(m_dx[i] + m_dx[i] + m_dx[i-1]) / m_s[i-1]
				+ (m_dx[i] + m_dx[i-1] + m_dx[i-1]) / m_s[i]);
		}
		m_f_dot.push_back(f_dot_val);
	}

	// the high boundary
	number_t f_dot_high = m_s[m_length-2];
	m_f_dot.push_back(f_dot_high);

}


template<class number_t>
void interpolant<number_t>::compute_c2_coefficients() {

	if (m_skip_monotonicity_constrain || m_constrain_method.empty()) {

		// only copy inital approximate derivatives to c2 coefficient vector
		m_c2 = m_f_dot;
		return;
	}

	if (m_constrain_method == "hyman-1983") {
		compute_hyman_1983_c2_coefficients();
		return;
	}
	if (m_constrain_method == "hyman-1989") {
		compute_hyman_1989_c2_coefficients();
		return;
	}

	// if we are here, the requested method for the monotonicity constraining
	// is not defined;  throw an error and exit
	std::ostringstream oss1, oss2;
	oss1 << "Don't know this method for applying monotonicity constraints";
	oss2 << "in fanoman::interpolant::compute_c2_coefficients(): '"
		<< m_constrain_method << "'.";
	m_printer.error(oss1.str());
	m_printer.error(oss2.str());
	throw std::runtime_error("Error in fanoman::interpolant::compute_c2_coefficients(): See message"
		" above for details.");
}


template<class number_t>
void interpolant<number_t>::compute_hyman_1983_c2_coefficients() {

	// computes the c2 coefficients as given in {hyman:1983}, eq. 2.6

	// compute the c2 coefficients following eq. 2.6
	for (size_t i = 0; i != m_length; ++i) {
		// the loop control variable i corresponds to n+\frac12
		number_t s_i_minus_half, s_i_plus_half;
		if (i == 0) {
			// first element, take S_{-\frac12} = S_{\frac12}
			s_i_plus_half = m_s[i];
			s_i_minus_half = s_i_plus_half;
		} else {
			if (i == m_length-1) {
				// last element, take S_{n+\frac12} = S_{n-\frac12}
				s_i_minus_half = m_s[i-1];
				s_i_plus_half = s_i_minus_half;
			} else {
				// elements in between
				s_i_plus_half = m_s[i];
				s_i_minus_half = m_s[i-1];
			}
		}

		// distinguish between the two cases of eq. 2.6
		if (m_f_dot[i] >= number_t(0.)) {
			// first line of eq. 2.6, i. e. \sigma \geq 0
			number_t c_val = std::min<number_t>(
				std::max<number_t>(number_t(0.), m_f_dot[i]),
				number_t(3.) * std::min<number_t>(
				boost::multiprecision::fabs(s_i_minus_half),
				boost::multiprecision::fabs(s_i_plus_half)));
			m_c2.push_back(c_val);
		} else {
			// second line of eq. 2.6, i. e. \sigma < 0
			number_t c_val = std::max<number_t>(
				std::min<number_t>(number_t(0.), m_f_dot[i]),
				number_t(-3.) * std::min<number_t>(
				boost::multiprecision::fabs(s_i_minus_half),
				boost::multiprecision::fabs(s_i_plus_half)));
			m_c2.push_back(c_val);
		}
	}
}


template<class number_t>
void interpolant<number_t>::compute_hyman_1989_c2_coefficients() {

	// computes the c2 coefficients as stated in {dougherty:1989},
	// using eqs. 4.4 and algorithm on p. 478

	// handle low bound
	number_t c2_low;
	if (sign(m_f_dot[0]) == sign(m_s[0])) {
		c2_low = number_t(sign(m_f_dot[0])) * std::min<number_t>(
			boost::multiprecision::fabs(m_f_dot[0]),
			number_t(3.) * boost::multiprecision::fabs(m_s[0]));
	} else {
		c2_low = number_t(0.);
	}
	m_c2.push_back(c2_low);

	// care about all values between i=1 and including i=m_length-2
	for (size_t i = 1; i != m_length-1; ++i) {

		number_t p_i_0 = m_s[i-1] * m_dx[i] + m_s[i] * m_dx[i-1];
		p_i_0 /= (m_dx[i-1] + m_dx[i]);

		number_t p_i_minus_1, p_i_plus_1;
		if (i > 1) {
			p_i_minus_1 = m_s[i-1] * (m_dx[i-1] + m_dx[i-1] + m_dx[i-2]);
			p_i_minus_1 -= m_s[i-2] * m_dx[i-1];
			p_i_minus_1 /= (m_dx[i-2] + m_dx[i-1]);
		}
		if (i < m_length-2) {
			p_i_plus_1 = m_s[i] * (m_dx[i] + m_dx[i] + m_dx[i+1]);
			p_i_plus_1 -= m_s[i+1] * m_dx[i];
			p_i_plus_1 /= (m_dx[i] + m_dx[i+1]);
		}

		number_t m_i = number_t(3.) * std::min<number_t>(
			boost::multiprecision::fabs(p_i_0), std::min<number_t>(
			boost::multiprecision::fabs(m_s[i-1]),
			boost::multiprecision::fabs(m_s[i])));

		if (i > 1) {
			// check if the expressions in the sum all have the same sign
			number_t abs_sum = boost::multiprecision::fabs(
				p_i_0 + p_i_minus_1 + (m_s[i-1] - m_s[i-2]) + (m_s[i] - m_s[i-1]));
			number_t sum_abs =
				boost::multiprecision::fabs(p_i_0)
				+ boost::multiprecision::fabs(p_i_minus_1)
				+ boost::multiprecision::fabs(m_s[i-1] - m_s[i-2])
				+ boost::multiprecision::fabs(m_s[i] - m_s[i-1]);
			if (abs_sum == sum_abs) {
				// all numbers have the same sign
				m_i = std::max<number_t>(m_i, number_t(1.5) * std::min<number_t>(
					boost::multiprecision::fabs(p_i_0),
					boost::multiprecision::fabs(p_i_minus_1)));
			}
		}
		
		if (i < m_length-2) {
			// check if the expressions in the sum all have the same sign
			number_t abs_sum = boost::multiprecision::fabs(
				(m_s[i] - m_s[i-1]) + (m_s[i+1] - m_s[i]) - p_i_0 - p_i_plus_1);
			number_t sum_abs =
				boost::multiprecision::fabs(-p_i_0)
				+ boost::multiprecision::fabs(-p_i_plus_1)
				+ boost::multiprecision::fabs(m_s[i] - m_s[i-1])
				+ boost::multiprecision::fabs(m_s[i+1] - m_s[i]);
			if (abs_sum == sum_abs) {
				// all numbers have the same sign
				m_i = std::max<number_t>(m_i, number_t(1.5) * std::min<number_t>(
					boost::multiprecision::fabs(p_i_0),
					boost::multiprecision::fabs(p_i_plus_1)));
			}
		}

		number_t c2_val(0.);
		if (sign(m_f_dot[i]) == sign(p_i_0)) {
			c2_val = number_t(sign(m_f_dot[i])) * std::min<number_t>(
				boost::multiprecision::fabs(m_f_dot[i]), m_i);
		}

		m_c2.push_back(c2_val);
	}

	// handle high bound
	number_t c2_high;
	if (sign(m_f_dot[m_length-1]) == sign(m_s[m_length-2])) {
		c2_high = number_t(sign(m_f_dot[m_length-1])) * std::min<number_t>(
			boost::multiprecision::fabs(m_f_dot[m_length-1]),
			number_t(3.) * boost::multiprecision::fabs(m_s[m_length-2]));
	} else {
		c2_high = number_t(0.);
	}
	m_c2.push_back(c2_high);

}


template<class number_t>
void interpolant<number_t>::compute_c3_c4_coefficients() {

	// compute c3 and c4 coefficients following expressions below eq. 2.1
	for (size_t i = 0; i != m_length-1; ++i) {

		number_t c3_value;
		c3_value = number_t(3.) * m_s[i] - m_c2[i+1] - number_t(2.) * m_c2[i];
		c3_value /= m_dx[i];
		m_c3.push_back(c3_value);

		number_t c4_value;
		c4_value = m_c2[i+1] + m_c2[i] - number_t(2.) * m_s[i];
		c4_value /= boost::multiprecision::pow(m_dx[i], 2);
		m_c4.push_back(c4_value);
	}

}


template<class number_t>
number_t interpolant<number_t>::interpolate(
		const number_t &x0, bool &is_extrapolated) {

	is_extrapolated = false;
	size_t k = 0;

	// check the bounds
	if (x0 < m_x.front() || x0 > m_x.back()) {
		// we have to extrapolate the y0 value because x0 is not in the range
		// covered by m_x
		if (x0 < m_x.front()) {
			// extrapolate below the low end of m_x
			k = 0;
		} else {
			// extrapolate above the high end of m_x
			k = m_length - 2;
		}
		is_extrapolated = true;

		{
			std::ostringstream oss1, oss2, oss3;
			oss1 << "Warning: the passed x value to interpolate an y value for lies";
			oss2 << "outside the grid: x0 = " << x0 << " not in [" << m_x.front()
				<< "," << m_x.back() << "]";
			oss3 << "The returned value will be extrapolated!";
			m_printer.debug(4, oss1.str());
			m_printer.debug(4, oss2.str());
			m_printer.debug(4, oss3.str());
		}
	} else {
		// we don't have to extrapolate
		// find the index k, for which x_k <= x0 <= x_{k+1}
		for (; k != m_length-1; ++k) {
	
			if (m_x[k] <= x0 && x0 <= m_x[k+1]) {
				// since we checked above that x0 lies within the defined grid of m_x
				// values, the following break is always hit;  The value then stored in k
				// is the index looked for
				break;
			}
		}
	}

	// calculate the interpolated y value at x0 using eq. 2.1
	number_t y0 = m_y[k];
	y0 += (x0 - m_x[k]) * m_c2[k];
	y0 += boost::multiprecision::pow((x0 - m_x[k]), 2) * m_c3[k];
	y0 += boost::multiprecision::pow((x0 - m_x[k]), 3) * m_c4[k];

	return y0;
}


template<class number_t>
int interpolant<number_t>::sign(number_t val) {

	return (number_t(0.) < val) - (val < number_t(0.));
}


} // namespace fanoman


#endif // FANOMAN_STIELTJES_INTERPOLANT_H

