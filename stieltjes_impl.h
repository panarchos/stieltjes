#ifndef FANOMAN_STIELTJES_STIELTJES_IMPL_H
#define FANOMAN_STIELTJES_STIELTJES_IMPL_H

#include <exception>
#include <sstream>
#include <algorithm> // for std::swap, until C++11
#include <utility>   // for std::swap, starting from C++11

#include <boost/multiprecision/number.hpp> // will provide fabs, sqrt and pow

#include "ql.h"
#include "interpolant.h"
//#include "util.h"


namespace fanoman {


template<class number_t>
void stieltjes<number_t>::init() {

	check_input();
	init_e_bounds();

	// set the maximal available/generated order of Q polynomials
	m_max_poly_order = std::min(m_max_allowed_poly_order, m_num_points);
}


template<class number_t>
void stieltjes<number_t>::check_input() {

	bool do_crash = false;

	// check if the numbers of passed energy differences and gammas match
	size_t n_e = m_e_diffs.size();
	size_t n_g = m_gammas.size();

	if (n_e != n_g) {
		std::ostringstream oss1, oss2, oss3, oss4;
		oss1 << "Error in stieltjes<number_t>::check_input():";
		oss2 << "Numbers of passed energy and gamma points do not match:";
		oss3 << "Number of energy points passed: " << n_e;
		oss4 << "Number of gamma points passed: " << n_g;
		m_printer.error(oss1.str());
		m_printer.error(oss2.str());
		m_printer.error(oss3.str());
		m_printer.error(oss4.str());

		do_crash = true;
	}

	// set the number of input points
	m_num_points = n_e;

	// check if the passed number of points is in the allowed range
	if (m_num_points < m_min_num_points) {
		std::ostringstream oss1, oss2, oss3;
		oss1 << "Error in stieltjes<number_t>::check_input():";
		oss2 << "Got less energy points than needed.  Minimum number of points"
			<< " required: " << m_min_num_points;
		oss3 << "Number of points passed: " << m_num_points;
		m_printer.error(oss1.str());
		m_printer.error(oss2.str());
		m_printer.error(oss3.str());

		do_crash = true;
	}
	if (m_num_points > m_max_num_points) {
		std::ostringstream oss1, oss2, oss3;
		oss1 << "Error in stieltjes<number_t>::check_input():";
		oss2 << "Got more energy points than allowed.  Maximum number of points"
			<< " allowed: " << m_max_num_points;
		oss3 << "Number of points passed: " << m_num_points;
		m_printer.error(oss1.str());
		m_printer.error(oss2.str());
		m_printer.error(oss3.str());

		do_crash = true;
	}

	if (do_crash)
		throw std::runtime_error("Error in stieltjes<number_t>::check_input(): See message above for details.");

}


template<class number_t>
void stieltjes<number_t>::init_e_bounds() {
	// run through the energy differences passed to the constructor and 
	// store the extrema in m_e_min and m_e_max, respectively

	m_e_min = m_e_diffs[0];
	m_e_max = m_e_diffs[0];

	for (typename std::vector<number_t>::iterator e_it = m_e_diffs.begin();
			e_it != m_e_diffs.end(); ++e_it) {

		m_e_min = std::min(m_e_min, *e_it);
		m_e_max = std::max(m_e_max, *e_it);
	}

	{
		std::ostringstream oss;
		oss << "Got emax: " << m_e_max << ", emin: " << m_e_min;
		m_printer.debug(5, oss.str());
	}

}


template<class number_t>
double stieltjes<number_t>::compute() {
	// compute the decay width gamma according to:
	//
	// Florian M{\"}uller-Plathe, Geerd H. F. Diercksen, Molecular photoionization
	// cross sections by moment theory. An introduction, in: Sylvio Canuto,
	// Jos{\'e} D'Albuquerque e Castro, Fernando J. Pix{\~}ao (Eds.), Electronic
	// structure of Atoms, Molecules and Solids, Proceedings of the II Escola
	// Brasiliera de Estruture Eletr{\hat}onica, Olinda, Brazil, July 17-22, 1989
	// (Workd Scientific, Singapore, 1990), p. 1-29 {mueller-plathe:1990}
	//

	if (m_e_req < m_e_min || m_e_req > m_e_max) {
		// return the mean gamma value because the stieltjes procedure will not work
		// if the requested energy lies outside the range of input energies
		return mean_gamma();
	}

	// shift energies to avoid problems arising from small denominators
	shift_energies();

	// initiate the a and b coefficients following {mueller-plathe:1990}, eqs.
	// 3.3.20-3.3.23
	// a coefficients run from 1 to m_max_allowed_poly_order, a_0 will be left
	// 0.0; one could shift the index of a so that it runs from 0 to
	// m_max_allowed_poly_order-1; however, not doing so preserves the indices as
	// given in literature, leading to easier comparable code
	std::vector<number_t> a_coeffs(m_max_allowed_poly_order+1, 0.);
	// b coefficients run from 0 to m_max_allowed_poly_order-1
	std::vector<number_t> b_coeffs(m_max_allowed_poly_order, 0.);
	// we will generate the m_max_allowed_poly_order'th polynomial for orthogonality
	// checking, so we need orders from 0 to m_max_allowed_poly_order, i. e.
	// m_max_allowed_poly_order+1 polynomials
	std::vector< std::vector<number_t> > q_polynomials(
		m_max_allowed_poly_order+1, std::vector<number_t>(m_num_points, 0.));
	generate_polynomials(a_coeffs, b_coeffs, q_polynomials);

	// determine which of the generated polynomials are usable
	// in addition, get_polynomial_orders will set m_do_convergence_check
	// according to the maximal usable polynomial order
	std::pair<size_t, size_t> usable_poly_orders;
	get_polynomial_orders(q_polynomials, usable_poly_orders);

	std::vector<number_t> gamma_histogram(usable_poly_orders.second+1);
	calc_gamma_histogram(a_coeffs, b_coeffs,
		q_polynomials, usable_poly_orders, gamma_histogram);

	number_t gamma = find_converged_gamma(gamma_histogram, usable_poly_orders);

	{
		std::ostringstream oss;
		oss << "Got final value for gamma: " << gamma;
		m_printer.result(oss.str());
	}
	return gamma.template convert_to<double>();
} 


template<class number_t>
number_t stieltjes<number_t>::find_converged_gamma(
		std::vector<number_t> &gamma_histogram,
		std::pair<size_t, size_t> &usable_poly_orders) {

	if (!m_do_convergence_check) {
		// only return the highest order gamma
		{
			m_printer.info("Not going for a convergence check of gamma because there are");
			m_printer.info("not enough polynomial orders available.");
			std::ostringstream oss;
			oss << "Available polynomial orders: " << usable_poly_orders.first
				<< " to " << usable_poly_orders.second;
			m_printer.info(oss.str());
		}
		return gamma_histogram[usable_poly_orders.second];
	}

	// do the convergence check and return appropriate gamma
	// m_do_convergence_check was only set to true in get_polynomial_orders()
	// if there are at least 3 orders available
	{
		m_printer.info("Going for a convergence check of gamma.");
		std::ostringstream oss;
		oss << "Available polynomial orders: " << usable_poly_orders.first
			<< " to " << usable_poly_orders.second;
		m_printer.info(oss.str());
	}

	number_t conv_threshold = m_conv_threshold_initial;
	for (size_t conv_iter = 0; conv_iter != m_conv_max_iter;
			++conv_iter, conv_threshold *= m_conv_threshold_progress) {

		// try to converge within m_conv_max_iter iterations, while softening the
		// convergence criterion at the end of each iteration by the factor
		// m_conv_threshold_progress
		for (size_t n_max = usable_poly_orders.second;
				n_max >= usable_poly_orders.first + 2; --n_max) {
			// run through available orders so that we can always check convergence
			// using order n_max and the next two lower ones
			number_t gamma_diff_12 = boost::multiprecision::fabs(
				gamma_histogram[n_max] - gamma_histogram[n_max-1]);
			number_t gamma_diff_13 = boost::multiprecision::fabs(
				gamma_histogram[n_max] - gamma_histogram[n_max-2]);
			number_t gamma_diff_23 = boost::multiprecision::fabs(
				gamma_histogram[n_max-1] - gamma_histogram[n_max-2]);
			number_t gamma_diff_mean = (
				gamma_diff_12 + gamma_diff_13 + gamma_diff_23) / number_t(3.);
			number_t gamma_mean = (
				gamma_histogram[n_max] + gamma_histogram[n_max-1] + gamma_histogram[n_max-2]) / number_t(3.);

			// set up the convergence criterion
			number_t conv_criterion = gamma_mean * conv_threshold;
			// modify the convergence criterion so that higher/lower orders are
			// preferred according to the value of m_conv_fac; 3 possible cases:
			// m_conv_fac == 1.0 -> nothing happens, convergence criterion does not
			//                      depend on the value of n_max
			// m_conv_fac  > 1.0 -> the smaller the value of n_max, the higher the
			//                      value of conv_criterion, so preference for lower
			//                      orders
			// m_conv_fac  < 1.0 -> the higher the value of n_max, the lower the
			//                      value of conv_criterion, so preference for higher
			//                      orders
			conv_criterion *= boost::multiprecision::pow(
				m_conv_fac, usable_poly_orders.second-n_max);
			if (gamma_diff_mean <= conv_criterion) {
				// we are converged, so return the converged value gamma_mean
				{
					std::ostringstream oss1, oss2;
					oss1 << "Detected convergence of gamma in " << n_max << "th order.";
					oss2 << "Convergence threshold in iteration " << conv_iter+1
						<< ": " << 100. * conv_threshold << "%";
					m_printer.result(oss1.str());
					m_printer.result(oss2.str());
				}
				return gamma_mean;
			}
		}
	}

	// if we are here, no convergence could be found within m_conv_max_iter
	// iterations.  If so, we return the mean value of the three highest available
	// orders as final gamma
	{
		std::ostringstream oss;
		oss << "No convergence found within " << m_conv_max_iter << " iterations.";
		m_printer.result(oss.str());
		m_printer.result("Using the mean gamma value of the three highest available orders.");
	}
	number_t gamma(0.);
	for (size_t i = 0; i <= 2; ++i)
		gamma += gamma_histogram[usable_poly_orders.second-i];
	gamma /= number_t(3.);

	return gamma;
}


template<class number_t>
void stieltjes<number_t>::calc_gamma_histogram(
		std::vector<number_t> &a_coeffs,
		std::vector<number_t> &b_coeffs,
		std::vector< std::vector<number_t> > &q_polynomials,
		std::pair<size_t, size_t> &usable_poly_orders,
		std::vector<number_t> &gamma_histogram) {

	for (size_t n_max = usable_poly_orders.first;
			n_max <= usable_poly_orders.second; ++n_max) {

		// fill in diagonal elements: diag[n] = a[n+1]
		std::vector<number_t> diag(n_max, 0.);
		for (size_t n = 0; n != n_max; ++n) {

			diag[n] = a_coeffs[n+1];
		}

		// fill in off-diagonal elements: offdiag[0] = 0, offdiag[n!=0] = -sqrt(b[n])
		std::vector<number_t> offdiag(n_max, 0.);
		for (size_t n = 1; n != n_max; ++n) {

			offdiag[n] = number_t(-1.) * boost::multiprecision::sqrt(b_coeffs[n]);
		}

		// set up diagonal n_max x n_max matrix (i. e. vector of vector)
		std::vector< std::vector<number_t> > ab_vec(n_max,
			std::vector<number_t>(n_max, 0.));
		for (size_t i = 0; i != n_max; ++i) {
			ab_vec[i][i] = 1.;
		}

		//int tql_error_status;
		//// This QL implementation returns the eigenvectors in ab_vec as vector of
		//// eigenvectors, i. e. the i'th element in ab_vec represents the i'th
		//// eigenvector corresponding to the i'th eigenvalue which is stored in diag
		//// as i'th element
		//tql2(diag, offdiag, ab_vec, tql_error_status);

		// This QL implementation returns the eigenvectors in ab_vec as vector of
		// eigenvectors, i. e. the i'th element in ab_vec represents the i'th
		// eigenvector corresponding to the i'th eigenvalue which is stored in diag
		// as i'th element
		int tql_error_status = ql<number_t>(m_printer).diagonalize(diag, offdiag, ab_vec);

		{
			// print out diag, offdiag, ab_vec for debugging
			m_printer.debug(6, "Got result from QL:");
			print_vec(diag, "The diagonal: ", 0, 6);
			size_t row_cnt = 0;
			for (typename std::vector< std::vector<number_t> >::iterator
					it = ab_vec.begin(); it != ab_vec.end(); ++it, ++row_cnt) {
				std::ostringstream oss;
				oss << "ab_vec row " << row_cnt << ": ";
				print_vec(*it, oss.str(), 0, 6);
			}
		}

		if (tql_error_status > 0) {
			std::ostringstream oss;
			oss << "QL: Failed to converge eigenvalue no. " << tql_error_status;
			m_printer.warning(oss.str());
		}
		if (tql_error_status < 0) {
			std::ostringstream oss;
			oss << "QL: Can't handle parameter no. " << tql_error_status
				<< " passed to stieltjes::tql2().";
			m_printer.warning(oss.str());
		}

		// transform the results of ql algorithm;  because the eigenvalues returned
		// from the QL procedure are inverse energies, sorting doesn't make sense if
		// there are both positive and negative eigenvalues because the ordering
		// would change again once the inverse value is formed.  Thus, sorting is
		// done again after setting up the energy and gamma vectors
		std::vector<number_t> e_vec(n_max);
		std::vector<number_t> gamma_vec(n_max);
		for (size_t i = 0; i != n_max; ++i) {

			e_vec[i] = number_t(1.) / diag[i];
			// calculate gamma (oscillator strength) values using eq. 3.4.17
			gamma_vec[i] = b_coeffs[0] * ab_vec[i][0] * ab_vec[i][0]; 
		}

		// when inverting energies, the ordering may get lost if there are positive
		// and negative elements in e_vec
		// so we sort the energy values, and according to these, the gamma values
		// again
		sort_ql_energies_and_gammas(e_vec, gamma_vec);

		// print out some additional information, i. e. the cumulative oscillator
		// strength distribution step function
		number_t cumulative_gamma = 0.;
		std::ostringstream oss_step;
		oss_step << "step-wise cumulative oscillator strength distribution for"
			<< " order " << n_max;
		m_printer.announce(oss_step.str(), '-', 6);
		for (size_t i = 0; i < n_max; ++i) {

			cumulative_gamma += gamma_vec[i];
			std::ostringstream oss_step_i;
      oss_step_i << e_vec[i] - m_e_offset + m_e_min << " " << cumulative_gamma;
			m_printer.debug(6, oss_step_i.str());
		}
		m_printer.debug(6, std::string(60, '-'));


		std::vector<number_t> e_vec_0(n_max);
		std::vector<number_t> gamma_vec_0(n_max);

		if (m_use_heller_derivatives) {

			// calculate gamma values using Heller derivatives as defined by eq. 40 in
			// {reinhardt:1979}

			std::vector<number_t> interpolation_points(n_max);
			number_t cur_x = number_t(0.);
			for (size_t i = 0; i < n_max; ++i) {

        cur_x += number_t(1.);
				interpolation_points[i] = cur_x;
			}

			std::vector<number_t> gamma_vec_eq = interpolant<number_t>(
				interpolation_points, e_vec, m_printer).get_c2_coefficients();

			for (size_t i = 0; i < n_max; ++i) {

				gamma_vec_0[i] = gamma_vec[i] / gamma_vec_eq[i];
				e_vec_0[i] = e_vec[i] - m_e_offset + m_e_min;

				std::ostringstream oss;
				oss << "e_0[" << i << "] = " << e_vec_0[i] << ", gamma_0["
					<< i << "] = " << gamma_vec_0[i];
				m_printer.debug(6, oss.str());
			}
		} else {

			// calculate gamma and energy values at each gamma_vec[i], gamma_vec[i+1]
			// and e_vec[i], e_vec[i+1] pair, respectively, by simple numerical
			// differentiation according to eq. 3.5.3

			for (size_t i = 0; i != n_max-1; ++i) {

				e_vec_0[i] = number_t(.5) * (e_vec[i+1] + e_vec[i]);
				// care about the possibly previously introduced energy shift by
				// m_e_offset
				e_vec_0[i] = e_vec_0[i] - m_e_offset + m_e_min;

				gamma_vec_0[i] = number_t(.5)
					* (gamma_vec[i+1] + gamma_vec[i])
					/ (e_vec[i+1] - e_vec[i]);

				std::ostringstream oss;
				oss << "e_0[" << i << "] = " << e_vec_0[i] << ", gamma_0["
					<< i << "] = " << gamma_vec_0[i];
				m_printer.debug(6, oss.str());
			}

      e_vec_0.pop_back();
			gamma_vec_0.pop_back();
		}

		// check the interpolated energies to be compliant with m_e_req, the
		// requested energy to compute gamma for, and put an appropriate value
		// into gamma_histogram
		if (m_e_req < e_vec_0.front() || m_e_req > e_vec_0.back()) {
			if (m_e_req < e_vec_0.front()) {
				std::ostringstream oss;
				oss << "Warning: the requested energy lies below the lowest grid point: "
					<< m_e_req << " < " << e_vec_0.front();
				m_printer.warning(oss.str());
				m_printer.warning("Large inaccuracy expected in the calculation of gamma!");
				// put the appripriate value into gamma_histogram (cf. 3.5.3)
				gamma_histogram[n_max] = number_t(.5) * gamma_vec.front() / e_vec.front();
				continue;
			} else {
				std::ostringstream oss;
				oss << "Warning: the requested energy lies above the largest grid point: "
					<< m_e_req << " > " << e_vec_0.back();
				m_printer.warning(oss.str());
				// put the appropriate value into gamma_histogram (cf. 3.5.3)
				gamma_histogram[n_max] = 0.;
			}
		}

		// interpolate the result using the monotonicity-preserving piecewise
		// cubic Hermite interpolant
		bool is_extrapolated = false;
		interpolant<number_t> interp(e_vec_0, gamma_vec_0, m_printer);
		gamma_histogram[n_max] = interp.interpolate(m_e_req, is_extrapolated);

		bool plot_width_function = false;
		if (plot_width_function) {

			// print out the interpolated width function
			std::ostringstream oss_gamma_plot;
			oss_gamma_plot << "interpolation of Gamma(E) for order " << n_max;
			m_printer.announce(oss_gamma_plot.str(), '-', 6);

			number_t cur_plot_e = e_vec_0.front();
			number_t plot_step_width = number_t(0.002);
			while (cur_plot_e < e_vec_0.back()) {

				number_t cur_plot_gamma = interp.interpolate(cur_plot_e, is_extrapolated);
				std::ostringstream oss_gamma_plot_line;
				oss_gamma_plot_line << cur_plot_e << " " << cur_plot_gamma;
				m_printer.debug(6, oss_gamma_plot_line.str());

				cur_plot_e += plot_step_width;
			}
			m_printer.debug(6, std::string(60, '-'));
		}

		// if we had to extrapolate, this is indicated in <is_extrapolated>
		std::ostringstream oss;
		if (is_extrapolated) {
			oss << "Extrapolated ";
		} else {
			oss << "Interpolated ";
		}
		oss << "gamma for order " << n_max << ": " << gamma_histogram[n_max];
		m_printer.debug(2, oss.str());
		
	}

	// we are done here looping through the usable polynomial orders, and gamma_histogram
	// now holds values at indices between and including usable_poly_orders.first
	// and usable_poly_orders.second
	return;
}


template<class number_t>
void stieltjes<number_t>::sort_ql_energies_and_gammas(
		std::vector<number_t> &e_vec,
		std::vector<number_t> &gamma_vec) {

	// sorts energy and, according to the latter, gamma vectors so that the energy
	// values in e_vec are in an ascending order

	for (size_t i = 0; i != e_vec.size(); ++i) {

		size_t k = i;
		number_t e = e_vec[k];
		for (size_t j = i+1; j != e_vec.size(); ++j) {

			if (e_vec[j] < e) {
				k = j;
				e = e_vec[j];
			}
		}

		if (k != i) {
			std::swap<number_t>(e_vec[i], e_vec[k]);
			std::swap<number_t>(gamma_vec[i], gamma_vec[k]);
		}
	}

	for (size_t i = 0; i != e_vec.size(); ++i) {
		std::ostringstream oss;
		oss << "Sorted QL values for order=" << i << ": e = " << e_vec[i]
			<< ", gamma = " << gamma_vec[i];
		m_printer.debug(6, oss.str());
	}
}


template<class number_t>
void stieltjes<number_t>::get_polynomial_orders(
		std::vector< std::vector<number_t> > &q_polynomials,
		std::pair<size_t, size_t> &usable_poly_orders) {

	// determines the minimal and maximal polynomial orders to be used by
	// checking the orthogonality/overlap between each pair of neighbouring
	// polynomials;  the determined orders are then returned in usable_poly_orders
	//
	// Q polynomials were previously built up from order 0 up to and including
	// m_max_poly_order;  we let n run from 0 to m_max_poly_order-1 and check
	// orthogonality between qpol[n+1] and qpol[n]
	m_printer.info("Checking orthogonality of Q polynomials.");

	size_t max_usable_poly_order = m_max_poly_order - 1;

	for (size_t n = 0; n != m_max_poly_order; ++n) {

		number_t poly_norm(0.), overlap(0.);
		for (size_t i = 0; i != m_num_points; ++i) {

			poly_norm += q_polynomials[n+1][i] * q_polynomials[n+1][i] * m_gammas[i];
			overlap += q_polynomials[n+1][i] * q_polynomials[n][i] * m_gammas[i];
		}

		{
			std::ostringstream oss;
			oss << "Order n=" << n << ", norm=" << poly_norm << ", overlap=" << overlap;
			m_printer.debug(5, oss.str());
		}

		// set the overlap to something small > 0 to avoid zero division-related
		// problems
		//overlap = std::max<number_t>(boost::multiprecision::fabs(overlap), number_t(1.0e-50));
		if (poly_norm / boost::multiprecision::fabs(overlap) < m_max_poly_overlap) {
			// the polynomial failing the orthogonality check (order n+1) should not
			// be used, so we take n here
			max_usable_poly_order = n;
			break;
		}
	}

	size_t min_usable_poly_order;

	if (max_usable_poly_order < 5) {
		min_usable_poly_order = max_usable_poly_order;
		m_printer.warning("Only very low-order approximation is available!");
	} else {
		min_usable_poly_order = 5;
	}

	usable_poly_orders.first = min_usable_poly_order;
	usable_poly_orders.second = max_usable_poly_order;

	// indicate that the approximated gamma should checked for convergence later
	// on if the highest usable order is sufficiently large and enough orders are
	// usable
	if (max_usable_poly_order > 6
			&& (max_usable_poly_order - min_usable_poly_order) >= 2) {
		m_do_convergence_check = true;
	}

	std::ostringstream oss;
	oss << "Using polynomial orders from " << min_usable_poly_order
		<< " up to " << max_usable_poly_order;
	m_printer.info(oss.str());
}


template<class number_t>
void stieltjes<number_t>::generate_polynomials(
		std::vector<number_t> &a_coeffs,
		std::vector<number_t> &b_coeffs,
		std::vector< std::vector<number_t> > &q_polynomials) {

	// initialize the a, b coefficients and Q polynomials for n=0 and n=1
	init_polynomials(a_coeffs, b_coeffs, q_polynomials);

	// build up the higher orders of the Q polynomials using the recurrence
	// relation as defined in {mueller-plathe:1990}, eq. 3.3.22
	number_t a_coeff_sum = a_coeffs[0]; // i. e. a_coeff_sum = 0.
	number_t b_coeff_prod = b_coeffs[0];
	/*
	// set up a vector holding the (1/{\omega_i})^n powers
	std::vector<number_t> inverse_e_powers(m_num_points, 1.);
	for (size_t i = 0; i != m_num_points; ++i) {
		inverse_e_powers[i] /= m_e_diffs[i];
	}
	*/

	for (size_t n = 2; n <= m_max_poly_order; ++n) {
		// n: order of generated polynomial
		// b_coeff_prod holds the product from orders 0 up to n-2
		// a_coeff_sum holds the sum from orders 0 up to n-2
		// inverse_e_powers holds 1/{\omega_i} values to the power of n-1

		// compute b_{n-1} (cf. eq. 3.3.21)
		for (size_t i = 0; i != m_num_points; ++i) {
			//b_coeffs[n-1] += inverse_e_powers[i] * q_polynomials[n-1][i] * m_gammas[i];
			b_coeffs[n-1] += q_polynomials[n-1][i] * m_gammas[i]
				/ boost::multiprecision::pow(m_e_diffs[i], n-1);
		}
		b_coeffs[n-1] /= b_coeff_prod;

		// compute a_n (cf. eq. 3.3.20)
		// update the sum over a coefficients according to the current order
		a_coeff_sum += a_coeffs[n-1];
		// update the product over b coefficients according to the current order
		b_coeff_prod *= b_coeffs[n-1];
		// update the inverse energies to the power n and calculate a_n
		for (size_t i = 0; i != m_num_points; ++i) {
			//inverse_e_powers[i] /= m_e_diffs[i];
			//a_coeffs[n] += inverse_e_powers[i] * q_polynomials[n-1][i]  * m_gammas[i];
			a_coeffs[n] += q_polynomials[n-1][i]  * m_gammas[i]
				/ boost::multiprecision::pow(m_e_diffs[i], n);
		}
		a_coeffs[n] /= b_coeff_prod;
		a_coeffs[n] -= a_coeff_sum;

		// now we have everything to compute q_polynomials[n] (cf. eq. 3.3.22)
		for (size_t i = 0; i != m_num_points; ++i) {
			q_polynomials[n][i] = (number_t(1.) / m_e_diffs[i] - a_coeffs[n]) *
				q_polynomials[n-1][i];
			q_polynomials[n][i] -= b_coeffs[n-1] * q_polynomials[n-2][i];
		}

		{
			std::ostringstream oss;
			oss << "Computed a coefficient for n=" << n << ": " << a_coeffs[n]
				<< ", b coefficient for n-1: " << b_coeffs[n-1];
			m_printer.debug(5, oss.str());
		}
	} // polynomial order n
}


template<class number_t>
void stieltjes<number_t>::init_polynomials(
		std::vector<number_t> &a_coeffs,
		std::vector<number_t> &b_coeffs,
		std::vector< std::vector<number_t> > &q_polynomials) {
	// generates the initial Q polynomials, i. e. n=0 and n=1, as given in
	// {mueller-plathe:1990}, eq. 3.3.23

	// a coefficient starts with n=1, b with n-1=0
	for (size_t i = 0; i != m_num_points; ++i) {
		a_coeffs[1] += m_gammas[i] / m_e_diffs[i]; // as in 3.3.20
		b_coeffs[0] += m_gammas[i]; // this line is not clear from eq. 3.3.21 since
		                            // the resulting sum should be multiplied by the
		                            // product of b coefficients b_0 ... b_{n-1}
		                            // which is nowhere defined
	}
	a_coeffs[1] /= b_coeffs[0]; // results with 3.3.20 if one takes b_0 as the sum
	                            // over all oscillator strengths as defined above

	// now set up the q polynomials for n=0 (cf. 3.3.23,
	//
	// \f$ Q_0\left( \frac{1}{\omega_i} \right) = 1 \f$
	//
	// and n=1
	//
	// \f$ Q_1\left( \frac{1}{\omega_i} \right) = \frac{1}{\omega_i} - a_1 \f$
	for (size_t i = 0; i != m_num_points; ++i) {
		q_polynomials[0][i] = 1.;
		q_polynomials[1][i] = 1. / m_e_diffs[i] - a_coeffs[1];
	}
}


template<class number_t>
void stieltjes<number_t>::shift_energies() {
	// shift the input energy differences to avoid problems arising from small
	// denominators
	for (typename std::vector<number_t>::iterator it = m_e_diffs.begin();
			it != m_e_diffs.end(); ++it) {

		*it += m_e_offset - m_e_min;
	}
}


template<class number_t>
double stieltjes<number_t>::mean_gamma() {
		// The requested energy doesn't lie within the range of passed energies.
		// In this case, the mean value of the passed input gammas is returned.

	std::ostringstream oss1, oss2, oss3;
	oss1 << "Requested energy is outside the borders of passed energies:";
	oss2 << "Requested " << m_e_req << ", not in range ["
		<< m_e_min << ", " << m_e_max << "]";
	oss3 << "The value calculated will be the mean of the input gammas!";
	m_printer.warning(oss1.str());
	m_printer.warning(oss2.str());
	m_printer.warning(oss3.str());

	number_t gsum(0.0);

	for (typename std::vector<number_t>::const_iterator it = m_gammas.begin();
			it != m_gammas.end(); ++it) {

		gsum += *it;
	}

	number_t gmean(gsum / m_num_points);

	return gmean.template convert_to<double>();
}


template<class number_t>
void stieltjes<number_t>::print_vec(
		const std::vector<number_t> &vec,
		std::string pre, size_t n_max, size_t debug_level) {

	// print n_max elements of the vector; if n_max == 0, all elements shall be
	// printed (default argument); if so, set n_max to a proper value, i. e. the
	// size of the vector, first
	if (n_max == 0) {
		n_max = vec.size();
	}

	std::ostringstream oss;
	oss << pre;
	oss << "length " << vec.size() << ", printing " << n_max << " elements: ";

	size_t el_cnt = 0;
	for (typename std::vector<number_t>::const_iterator it = vec.begin();
		it != vec.end() && el_cnt != n_max; ++it, ++el_cnt) {

		if (it != vec.begin())
			oss << ", ";
		oss << *it;
	}

	m_printer.debug(debug_level, oss.str());

}


} // namespace fanoman


#endif // FANOMAN_STIELTJES_STIELTJES_IMPL_H

