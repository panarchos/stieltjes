#ifndef FANOMAN_STIELTJES_STIELTJES_H
#define FANOMAN_STIELTJES_STIELTJES_H

#include "printer.h"
#include <vector>
#include <utility> // for pair


namespace fanoman {


template<class number_t>
class stieltjes {
private:
	std::vector<number_t> m_e_diffs;
	const std::vector<number_t> m_gammas;
	bool m_use_heller_derivatives;
	printer &m_printer;

	const size_t m_min_num_points, m_max_num_points;
	double m_e_req; //!< The requested energy, to get the gamma for; for us: always 0
	number_t m_e_offset; //!< The offset energies are shifted with to avoid
	                     //   problems with small denominators; shouldn't affect
	                     //   the result
	size_t m_max_allowed_poly_order;
	size_t m_max_poly_order;
	number_t m_max_poly_overlap;
	bool m_do_convergence_check;
	size_t m_max_ql_iterations;
	number_t m_conv_threshold_initial;
	number_t m_conv_threshold_progress;
	size_t m_conv_max_iter;
	number_t m_conv_fac;

	size_t m_num_points;
	number_t m_e_min, m_e_max;

public:
	stieltjes(
		const std::vector<double> &e_diffs_,
		const std::vector<double> &gammas_,
		double e_req_,
		size_t max_allowed_poly_order_,
		double e_offset_,
		bool use_heller_derivatives_,
		printer &printer_) :
		m_e_diffs(e_diffs_.begin(), e_diffs_.end()),
		m_gammas(gammas_.begin(), gammas_.end()),
		m_use_heller_derivatives(use_heller_derivatives_),
		m_printer(printer_),
		m_min_num_points(4),
		m_max_num_points(15000),
		m_e_req(e_req_),
		m_e_offset(e_offset_),
		m_max_allowed_poly_order(max_allowed_poly_order_),
		m_max_poly_overlap(100.0),
		m_do_convergence_check(false),
		m_max_ql_iterations(30),
		m_conv_threshold_initial(0.05),
		m_conv_threshold_progress(1.2),
		m_conv_max_iter(10),
		m_conv_fac(1.0) {

		init();
	}

	double compute();

private:
	void init();
	void check_input();
	void init_e_bounds();

	double mean_gamma();

	void shift_energies();
	void generate_polynomials(
		std::vector<number_t> &a_coeffs, std::vector<number_t> &b_coeffs,
		std::vector< std::vector<number_t> > &q_polynomials);
	void init_polynomials(
		std::vector<number_t> &a_coeffs, std::vector<number_t> &b_coeffs,
		std::vector< std::vector<number_t> > &q_polynomials);
	void get_polynomial_orders(
		std::vector< std::vector<number_t> > &q_polynomials,
		std::pair<size_t, size_t> &usable_poly_orders);
	void calc_gamma_histogram(
		std::vector<number_t> &a_coeffs, std::vector<number_t> &b_coeffs,
		std::vector< std::vector<number_t> > &q_polynomials,
		std::pair<size_t, size_t> &usable_poly_orders,
		std::vector<number_t> &gamma_histogram);
	void sort_ql_energies_and_gammas(
		std::vector<number_t> &e_vec,
		std::vector<number_t> &gamma_vec);
	number_t find_converged_gamma(
		std::vector<number_t> &gamma_histogram,
		std::pair<size_t, size_t> &usable_poly_orders);

	void print_vec(const std::vector<number_t> &vec,
		std::string pre = "", size_t n_max = 0, size_t debug_level = 0);

};


} // namespace fanoman


#include "stieltjes_impl.h"


#endif // FANOMAN_STIELTJES_STIELTJES_H

