#ifndef FANOMAN_STIELTJES_QL_H
#define FANOMAN_STIELTJES_QL_H


#include <boost/math/tools/precision.hpp> // for epsilon
#include <boost/multiprecision/number.hpp> // sqrt, pow, fabs

#include <vector>
#include <algorithm> // for std::swap, until C++11
#include <utility>   // for std::swap, starting from C++11

#include "printer.h"


namespace fanoman {

	/** This is an adapted implementation of the tql2 algorithm as stated in
	 *
	 * Hilary Bowlder, R. S. Martin, C. Reinsch, J. H. Wilkinson, Numerische
	 * Mathematik 11, 1968, 293-306.
	 *
	*/

template<class number_t>
class ql {
private:
	size_t m_max_iter;
	number_t m_macheps;
	printer &m_printer;

public:
	ql(printer &printer_,
		size_t max_iter_ = 30) :
		m_printer(printer_),
		m_max_iter(max_iter_),
		m_macheps(boost::math::tools::epsilon<number_t>()) {
	}

	int diagonalize(
		std::vector<number_t> &d,
		std::vector<number_t> &e,
		std::vector< std::vector<number_t> > &z,
		bool do_sort = false);

private:
	void sort_eigenvalues_eigenvectors(
		std::vector<number_t> &d,
		std::vector< std::vector<number_t> > &z);

};


template<class number_t>
int ql<number_t>::diagonalize(
		std::vector<number_t> &d,
		std::vector<number_t> &e,
		std::vector< std::vector<number_t> > &z,
		bool do_sort) {

	int ierr=0;

	size_t n = d.size();

	if (n == 1) {
		m_printer.warning("Error: Dimension of vectors should be greater than 1.");
		ierr = -1;
		return ierr;
	}

	// shift the values of e to the low end of the vector
	for (size_t i = 1; i != n; ++i) {
		e[i-1] = e[i];
	}
	// set the last element in e to zero
	e.back() = 0.;

	number_t f(0.), tst1(0.);

	// find eigenvalue l
	for (size_t l = 0; l != n; ++l) {

		number_t h = m_macheps * (boost::multiprecision::fabs(d[l])
			+ boost::multiprecision::fabs(e[l]));
		tst1 = std::max<number_t>(tst1, h);

		bool is_converged = false;
		size_t m = l;
		for (; m != n; ++m) {
			// look for small sub-diagonal element

			if (boost::multiprecision::fabs(e[m]) <= tst1 && m == l) {
				// eigenvalue l converged, set is_converged true to indicate that we
				// want to continue in l with the next eigenvalue
				is_converged = true;
			}
		}
		// reset m
		--m;

		if (is_converged) {
			d[l] += f;
			continue; // continue in l loop with the next eigenvalue
		}

		// try to converge eigenvalue/eigenvector l within m_max_iter
		// iterations
		size_t j = 0;
		while(boost::multiprecision::fabs(e[l]) > tst1) {
			// increase j and check if we are still in the allowed number of
			// iterations

			if (j == m_max_iter) {
				// set ierr to the number of the eigenvalue which failed to converge
				// and return
				ierr = l;
				return ierr;
			}
			++j;

			// form shift
			number_t p = (d[l+1] - d[l]) / (number_t(2.) * e[l]);
			number_t r = boost::multiprecision::sqrt(
				boost::multiprecision::pow(p, 2) + number_t(1.));
			if (p < number_t(0.)) {
				h = d[l] - e[l] / (p - r);
			} else {
				h = d[l] - e[l] / (p + r);
			}
			for (size_t i = l; i != n; ++i) {
				d[i] -= h;
			}
			f += h;

			// QL transformation
			p = d[m];
			number_t c(1.);
			number_t s(0.);
			for (size_t ii = m; ii > l; --ii) {

				size_t i = ii - 1; // i runs from m-1 down to and including l in steps of -1
				number_t g = c * e[i];
				h = c * p;

				if (boost::multiprecision::fabs(p) >= boost::multiprecision::fabs(e[i])) {

					c = e[i] / p;
					r = boost::multiprecision::sqrt(boost::multiprecision::pow(c, 2)
						+ number_t(1.));
					e[i+1] = s * p * r;
					s = c / r;
					c = number_t(1.) / r;
				} else {

					c = p / e[i];
					r = boost::multiprecision::sqrt(boost::multiprecision::pow(c, 2)
						+ number_t(1.));
					e[i+1] = s * e[i] * r;
					s = number_t(1.) / r;
					c = c / r;
				}

				p = c * d[i] - s * g;
				d[i+1] = h + s * (c * g + s * d[i]);

				// form vector; since we return a vector of eigenvetors, the index
				// ordering is exchanged with respect to the original algorithm
				for (size_t k = 0; k != n; ++k) {

					h = z[i+1][k];
					z[i+1][k] = s * z[i][k] + c * h;
					z[i][k] = c * z[i][k] - s * h;
				}
			}

			e[l] = s * p;
			d[l] = c * p;
		}

	d[l] += f;
	}

	// order eigenvalues and eigenvectors
	if (do_sort)
		sort_eigenvalues_eigenvectors(d, z);

	return ierr;
}


template<class number_t>
void ql<number_t>::sort_eigenvalues_eigenvectors(
		std::vector<number_t> &d,
		std::vector< std::vector<number_t> > &z) {

	size_t n = d.size();

	for(size_t i = 0; i != n; ++i) {

		size_t k = i;
		number_t p = d[i];
		for (size_t j = i + 1; j != n; ++j) {
			// run through all left indices and check if we find a lower eigenvalue

			if (d[j] < p) {
				k = j;
				p = d[j];
			}
		}

		// k now holds the index corresponding to the smallest eigenvalue left
		// so swap values of d[i] and d[k], vectors z[i] and z[k]
		if (k != i) {
			std::swap<number_t>(d[i], d[k]);
			std::swap< std::vector<number_t> >(z[i], z[k]);
		}

	}
	return;
}


} // namespace fanoman


#endif // FANOMAN_STIELTJES_QL_H

