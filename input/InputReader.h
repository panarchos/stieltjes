#ifndef STIELTJES_INPUT_INPUTREADER_H
#define STIELTJES_INPUT_INPUTREADER_H


#include <fstream>
#include <string>
#include <vector>


namespace fanoman {

class InputReader {
private:
	std::string m_filename;
	int m_n_points;
	double m_g_scale_factor;

public:
	InputReader(
		const std::string &filename_,
		int n_points_,
		double g_scale_factor_) :
		m_filename(filename_),
		m_n_points(n_points_),
		m_g_scale_factor(g_scale_factor_)
	{ }

	void perform(std::vector<double> &gammas, std::vector<double> &e_diffs);
};


void InputReader::perform(
		std::vector<double> &gammas, std::vector<double> &e_diffs) {

	// reads the input from m_filename
	double gamma, e_diff;
	std::ifstream input(m_filename);
	std::vector<double> gammas_, e_diffs_;

	while (input >> e_diff, input >> gamma) {
		e_diffs_.push_back(e_diff);
		gammas_.push_back(m_g_scale_factor * gamma);
	}

	int n = m_n_points;
	if (n == 0 || n > gammas_.size()) {
		n = gammas_.size();
	}

	for (size_t i = 0; i != n; ++i) {
		gammas.push_back(gammas_[i]);
		e_diffs.push_back(e_diffs_[i]);
	}
}


} // namespace fanoman


#endif // STIELTJES_INPUT_INPUTREADER_H

