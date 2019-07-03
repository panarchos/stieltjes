#include <sstream>
#include <iostream>
using std::cout;
using std::endl;
#include <deque>
using std::deque;
#include <string>
using std::string;
#include <vector>
using std::vector;

#include "printer.h"
#include "OptionParser/OptionParser.h"
#include "input/InputReader.h"
#include "stieltjes.h"
#include <boost/multiprecision/float128.hpp>
#include <boost/multiprecision/number.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>

extern "C" {
#include <quadmath.h>
}

using namespace fanoman;


Options_t parseCommandline(int argc, char **argv, deque<string> &arguments) {
	OptionParser parser("%prog [options]");
	parser.add_string_option("inputFilename", "-i", "--input-filename",
		"The input file to read gamma values from.  Defaults to 'gamma.in',",
		"gamma.in", "INPUT_FILENAME");
	parser.add_string_option("outputFilename", "-o", "--output-filename",
		"The output filename.  Existing files will be overwritten.  Defaults to"
		" 'gamma.out'.", "gamma.out", "OUTPUT_FILENAME");
	parser.add_int_option("debugLevel", "-d", "--debug-level", "The debug level."
		" This defaults to 3.", 3, "DEBUG_LEVEL");
	parser.add_int_option("maxPolyOrder", "-m", "--max-poly-order", "The maximal"
		" polynomial order to use for Stieltjes imaging.  Defaults to 50.", 50,
		"MAX_POLY_ORDER");
	parser.add_int_option("cpp_dec_floatLen", "-l", "--cpp-dec-float-len",
		"Use boost::multiprecision::cpp_dec_float<DIGITS> as number type for"
		" the evaluation of gamma.  Try 50, 100, 200, 300, 500.  The default is 0"
		" which means use of boost::multiprecision::float128 as number type.",
		0, "DIGITS");
	parser.add_int_option("numEPoints", "-n", "--num-points", "Use so many input"
		" energy/gamma pairs for Stieltjes imaging.  The default is to use all"
		" available points.", 0, "NUM_POINTS");
	parser.add_float_option("scaleGammaBy", "-x", "--scale-gammas-by",
		"Scale the input gamma values by this factor during the Stieltjes procedure.",
		1.0, "SCALING_FACTOR");
	parser.add_bool_option("ev", "", "--ev", "Prescale input gamma values to be"
		" in units of eV.  This is a shortcut for '-x 27.2113845657'.", false);
	parser.add_float_option("shiftEBy", "-s", "--shift-energies-by",
		"Shift the input energies by this number during the Stieltjes procedure.",
		1.0, "ENERGY_SHIFT");
	parser.add_float_option("reqE", "-e", "--requested-energy",
		"Compute gamma at REQ_ENERGY.  This defaults to 0.0.", 0.0, "REQ_ENERGY");
	parser.add_bool_option("hellerDerivatives", "", "--heller-derivatives",
		"Use Heller derivatives instead of Stieltjes ones.", false);
	Options_t options = parser.parse_args(argc, argv, arguments);

	if (options["ev"].get_bool()) {
		options["scaleGammaBy"].set_double_value(27.2113845657);
	}

	return options;
}


void read_input(Options_t &options,
		std::vector<double> &gammas, std::vector<double> &e_diffs) {

	std::string inputFilename = options["inputFilename"].get_string();
	int nPoints = options["numEPoints"].get_int();
	double gScaleFactor = options["scaleGammaBy"].get_double();
	InputReader(inputFilename, nPoints, gScaleFactor).perform(gammas, e_diffs);
}

int main(int argc, char **argv) {
	deque<string> arguments;
	Options_t options = parseCommandline(argc, argv, arguments);

	std::ofstream output(options["outputFilename"].get_string());

	printer m_printer(options["debugLevel"].get_int(), output);

	std::ostringstream header_oss;
	header_oss << "Computing gamma for energy: " << options["reqE"].get_double();
	m_printer.announce(header_oss.str());
	
	std::vector<double> gammas, e_diffs;
	read_input(options, gammas, e_diffs);

	for (auto g : gammas) {
		std::ostringstream oss;
		oss << "Read gamma value: " << g;
		m_printer.info(oss.str());
	}

	for (auto e : e_diffs) {
		std::ostringstream oss;
		oss << "Read e_diff value: " << e;
		m_printer.info(oss.str());
	}

	const int cpp_dec_floatLen(options["cpp_dec_floatLen"].get_int());

	double gamma;
	
	if (cpp_dec_floatLen > 0) {

		std::ostringstream oss;
		oss << "Using boost::multiprecision::number< boost::multiprecision::cpp_dec_float<"
			<< cpp_dec_floatLen << "> > as number type.";
		m_printer.info(oss.str());

		if (cpp_dec_floatLen == 50) {
			stieltjes< boost::multiprecision::number<
					boost::multiprecision::cpp_dec_float<50> > > sti(
				e_diffs,
				gammas,
				options["reqE"].get_double(),
				options["maxPolyOrder"].get_int(),
				options["shiftEBy"].get_double(),
				options["hellerDerivatives"].get_bool(),
				m_printer);
			gamma = sti.compute();
		}
		if (cpp_dec_floatLen == 100) {
			stieltjes< boost::multiprecision::number<
					boost::multiprecision::cpp_dec_float<100> > > sti(
				e_diffs,
				gammas,
				options["reqE"].get_double(),
				options["maxPolyOrder"].get_int(),
				options["shiftEBy"].get_double(),
				options["hellerDerivatives"].get_bool(),
				m_printer);
			gamma = sti.compute();
		}
		if (cpp_dec_floatLen == 200) {
			stieltjes< boost::multiprecision::number<
					boost::multiprecision::cpp_dec_float<200> > > sti(
				e_diffs,
				gammas,
				options["reqE"].get_double(),
				options["maxPolyOrder"].get_int(),
				options["shiftEBy"].get_double(),
				options["hellerDerivatives"].get_bool(),
				m_printer);
			gamma = sti.compute();
		}
		if (cpp_dec_floatLen == 300) {
			stieltjes< boost::multiprecision::number<
					boost::multiprecision::cpp_dec_float<300> > > sti(
				e_diffs,
				gammas,
				options["reqE"].get_double(),
				options["maxPolyOrder"].get_int(),
				options["shiftEBy"].get_double(),
				options["hellerDerivatives"].get_bool(),
				m_printer);
			gamma = sti.compute();
		}
		if (cpp_dec_floatLen == 500) {
			stieltjes< boost::multiprecision::number<
					boost::multiprecision::cpp_dec_float<500> > > sti(
				e_diffs,
				gammas,
				options["reqE"].get_double(),
				options["maxPolyOrder"].get_int(),
				options["shiftEBy"].get_double(),
				options["hellerDerivatives"].get_bool(),
				m_printer);
			gamma = sti.compute();
		}

	} else {

		std::ostringstream oss;
		oss << "Using boost::multiprecision::float128 as number type.";
		m_printer.info(oss.str());

		stieltjes<boost::multiprecision::float128> sti(
			e_diffs,
			gammas,
			options["reqE"].get_double(),
			options["maxPolyOrder"].get_int(),
			options["shiftEBy"].get_double(),
			options["hellerDerivatives"].get_bool(),
			m_printer);
		gamma = sti.compute();
	}

	std::cout << "The final gamma value: " << gamma << std::endl;
	
	return 0;
}

