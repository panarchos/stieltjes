#include <sstream>

#include "printer.h"


namespace fanoman {


std::ostream& printer::print_(unsigned msg_level, const std::string &msg) {
	if (m_print_level >= msg_level) {
		m_out << msg << std::endl;
	}
	return m_out;
}


std::ostream& printer::print(unsigned msg_level, const std::string &msg) {
	std::ostringstream oss;
	oss << m_header << msg;
	return print_(msg_level, oss.str());
}


std::ostream& printer::debug(unsigned debug_level, const std::string msg) {
	// general debug information printer; debug_level starts from 0 (basic debug).
	// debug messages are printed, if FANOMAN_DEBUG_LEVEL + debug_level >=
	// m_print_level
	return print(FANOMAN_PRINT_LEVEL_DEBUG + debug_level, msg);
}


std::ostream& printer::debug(const std::string msg) {
	// prints basic debug information
	return debug(0, msg);
}


std::ostream& printer::info(const std::string msg) {
	return print(FANOMAN_PRINT_LEVEL_INFO, msg);
}


std::ostream& printer::warning(const std::string msg) {
	return print(FANOMAN_PRINT_LEVEL_WARNING, msg);
}


std::ostream& printer::error(const std::string msg) {
	return print(FANOMAN_PRINT_LEVEL_ERROR, msg);
}


std::ostream& printer::result(const std::string msg) {
	return print(FANOMAN_PRINT_LEVEL_RESULT, msg);
}


std::ostream& printer::announce(const std::string msg,
		const char sep, unsigned msg_level) {
	std::ostringstream oss;
	const std::string sep_string(80, sep);
	oss << sep_string << "\n";
	oss << m_header << msg << "\n";
	oss << sep_string;
	return print_(msg_level, oss.str());
}


} // namespace fanoman

