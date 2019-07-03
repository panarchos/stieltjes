#ifndef FANOMAN_BASE_PRINTER_H
#define FANOMAN_BASE_PRINTER_H

#include <iostream>
#include <string>

#define FANOMAN_PRINT_LEVEL_RESULT  0
#define FANOMAN_PRINT_LEVEL_ERROR   0
#define FANOMAN_PRINT_LEVEL_WARNING 0
#define FANOMAN_PRINT_LEVEL_INFO    1
#define FANOMAN_PRINT_LEVEL_DEBUG   2

namespace fanoman {


class printer {
private:
	unsigned m_print_level;
	std::ostream &m_out;
	std::string m_header;
public:
	printer(unsigned print_level = 0, std::ostream &out = std::cout,
			std::string header = "=== STIELTJES === ") :
		m_print_level(print_level), m_out(out), m_header(header) { }

	printer(const printer &printer_) :
		m_print_level(printer_.m_print_level),
		m_out(printer_.m_out),
		m_header(printer_.m_header) { }

	std::ostream& debug(unsigned debug_level, const std::string msg);
	std::ostream& debug(const std::string msg);
	std::ostream& info(const std::string msg);
	std::ostream& warning(const std::string msg);
	std::ostream& error(const std::string msg);
	std::ostream& result(const std::string msg);

	std::ostream& announce(const std::string msg, const char sep = '*',
		unsigned msg_level = FANOMAN_PRINT_LEVEL_INFO);

private:
	std::ostream& print(unsigned msg_level, const std::string &msg);
	std::ostream& print_(unsigned msg_level, const std::string &msg);
};


} // namespace fanoman


#endif // FANOMAN_BASE_PRINTER_H

