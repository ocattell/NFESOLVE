#ifndef DEBUGHEADERDEF
#define DEBUGHEADERDEF

#include <iostream>
#include <stdexcept>
#include <string>

inline void DebugCheck(const bool state, const std::string errorMsg)
{
	if (state)
	{
		std::cerr << std::endl << "error: " << errorMsg << std::endl;
		throw std::logic_error(errorMsg);
	}
}

#endif