#ifndef _HEADER_ERPORT_ERROR_H_
#define _HEADER_ERPORT_ERROR_H_

#include "jabTypes.h"
#include <iostream>

constexpr static bool print = true;

/**
* @brief Report error message
* @param message the error message
*/
void reportError(const char* message) noexcept;

template <typename Arg, typename... Args>
void JAB_REPORT_ERROR(Arg&& arg, Args&&... args) noexcept {
	if constexpr (print) {
		std::cout << std::forward<Arg>(arg);
		((std::cout << std::forward<Args>(args)), ...);
	}
}

#endif // !_HEADER_ERPORT_ERROR_H_
