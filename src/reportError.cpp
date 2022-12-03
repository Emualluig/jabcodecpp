#include "reportError.h"

void reportError(const char* message) noexcept {
	if constexpr (print) {
		std::cout << "JABCode Error: " << message << "\n";
	}
}
