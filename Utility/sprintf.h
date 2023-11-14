#include <string>
#include <stdarg.h>

namespace oyrke { namespace utility {
  std::string vsprintf(const char *format, va_list args);
  std::string vsnprintf(size_t maxlen, const char *format, va_list args);
  std::string sprintf(const char *format, ...);
  std::string snprintf(size_t maxlen, const char *format, ...);
}}
  
