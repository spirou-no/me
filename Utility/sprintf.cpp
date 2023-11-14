#include "sprintf.h"

#include <stdio.h>
#include <stdarg.h>


namespace oyrke { namespace utility {
  std::string
  vsprintf(const char *format, va_list args) {
    int maxlen = 1 + _vscprintf(format, args);
    char *buffer = new char[maxlen];
    ::vsprintf(buffer, format, args);
    std::string result(buffer);
    delete[] buffer;

    return result;
  }

  std::string
  vsnprintf(size_t maxlen, const char *format, va_list args) {
    char *buffer = new char[maxlen];
    ::_vsnprintf(buffer, maxlen, format, args);
    std::string result(buffer);
    delete[] buffer;
    return result;
  }

  std::string
  sprintf(const char *format, ...) {
    va_list args;
    va_start(args, format);
    std::string result = vsprintf(format, args);
    va_end(args);

    return result;
  }

  std::string
  snprintf(size_t maxlen, const char *format, ...) {
    va_list args;
    va_start(args, format);
    std::string result = vsnprintf(maxlen, format, args);
    va_end(args);

    return result;
  }
  
}}
