#include <Utility/stopwatch.h>
#include <windows.h>

namespace oyrke { namespace utility {

stopwatch::tick_t stopwatch::frequency_ = 0;


stopwatch::stopwatch() {
  init_static();
  restart();
}

stopwatch::stopwatch(Halted) {
  init_static();
  clear();
}

stopwatch::~stopwatch() {
}

void
stopwatch::start() {
  is_halted_   = false;
  start_ticks_ = ticks_now();
}

void
stopwatch::stop() {
  if (!is_halted_) {
    accumulated_ticks_ += last_split_ticks();
    start_ticks_        = 0;
    is_halted_          = true;
  }
}

void
stopwatch::clear() {
  is_halted_         = true;
  accumulated_ticks_ = 0;
  start_ticks_       = 0;
}

void
stopwatch::restart() {
  accumulated_ticks_ = 0;
  start();
}

bool
stopwatch::is_running() const {
  return !is_halted_;
}

stopwatch::tick_t
stopwatch::last_split_ticks() const {
  return ticks_now() - start_ticks_;
}

stopwatch::tick_t
stopwatch::elapsed_ticks() const {
  return is_halted_ ? accumulated_ticks_ : accumulated_ticks_ + last_split_ticks();
}

double
stopwatch::elapsed_time() const {
  return double(elapsed_ticks()) / frequency_;
}

void
stopwatch::init_static() {
  if (frequency_ == 0) {
    LARGE_INTEGER tmp;
    QueryPerformanceFrequency(&tmp);
    frequency_ = tmp.QuadPart;
  }
}

stopwatch::tick_t
stopwatch::ticks_pr_second() {
  return frequency_;
}

stopwatch::tick_t
stopwatch::ticks_now() {
  LARGE_INTEGER counter;
  QueryPerformanceCounter(&counter);
  return counter.QuadPart;
}
}}