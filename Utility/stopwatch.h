#pragma once

namespace oyrke { namespace utility {

class stopwatch {
public:
  enum Halted { HALTED };
  typedef __int64 tick_t;

  stopwatch();
  stopwatch(Halted);
  ~stopwatch();

  void start();
  void stop();
  void restart();
  void clear();

  tick_t  elapsed_ticks() const;
  double  elapsed_time() const;

  static tick_t ticks_now();
  static tick_t ticks_pr_second();

  bool is_running() const;

private:
  static void init_static();
  tick_t last_split_ticks() const;

  static tick_t     frequency_;
  tick_t            start_ticks_;
  tick_t            accumulated_ticks_;
  bool              is_halted_;
};

}}
