#include <string>

namespace oyrke { namespace test {

    typedef int (*testfunction_fn)();

    double run_performance_test(int n_loops, testfunction_fn f);
    void   run_adaptive_performance_test(double max_time, testfunction_fn f, int &n_iterations, double &test_time);
    void   run_adaptive_performance_test(const char *title, double max_time, testfunction_fn f);
    void   run_adaptive_performance_test(const std::string& title, double max_time, testfunction_fn f);

}}
