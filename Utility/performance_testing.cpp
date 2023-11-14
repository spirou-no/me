#include <Utility/performance_testing.h>
#include <Utility/algobase.h>
#include <Utility/stopwatch.h>

#include <algorithm>
#include <limits>

namespace oyrke { namespace test {

    //! Run a test a fixed number of iterations, return elapsed time in seconds.
    double 
    run_performance_test(
        int             n_loops,    //!< Number of iterations to run f
        testfunction_fn f           //!< Function to be evaluated pr iteration
    ) {
        oyrke::utility::stopwatch timer;
        int result = 0;
        for (int i = 0; i < n_loops; ++i) {
            result += f();
        }
        double elapsed = timer.elapsed_time();
        return elapsed;
    }



    //! Run testfunction for a given period of time.  
    // Return actual number of iterations and elapsed time.
    void
    run_adaptive_performance_test(
        double          max_time,       //!< How long to run the tests (in seconds)
        testfunction_fn f,              //!< Function to be evaluated
        int            &n_iterations,   //!< Out, actual number of iterations test was called
        double         &test_time       //!< Out, actual elapsed test time
    ) {
        n_iterations = 1;
        bool done    = false;
        bool has_overflowed = false;
        do {
            test_time = run_performance_test(n_iterations, f);
            double ratio = max_time / test_time;
            done = has_overflowed || ratio < 1.1;
            if (!done) {
                int test_time_ms = int(1000.0 * test_time);
                test_time_ms = oyrke::algorithm::clip(test_time_ms, 10, 1000);
                ratio = std::min(ratio, 10.0);
                int new_iters = int(n_iterations * ratio);
                if (new_iters < n_iterations) {
                    n_iterations = std::numeric_limits<int>::max();
                    has_overflowed = true;
                }
                else {
                    new_iters = std::max(n_iterations+1, new_iters);
                    n_iterations = new_iters;
                }
            }
        } while (!done);
    }

    void
    run_adaptive_performance_test(const char *title, double max_time, testfunction_fn f) {
        int iters = 0;
        double test_time = 0.0;

        run_adaptive_performance_test(max_time, f, iters, test_time);
        double time_pr_loop = test_time / double(iters);
        ::printf("%s %12.3f ms  %10d  %10.2f ns\n", title, 1e3*test_time, iters, 1e9 * time_pr_loop);
    }

    void
    run_adaptive_performance_test(const std::string& title, double max_time, testfunction_fn f) {
        run_adaptive_performance_test(title.c_str(), max_time, f);
    }

}}