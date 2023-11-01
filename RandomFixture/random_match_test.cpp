#include <RandomFixture/random_match.h>
#include <utility/stopwatch.h>
#include <boost/lambda/bind.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <numeric>

namespace {
    class fixture_stats {

        struct time_count_pair {
            double elapsed_time;
            int    count;

            time_count_pair(double t, int c) : elapsed_time(t), count(c) {}
            static bool compare(const time_count_pair &lhs, const time_count_pair &rhs) {
                return lhs.count < rhs.count;
            }
        };
        
        std::vector<time_count_pair> samples_;
        typedef std::vector<time_count_pair>::size_type size_type;

    public:
        fixture_stats() {}
        void add(double t, int count) {
            samples_.push_back(time_count_pair(t, count));
        }

        int maximum() const;
        int minimum() const;
        //int median() const;

        size_t samples() const { return samples_.size(); }
        double mean_time() const;
        double mean() const;
        double stddev() const { return moment(2); }
        double skew() const  { return moment(3); }
        double kurtosis() const  { return moment(4); }

        double moment(int m) const;
    };
    
    double
    fixture_stats::mean_time() const {
        double sum = 0.0;

        for (size_t i = 0; i < samples_.size(); ++i) {
            sum += samples_[i].elapsed_time;
        }
        /*
        double acc = std::accumulate(samples_.begin(), samples_.end(), 0.0, 
        boost::lambda::_1 + boost::lambda::bind(&time_count_pair::elapsed_time, boost::lambda::_2));
        */
        sum /= samples_.size();
        return sum;
    }


    double
    fixture_stats::mean() const {
        double sum = 0.0;

        for (size_t i = 0; i < samples_.size(); ++i) {
            sum += samples_[i].count;
        }
        
        sum /= samples_.size();
        return sum;
    }

    int 
    fixture_stats::minimum() const { 
        return std::min_element(samples_.begin(), samples_.end(), time_count_pair::compare)->count;
    }
    
    
    int 
    fixture_stats::maximum() const { 
        return std::max_element(samples_.begin(), samples_.end(), time_count_pair::compare)->count;
    }
/*
    int
    fixture_stats::median() const {
        size_type med_element = samples_.size()/2;
        std::nth_element(samples_.begin(), samples_.begin() + med_element, samples_.end(), time_count_pair::compare);
        return samples_[med_element].count;
    }
*/
    double
    fixture_stats::moment(int /* m */) const {
        //pow( sum (pow(x-mean,m) ), 1/m)
        double average = mean();
        return average;
        //std::accumulate()
    }
}

// local variable is initialized but not referenced
#pragma warning(disable: 4189)

void 
random_match_test(int n, int r, bool print) {
    oyrke::algorithm::fixture::random_fixture fixture(n);

    fixture.set_team_count(n);
    fixture.set_repetitions(r);
    fixture.prepare_draw();
    size_t draw_count = 0;
    std::vector<ptrdiff_t> delta;
    size_t prev_match_count = 0;
    size_t prev_draw_count = 0;
    size_t prev_draw_printed = 0;
    size_t prev_match_printed = 0;

    oyrke::utility::stopwatch timer;

    do {
        fixture.draw_next_match();
        ++draw_count;
        size_t new_match_count = fixture.current_drawn_matches();

        if (new_match_count > prev_match_count) {
            ptrdiff_t delta_match_count = new_match_count - prev_match_count;
            ptrdiff_t delta_draw_count = draw_count - prev_draw_count;
            delta.push_back(delta_draw_count);
            prev_match_count = new_match_count;
            prev_draw_count  = draw_count;
            if (timer.elapsed_time() > 1.0) {
                ptrdiff_t delta_draw_printed = draw_count - prev_draw_printed;
                ptrdiff_t delta_match_printed = new_match_count - prev_match_printed;
                double r1 = double(draw_count) / new_match_count;
                double r2 = double(delta_draw_printed) / delta_match_printed;

                std::cout //<< std::endl << std::endl 
                        << "Draw # = " << draw_count << " / " << new_match_count << " r = " << r1
                        << "  Delta = " << delta_draw_printed << " / " << delta_match_printed  << " r = " << r2
                        << std::endl;
                prev_draw_printed = draw_count;
                prev_match_printed = new_match_count;
                timer.restart();
            }
        }

        //fixture.print(std::cout);
    } while (!fixture.is_complete());

    std::cout << "Ratio is " << draw_count << " / " << fixture.number_of_matches() << " = " << double(draw_count)/fixture.number_of_matches();
    if (print) {
        std::cout << std::endl;
        fixture.print(std::cout);
    }
}

/*
Error when n=6
Fra trekning 10->11 får D 2 kamper.  D møter C etter trekning #10.  Matrisa viser da C-D og D-C i runde 3.
Trekning 11 setter B-D (og D-B), men bare C-D blir fjernet, dvs D-C ligger igjen i runde 3.

Legg til konsistens sjekk
o team
  - valider at array stemmer overens med unallocated_ variabelen.
  - valider at et lag bare er tilordnet en runde
  - valider at er lag møter et annet bare en gang
o day
  - valider at array stemmer overens med unallocated_ variabelen
  - hvert lag er repetert 0 eller 1 gang.
  - at runden er symmetrisk (vil kanskje endres når nøytral implementasjon blir endret).
o for fixture
  - at antall uallokerte lag stemmer overens med total_matches - current_alloced
  - at day og team strukturer er konsistente, dvs den ene kan lages fra den andre
*/