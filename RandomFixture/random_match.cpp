#include <RandomFixture/random_match.h>
#include <Utility/stopwatch.h>
#include <utility/algorithm.h>

#include <algorithm>
#include <functional>
#include <vector>
//#include <stdio.h>
#include <iostream>
#include <math.h>
#include <assert.h>

//#define heavy_assert assert
#define heavy_assert

/*
Here is the algorithm for generating fixtures at random. As you might
guess, it is a bit more complicated than just choosing matches at
random and hoping they fit together.

The algorithm has two parts, or "heuristics". You run the algorithm by
choosing one of the two heuristics at random (with equal probability),
and then repeating until all the fixtures have been found. (I will do
an example below.) All choices in this algorithm are made at random,
with all the possible choices having equal probability. The algorithm
works for any even number n of teams -- if you have an odd number of
teams, add a "dummy" team to make n even. (We talked about this before.)

H1:
------------
1. Choose a team, x, that has not played all its matches.
2. Choose a round, r, in which team x does not play.
3. Choose another team, y, that team x does not play.
4 (a) If team y does not play in round r, add the match "x vs y" to round r.
  (b) Otherwise, remove the match that team y previously played in round r, 
      and replace it with "x vs y".

H2:
------------
1. Choose a round, r, that has less than n/2 matches in it.
2. Choose two teams, x and y, that do not play in round r.
3. (a) If "x vs y" does not exist in another round, add "x vs y" to round r.
   (b) Otherwise, add "x vs y" to round r, and remove it from the round
       in which it previously existed.

Here is an example with n=4 teams. It takes 10 steps to arrange the
6 matches:

Step 1: Use H1. x=3, r=1, y=4. Add "3 vs 4" to round 1. (This was a 
        completely free choice).
Step 2: Use H1. x=2, r=1, y=4. (Another free choice.) Add "2 vs 4" to 
        round 1, and remove "3 vs 4".
Step 3: Use H2. r=3, x=1, y=4. (Another free choice, because we don't
        have any matches in round 3 yet.) Add "1 vs 4" to round 3.
Step 4: Use H2. r=2, x=1, y=3. Add "1 vs 3" to round 2.
(There seems to be a problem here. 2 vs 4 and 1 vs 3 must be in the
same round, but the algorithm puts them in different rounds.)
Step 5: Use H1. x=3, r=3 (r could not be 2 because team 3 already plays
in round 2), y=4 (could not be 1). Add "3 vs 4" to round 3, and remove
"1 vs 4".
Step 6: Use H1. x=1, r=3 (cannot be 2), y=4 (cannot be 3). Add "1 vs 4"
to round 3, and remove "3 vs 4".
(This completely undoes step 5! In theory, the algorithm might go on
doing this for ever, but, of course, in practice it will eventually
do something else.)
Step 7: Use H2. r=3, so that x=2 and y=3 (or x=3 and y=2) is the only
possible choice. Add "2 vs 3" to round 3.
Step 8: Use H1. x=2, so r=2 and y=1 is the only possible choice. Add
"2 vs 1" to round 2, and remove "1 vs 3". This fixes the problem caused
in Step 4.
Step 9: Use H2. r=2, so x=3 and y=4 (or x=4, y=3). Add "3 vs 4" to 
round 2.
Step 10: Use H1. x=1 (x had to be either 1 or 3), so r=1 and y=3. 
Add "1 vs 3" to round 1, and we are done.

Round 1: 2 vs 4, 1 vs 3.
Round 2: 2 vs 1, 3 vs 4.
Round 3: 1 vs 4, 2 vs 3.

If you check the heuristics, you will see that it is impossible to have
two teams meet twice, and it is also impossible to remove a match
without adding another. In practice, the algorithm always seems to 
reach the end without getting stuck. It's also easy to program; I found
it was useful to have an array a[i,j] where a[i,j]=0 if "i vs j" has
not been added yet, and a[i,j]=k if "i vs j" is currently in round k.
In this way, it is easy to check whether a match between any two teams
has been arranged yet, and if so, for which round.

Of course, once we have these "random" fixtures, we still have to arrange
who plays at home in each round, and we also have to make sure that
two teams from the same city do not play at home in the same round.

This probably does not help with the English league, though, because that
started in 1888 and the paper where I read about this algorithm was
published in 1987! (The reference, if you want to read about it yourselves,
is "A hill-climbing algorithm for the construction of one-factorizations
and Room squares", by J.H.Dinitz and D.R.Stinson, in the SIAM Journal of
Algebra and Discrete Mathematics (1987, vol. (Band?) 8, pp. (Seiten) 
430-438). Oyvind, was that "Room" squares you meant? I don't know about 
them myself, so maybe you could enlighten us.)

-- Ken Butler butler@cs.sfu.ca 
*/
/*
The "standard" algorithm applies to a single repetition played on neutral ground, or at least
nobody cares about who is playing home and away.

The following describes how to extend the algorithm to deal with
1. Repetitions, i.e. the teams play each other several times, but still on neutral ground
2. Home and away matches
3. Mix all, i.e. a mixture of home, away and neutral matches with repetitions.
   I.e. each team plays M home matches and N on neutral ground.

4. Lock some matches and/or days.  I.e. ensure that the 2 first and 2 last days all teams have
   one home match.
5. Lock m first days.  E.g. with dynamic match allocation, one will attempt to pair the most
   even teams towards the end of the tournament.  The rank changes as matches are played.
   Hence the m first days are fixed when redrawing the remaining days.

The travelling tournament problem and various optimizations
6. Minimize travel distance/cost.  Each team can have k consecutive away matches.

H1:
------------
1. Choose a team, x, that has not played all its matches.
2. Choose a round, r, in which team x does not play.
3. Choose another team, y, that has not played all matches against team x [that team x does not play.]
4 (a) If team y does not play in round r, add the match "x vs y" to round r.
  (b) Otherwise, remove the match that team y previously played in round r, 
      and replace it with "x vs y".

H2:
------------
1. Choose a round, r, that has less than n/2 matches in it.
2. Choose two teams, x and y, that do not play in round r.
3. (a) If x vs y is repeated less than N times, add x vs y to round r
   (b) If x vs y is repeated N times,add x vs y to round r, and
       remove one of the previous matches.
[3. (a) If "x vs y" does not exist in another round, add "x vs y" to round r.
    (b) Otherwise, add "x vs y" to round r, and remove it from the round
        in which it previously existed.]


Solution for 1-3 is quite simple.
a. Draw as for neutral matches, i.e. don't decide home & away yet.
b. Assign home & neutral matches (not in here).

To support multiple repetitions, make changes as follows
day: no changes
team: change open_days_ to open_teams_, i.e. open_teams_[day_ix] = team_ix tells which other
      team I meet on a given day

      === team ===
      insert_match        -- as is
remove_match        -- remove an arbitrary match against given team
                       OR provide day index instead of team index
matched_day         -- --> matched_team(day_ix)???
is_matched          -- return # matches?
teams_not_allocated -- as is, change impl
open_days           -- as is

=== heuristic A ===


=== heuristic B ===

*/


char
convert_team(size_t x) {
    char result = 0;
    if (x < 26) {
        result = char(x+'A');
    }
    else if (x < 52) {
        result = char(x - 26 + 'a');
    }
    else if (x < 62) {
        result = char(x - 52 + '0');
    }
    return result;
}


size_t
random(size_t max_num) {
    /*
    int x = rand();
    double d = double(x) / RAND_MAX;
    d *= max_num;
    int y = int(d);
    assert(y < max_num);
    */
    // rand returns inclusive RAND_MAX
    return int((double(rand()) / (1+RAND_MAX)) * max_num);
}

// TODO: initialize random generator
// better: Use boost random generator.

namespace oyrke { namespace algorithm { namespace fixture { 

    namespace detail {
        bool is_special(size_t x) {
            return x == UNASSIGNED || x == UNUSED;
        }


        size_t
        count_duplicates(const std::vector<int> &to_test) {
            std::vector<int> tmp(to_test);
            std::vector<int>::iterator it_end = std::remove_if(tmp.begin(), tmp.end(), is_special);
            std::sort(tmp.begin(), it_end);

            // unique doesn't remove anything, it simply moves the duplicates to the end.
            size_t unique_items = std::unique(tmp.begin(), it_end) - tmp.begin();
            size_t duplicates = (it_end - tmp.begin()) - unique_items;
            
            return duplicates;
        }

        size_t
        count_value(const std::vector<int> &to_test, int x) {
            size_t c = std::count(to_test.begin(), to_test.end(), x);
            return c;
        }


        template <typename FwdIter, typename T>
        FwdIter
        nth_element_equal(FwdIter first, FwdIter last, T value, size_t n) {
            FwdIter found_iter = last;
            for (; first != last; ++first) {
                if (*first == value) {
                    if (n == 0) {
                        found_iter = first;
                        break;
                    }
                    --n;
                }
            }

            return found_iter;
        }

        template <typename FwdIter>
        FwdIter
        nth_element_pdf(FwdIter first, FwdIter last, size_t n) {
            FwdIter found_iter = last;
            for (; first != last; ++first) {
                size_t bin_pdf = *first;
                if (n < bin_pdf) {
                    found_iter = first;
                    break;
                }
                else {
                    n -= bin_pdf;
                }
            }

            return found_iter;
        }


        template <typename FwdIter, typename T, typename SIZE>
        FwdIter
        find_nth_element(const FwdIter &first, const FwdIter &last, const T &value, SIZE n) {
            FwdIter iter = first;
            for (; iter != last; ++iter) {
                if (*iter == value) {
                    if (n == 0) {
                        break;
                    }
                    else {
                        --n;
                    }
                }
            }
            return iter;
        }

        template <typename FwdIter, typename UnaryPred, typename SIZE>
        FwdIter
        find_nth_element_if(const FwdIter &first, const FwdIter &last, const UnaryPred &pred, SIZE n) {
            FwdIter iter = first;
            for (; iter != last; ++iter) {
                if (pred(*iter)) {
                    if (n == 0) {
                        break;
                    }
                    else {
                        --n;
                    }
                }
            }
            return iter;
        }



        team::team() { }

        team::team(team_index_t my_ix, size_t day_count, size_t repetitions) 
            : key_(my_ix)
            , unallocated_days_(day_count)
            , repetitions_(repetitions) { 
            open_teams_.resize(day_count, (size_t)UNASSIGNED);
            //open_teams_[key_] = UNUSED;
        }
    
        size_t
        team::match_count(team_index_t other_team) const {
            size_t num = std::count(open_teams_.begin(), open_teams_.end(), other_team);
            return num;
        }

        bool
        team::is_matched(day_index_t day_ix) const {
            return !is_special(open_teams_[day_ix]);
        }

        day_index_t
        team::matched_day(team_index_t team_ix, size_t repetition_ix) const {
            team_index_list_t::const_iterator it = find_nth_element(open_teams_.begin(), open_teams_.end(), team_ix, repetition_ix);
            day_index_t day_ix = it == open_teams_.end() ? UNASSIGNED : it - open_teams_.begin();
            return day_ix;
        }

        team_index_t
        team::matched_team(day_index_t day_ix) const {
            return open_teams_[day_ix];
        }

        size_t
        team::teams_not_allocated() const {
            return unallocated_days_; // MULTI_REP
        }

        size_t
        team::open_days() const {
            return unallocated_days_;
        }

        bool
        team::has_open_matches() const {
            return unallocated_days_ != 0;
        }

        void
        team::insert_match(day_index_t d, team_index_t t) {
            assert(open_teams_[d] == UNASSIGNED);
            open_teams_[d] = t;
            --unallocated_days_;
            assert(unallocated_days_ >= 0);
            //heavy_assert(count_duplicates(open_teams_) == 0);
            // heavy_assert(max_duplicates(open_teams_) <= num_repetitions;
            //heavy_assert(count_value(open_teams_, UNASSIGNED) == unallocated_days_);
        }

        void
        team::remove_match(day_index_t d) {
            assert(open_teams_[d] != UNASSIGNED);
            open_teams_[d] = UNASSIGNED;
            ++unallocated_days_;
            assert(unallocated_days_ <= open_teams_.size());
            //heavy_assert(count_duplicates(open_teams_) == 0);
            //heavy_assert(max_duplicates(open_teams_) <= num_repetitions;
            //heavy_assert(count_value(open_teams_, UNASSIGNED) == unallocated_days_);
        }

        /*! Find the team index corresponding to the "open index".
        open index refers to an open match index, and not a team index!
        This number identifies a team, where a team has as many slots as it has unplayed matches.
        */
        team_index_t
        team::get_open_team(size_t open_no) const {
            // return the open_no'th UNASSIGNED index
            //return nth_element_equal(open_days_.begin(), open_days_.end(), UNASSIGNED, open_no) - open_days_.begin();
            // MULTI_REP
			size_t team_count = 1 + open_teams_.size() / repetitions_;
            std::vector<team_index_t> team_alloc(team_count, repetitions_);
            for (size_t i = 0; i < open_teams_.size(); ++i) {
                team_index_t team_ix = open_teams_[i];
                if (team_ix != UNASSIGNED) {
                    --team_alloc[team_ix];
                }
            }

            team_alloc[key_] = 0;
            team_index_list_t::const_iterator it = nth_element_pdf(team_alloc.begin(), team_alloc.end(), open_no);
            team_index_t team_ix = it - team_alloc.begin();

            assert(static_cast<size_t>(team_ix) < team_count);
            return team_ix;
        }


        day_index_t
        team::get_open_day(size_t day_no) const {
            day_index_t day = nth_element_equal(open_teams_.begin(), open_teams_.end(), UNASSIGNED, day_no) - open_teams_.begin();
            assert(day < static_cast<int>(open_teams_.size()));
            return day;
        }




        day::day() { }

        day::day(day_index_t ix, size_t team_count, size_t /*day_count*/) : key_(ix), unallocated_(team_count) {
            allocated_team_.resize(team_count, UNASSIGNED); // day_count
        }

        team_index_t 
        day::matched_team(team_index_t t) const { 
            return allocated_team_[t]; 
        } 

        bool
        day::is_matched(team_index_t t) const {
            return !is_special(allocated_team_[t]);
        }
        
        size_t
        day::teams_not_allocated() const {
            return unallocated_;
        }

        bool
        day::has_open_matches() const {
            return unallocated_ != 0;
        }

        void
        day::insert_match(team_index_t ta, team_index_t tb) {
            assert(allocated_team_[ta] == UNASSIGNED);
            assert(allocated_team_[tb] == UNASSIGNED);
            assert(unallocated_ > 0);

            allocated_team_[ta] = tb;
            allocated_team_[tb] = ta;
            unallocated_ -= 2;
            assert(unallocated_ >= 0);
            //heavy_assert(count_duplicates(allocated_team_) == 0);
            //heavy_assert(count_value(allocated_team_, UNASSIGNED) == unallocated_);
        }

        void
        day::remove_match(team_index_t ta, team_index_t tb) {
            assert(allocated_team_[ta] != UNASSIGNED);
            assert(allocated_team_[tb] != UNASSIGNED);
            assert(unallocated_ < allocated_team_.size());

            allocated_team_[ta] = UNASSIGNED;
            allocated_team_[tb] = UNASSIGNED;
            unallocated_ += 2;
            //heavy_assert(count_duplicates(allocated_team_) == 0);
            //heavy_assert(count_value(allocated_team_, UNASSIGNED) == unallocated_);
        }

        void
        day::remove_team(team_index_t t) {
            if (allocated_team_[t] != UNASSIGNED) {
                allocated_team_[t] = UNASSIGNED;
                ++unallocated_;
            }

            assert(unallocated_ < allocated_team_.size());
            assert(unallocated_ >= 0);
            //heavy_assert(count_duplicates(allocated_team_) == 0);
            //heavy_assert(count_value(allocated_team_, UNASSIGNED) == unallocated_);
        }

        team_index_t
        day::get_open_team(size_t team_no) const {
            team_index_t team = nth_element_equal(allocated_team_.begin(), allocated_team_.end(), UNASSIGNED, team_no) - allocated_team_.begin();
            assert(team < static_cast<int>(allocated_team_.size()));
            return team;
        }
    }


    random_fixture::random_fixture(size_t number_of_teams) 
    : num_teams_(number_of_teams), num_repetitions_(1), current_match_count_(0), retry_count_(0) {
    }

    void
    random_fixture::set_team_count(size_t n_teams) {
        num_teams_ = n_teams;
    }


    //! Set number of times the teams meet.  
    void
    random_fixture::set_repetitions(size_t repetitions) {
        num_repetitions_ = repetitions;
    }

    void
    random_fixture::draw_next_match() {
        if (random(2) == 0) {
            run_heuristic_a();
        }
        else {
            run_heuristic_b();
        }
    }

    void
    random_fixture::prepare_draw() {
        // TODO if num_teams is an odd number, add a dummy team.
        size_t num_days = number_of_days();
        days_.resize(num_days);
        teams_.resize(num_teams_);

        for (size_t i = 0; i < num_days; ++i) {
            days_[i] = detail::day(i, num_teams_, num_days);
        }

        for (size_t i = 0; i < num_teams_; ++i) {
            teams_[i] = detail::team(i, num_days, num_repetitions_);
        }
        current_match_count_ = 0;
        retry_count_ = num_teams_ * num_repetitions_;
    }

    size_t
    random_fixture::number_of_days() const {
        size_t num_days = (num_teams_ - 1) * num_repetitions_;
        return num_days;
    }


    //! Total number in fixture (when it is completed)
    size_t
    random_fixture::number_of_matches() const {
        size_t num_matches = number_of_days() * num_teams_ / 2;
        return num_matches;
    }

    //! Number of teams participating in fixture
    size_t
    random_fixture::number_of_teams() const {
        return num_teams_;
    }

    //! Number of repetitions (on neutral ground).
    size_t
    random_fixture::number_of_repetitions() const {
        return num_repetitions_;
    }

    //! Number of drawn matches (currently).  
    // Useful to monitor progress when drawing one match at a time.
    size_t
    random_fixture::current_drawn_matches() const {
        return current_match_count_;
    }

    //! Return true if all matches in fixture are assigned
    bool
    random_fixture::is_complete() const {
        return current_match_count_ == number_of_matches();
    }

    //! Make full fixture
    void 
    random_fixture::make_full_fixture(/* progress_reporter */) {
        prepare_draw();

        size_t num_matches = number_of_matches();
        while (current_match_count_ < num_matches) {
            draw_next_match();
        }
    }

     
    void 
    random_fixture::get_fixture(/* fixture */) {
        // TODO
    }


    //! Return the number of days t still has to play
    size_t 
    random_fixture::unallocated_days_for_team(const detail::team *t) const {
        return t->teams_not_allocated();
    }

    //! Return the number of teams t still has to play
    size_t 
    random_fixture::unallocated_teams_for_team(const detail::team *t) const {
        return t->teams_not_allocated();
    }

    //! Return number of teams who are still unallocated for given round
    size_t 
    random_fixture::unallocated_teams_for_day(const detail::day *d) const {
        return d->teams_not_allocated();
    }


    void 
    random_fixture::run_heuristic_a() {
        // 1st pick a random team with open matches
        detail::team *team_1 = random_open_team();

        // pick a day team_1 doesn't play yet
        size_t to_play = team_1->open_days();
        size_t use_day = random(to_play);
        day_index_t day_ix = team_1->get_open_day(use_day);  
        detail::day *current_day = &days_[day_ix];
         
        // pick which other team team_1 will meet in this round
        size_t open_ix = random(team_1->teams_not_allocated());
        team_index_t team_ix_2 = team_1->get_open_team(open_ix);
        detail::team *team_2 = &teams_[team_ix_2];

        assert(team_1->key() != team_ix_2);

        // do some cleanup.  team_2 may have been drawn for this day already.  Remove it.
        team_index_t remove_team = current_day->matched_team(team_ix_2);
        if (remove_team != UNASSIGNED) {
            current_day->remove_match(team_ix_2, remove_team);  // removes matches from day involving team_2
            teams_[remove_team].remove_match(day_ix);        // MULTI_REP remove_match(day_ix)
            team_2->remove_match(day_ix);                  // MULTI_REP remove_match(day_ix)
        }
        else {
            ++current_match_count_;
        }
        // ok, now we got which 2 teams will meet in the computed day.
        // update datastructures.
        team_1->insert_match(day_ix, team_ix_2);
        team_2->insert_match(day_ix, team_1->key());
        current_day->insert_match(team_1->key(), team_ix_2);
    }

        
    void 
    random_fixture::run_heuristic_b() {
        // pick a random, open day (i.e. a day where not all teams play yet.
        detail::day *current_day = random_open_day();

        // pick 2 teams for that day.  Pick any 2 teams that doesn't play on this day already.
        size_t open_teams = current_day->teams_not_allocated();
        size_t open_1 = random(open_teams);
        size_t open_2 = random(open_teams-1);
        if (open_2 >= open_1) {
            ++open_2;
        }

        team_index_t team_ix_1 = current_day->get_open_team(open_1);
        team_index_t team_ix_2 = current_day->get_open_team(open_2);

        detail::team *team_1 = &teams_[team_ix_1];
        detail::team *team_2 = &teams_[team_ix_2];

        // do cleanup, team_1 and team_2 may play on another day already.
#if 0
        if (team_1->is_matched(team_ix_2)) {
            day_index_t remove_day = team_1->matched_day(team_ix_2);
            days_[remove_day].remove_match(team_ix_1, team_ix_2);

            team_1->remove_match(team_ix_2);
            team_2->remove_match(team_ix_1);
        }
        else {
            ++current_match_count_;
        }
#endif
        if (team_1->match_count(team_ix_2) == num_repetitions_) {
            size_t repetition_ix = random(num_repetitions_);  // which existing match to remove
            day_index_t remove_day = team_1->matched_day(team_ix_2, repetition_ix);
            days_[remove_day].remove_match(team_ix_1, team_ix_2);

            team_1->remove_match(remove_day); //(team_ix_2);
            team_2->remove_match(remove_day); //(team_ix_1);
        }
        else {
            ++current_match_count_;
        }
        // ok, now we got which 2 teams will meet in the computed day.
        // update datastructures.
        day_index_t day_ix = current_day->key();
        team_1->insert_match(day_ix, team_ix_2);
        team_2->insert_match(day_ix, team_ix_1);
        current_day->insert_match(team_ix_1, team_ix_2);
    }


    // TODO handle many repetitions
    void
    random_fixture::print(std::ostream &os) const {
        int width = 2 + int(::log10(double(1+number_of_days())));
        size_t mid_reps = (number_of_repetitions() + 1) / 2;
        os << "   ";
        for (size_t i = 0; i < teams_.size(); ++i) {
            os.width(width);
            os << convert_team(i);
            for (size_t j = 1; j < mid_reps; ++j) {
                os.width(width);
                os << ' ';
            }
        }
        os << std::endl;

        for (size_t i = 0; i < teams_.size(); ++i) {
            os.width(3);
            os << convert_team(i);
            for (size_t j = 0; j < teams_.size(); ++j) {
                char sep = ' ';
                size_t start_reps = i < j ? 0 : mid_reps;
                size_t end_reps   = i < j ? mid_reps : 2*mid_reps; // 2x in case repetitions is odd
                for (size_t k = start_reps; k < end_reps; ++k) {
                    bool get_it = k < number_of_repetitions() && i < j && k < mid_reps || i > j && k >= mid_reps;
                    day_index_t day = get_it ? teams_[i].matched_day(j, k) : UNASSIGNED;
                    if (day != UNASSIGNED) {
                        os << sep;
                        os.width(width-1);
                        os << day;
                        sep = '-';
                    }
                    else {
                        os.width(width);
                        os << ' ';
                    }
                }
            }
            os << std::endl;
        }
    }

    // TODO merge random_open_days and teams.  move random try to own function.
    // Rationale: all teams and days will have some open days/teams unless for the last
    // few draws.  Hence it is more effective to randomly pick an element, and then check if it
    // ok (as it will be most of the time).
    // only if the random picking fails, we will resort to the method of actally searching all elements
    // using random_1_if.
    // Not sure if we need the retry_count_.  If say size/4 random picks fails, then use the fallback.
    // Speedup for 1000 teams is about 6 times!
    detail::day *
    random_fixture::random_open_day() {
        detail::day *draw = 0;
        size_t try_count = retry_count_ > 0  ?  days_.size() / 4 : 0;

        while (draw == 0 && try_count != 0) {
            day_index_t ix = random(days_.size());
            draw = &days_[ix];
            if (!draw->has_open_matches()) {
                draw = 0;
                --try_count;
            }
        }
        if (draw == 0) {
            retry_count_ -= retry_count_ != 0 ? 1 : 0;
            draw = &(*random_1_if(days_.begin(), days_.end(), std::mem_fun_ref(&detail::day::has_open_matches)));
        }
        return draw;
    }


	//! Find a random team which has not played all its matches yet.
    detail::team * 
    random_fixture::random_open_team() {
        detail::team *draw = 0;
        size_t try_count = retry_count_ > 0  ?  teams_.size() / 4 : 0;

		// try a few rounds of picking a team at random, will succeed most of the time...
        while (draw == 0 && try_count != 0) {
            team_index_t ix = random(teams_.size());
            draw = &teams_[ix];
            if (!draw->has_open_matches()) {
                draw = 0;
                --try_count;
            }
        }

		// ... only when a few rounds remain, we will run out of "free" teams, and 
		// we resort to the slow but safe process of searching through all remaining teams.
        if (draw == 0) {
            retry_count_ -= retry_count_ != 0 ? 1 : 0;
            draw = &(*random_1_if(teams_.begin(), teams_.end(), 
				                  std::mem_fun_ref(&detail::team::has_open_matches)));
        }
        return draw;
    }


    void
    random_fixture::remove_match(detail::day *d, detail::team *home, detail::team *away) {
        team_index_t home_ix = home->key();
        team_index_t away_ix = away->key();

		assert(home->matched_team(d->key()) == away->key());
		assert(d->matched_team(home->key()) == away->key());
		
        d->remove_match(home_ix, away_ix);
        home->remove_match(away_ix);
        away->remove_match(home_ix);
    }



    void
    random_fixture::remove_team(detail::day *d, detail::team *t) {
        team_index_t team_ix = t->key();
        if (d->is_matched(team_ix)) {
            team_index_t other_team = d->matched_team(team_ix);
            remove_match(d, t, &teams_[other_team]);
        }
        // TODO
    }

#if 0
    void 
    random_fixture::insert_match(detail::day *d, detail::team *home, detail::team *away) {
        // if home and away has played all their matches, remove one of them.
        // remove matches played by home and away in d.
        // if home vs away already plays on d, then no-op
        // return how many matches where removed
        if (d->matched_team(home->key() != away->key()) {
            // ok, they're not playing in d 
            int open_count = open_day_count(home, away);
            if (open_count == 0) {
                // remove random repetition
                ++removed_count;
            }
            if (d->is_matched(home->key()) {
                d->remove_match(home->key(), d->matched_team(home->key()));
                ++removed_count;
            }
            if (d->is_matched(away->key()) {
                d->remove_match(away->key(), d->matched_team(away->key()));
                ++removed_count;
            }
        }
        // 
        // TODO
        // int open_day_count(team)
        // int open_day_count(home, away)
        // bool is_playing(team, day)
    }


    int
    random_fixture::open_day_count(const detail::team *team) const {
        // return how many days team is not assigned to yet
        // same as unallocated_days_for_team
    }

    int
    random_fixture::open_day_count(const detail::team *ta, const detail::team *tb) const {
        // return how many days ta vs tb has to meet, i.e. not assigned yet.
        return number_of_repetitions() - ta->match_count(tb->key());
    }

    bool
    random_fixture::is_playing(const detail::day *day, const detail::team *team) const {
        // return true if team is playing on given day
        return day->is_matched(team->key());
    }
#endif
}}}
