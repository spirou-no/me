#pragma once

#include <vector>
#include <iostream>

namespace oyrke { namespace algorithm { namespace fixture {

    typedef ptrdiff_t day_index_t;
	typedef ptrdiff_t team_index_t;
    typedef std::vector<team_index_t> team_index_list_t;
    typedef std::vector<day_index_t>  day_index_list_t;
    static const int UNASSIGNED = -1;
    static const int UNUSED = -2;
    namespace detail {
        class team {
        public:
            team();
            team(team_index_t my_ix, size_t day_count, size_t repetitions);

            void insert_match(day_index_t d, team_index_t t);
            void remove_match(day_index_t t);

            bool         is_matched(day_index_t day_ix) const;
            day_index_t  matched_day(team_index_t team_ix, size_t repetition_ix) const;
            team_index_t matched_team(day_index_t day_ix) const;
            size_t       match_count(team_index_t other_team) const;

            size_t       teams_not_allocated() const;    // REmOVE???
            size_t       open_days() const;
            day_index_t  get_open_day(size_t use_day) const;
            team_index_t get_open_team(size_t use_team) const;

            bool         has_open_matches() const;
            team_index_t key() const { return key_; }

        private:
            team_index_t      key_;
            team_index_list_t open_teams_;  // v[day] = team, i.e. I meet other team @ given day
            size_t            unallocated_days_; // unallocated days
            size_t            repetitions_;
        };

        class day {
        public:
            day();
            day(day_index_t my_ix, size_t team_count, size_t day_count);
            void insert_match(team_index_t ta, team_index_t tb);
            void remove_match(team_index_t ta, team_index_t tb);
            void remove_team(team_index_t t);

            team_index_t matched_team(team_index_t t) const;
            bool         is_matched(team_index_t t) const;

            size_t       teams_not_allocated() const;

            bool         has_open_matches() const;
            team_index_t get_open_team(size_t use_team) const;

            day_index_t key() const { return key_; }

        private:
            day_index_t       key_;
            team_index_list_t allocated_team_;  // v[team_1] = team_2 team_1 and team_2 meet on my day
            size_t            unallocated_;
        };
    }


    /*
    Statistics
    o number of teams
    o number of days
    o number of heuristic moves
    o time used
    o number of draws from match-to-match and day-to-day
    o change in fixture from one match to the next (and day to day).
    */
    class random_fixture {
    public:
        random_fixture(size_t number_of_teams);

        //! Set number of times the teams meet.  
        void set_team_count(size_t n_teams);
        void set_repetitions(size_t repetitions);
        void make_full_fixture(/* progress_reporter */);
        void draw_next_match();
        void prepare_draw();

        void get_fixture(/* fixture */);

        bool is_complete() const;

        size_t number_of_days() const;
        size_t number_of_matches() const;
        size_t number_of_teams() const;
        size_t number_of_repetitions() const;

        size_t current_drawn_matches() const;

        void print(std::ostream& os) const;

    private:
        typedef std::vector<detail::team>  team_list_t;
        typedef std::vector<detail::day>   day_list_t;

        day_list_t   days_;
        team_list_t  teams_;

        size_t current_match_count_;
        size_t num_teams_;
        size_t num_repetitions_;
        size_t retry_count_;
        // num_days
        // num_matches

        //! Return the number of days t still has to play
        size_t unallocated_days_for_team(const detail::team *t) const;

        //! Return the number of teams t still has to play
        size_t unallocated_teams_for_team(const detail::team *t) const;

        //! Return number of teams who are still unallocated for given round
        size_t unallocated_teams_for_day(const detail::day *r) const;

        //! Remove the match between t1 and t2 from day r.
        void remove_match(detail::day *d, detail::team *home, detail::team *away);

        //! Remove given team from day
        void remove_team(detail::day *d, detail::team *t);

        void insert_match(detail::day *d, detail::team *home, detail::team *away);

        void run_heuristic_a();
        void run_heuristic_b();

        detail::day*  random_open_day();
        detail::team* random_open_team();
    };

}}}
