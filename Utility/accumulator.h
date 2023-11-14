#include <list>
#include <utility>

template <typename T, typename BinaryOp>
class fixed_cost_accumulator {

public:
    template <bool has_inplace>
    struct combiner {
        void operator()(T& lhs, const T& rhs) { compute_traits::inplace_combine(lhs, rhs); }
    };
    
    template <>
    struct combiner<false> {
        void operator()(T& lhs, const T& rhs) { T tmp; compute_traits::combine(lhs, rhs, tmp); compute_traits::move(tmp, lhs); }
    };
    
    typedef combiner<typename compute_traits::has_inplace_combine> combiner_type;
    
    
    struct cost_value_pair {
        size_t              cost;
        compute_value_type  value;
    };
    
    fixed_cost_accumulator();
    
    accumulate(InIt first, InIt last) {
        for (InIt it = first; it != last; ++it) {
            accumulate(*it);
        }
    }
    
    accumulate(const_reference x);
    result(compute_value_type& r);
    
private:
    typedef std::list<cost_value_pair>  lod_stack_t;
    typedef lod_stack_t::iterator       lod_iterator;
    
    hierarchy_cost_t state_;
};

/*
o input valuetype
o compute/result valuetype
o compute_traits
  * move
  * has_inplace_combine
    - inplace_combine
    - combine
*/

void
fixed_cost_accumulator::accumulate_it(FwdIt it) {
    
    if (!is_top_empty()) {
        combine_impl::combine(top(), *it, state_.front();
        clear_top();
    }
    else {
        // merge 2 elements with same cost
        keep_top_reference(1, it);
    }
}


// ImplTraits<value_type and aggregate_type are different>
combine(const_reference x, const_reference y, aggregate_reference result) {
    compute_traits::assign(x, aggr_temp_1);
    compute_traits::assign(y, aggr_temp_2);
    combine(temp_1, temp_2, result);
}

combine(const_aggregate_reference x, const_reference y, aggregate_reference result) {
    compute_traits::assign(y, aggr_temp_2);
    combine(x, temp_2, result);
}

combine(const_reference x, const_aggregate_reference y, aggregate_reference result) {
    compute_traits::assign(x, aggr_temp_1);
    combine(temp_1, y, result);
}

combine(const_aggregate_reference x, const_aggregate_reference y, aggregate_reference result) {
    compute_traits::combine(x, y, result);
}

// ImplTraits<value_type and aggregate_type are same>
combine(const_reference, const_reference, reference) {}

void
fixed_cost_accumulator::accumulate_value(const_reference x) {
    
    if (!is_top_empty()) {
        combine_impl::combine(top(), x, state_.front());
        clear_top();
    }
    else {
        // merge 2 elements with same cost
        keep_top_reference(1, x);
    }
}

keep_top_reference(cost_type cost, FwdIt it) {
    top_.cost = cost;
    top_.iterator = it;
    top_.type = ITERATOR;
}

keep_top_reference(cost_type cost, const_reference x) {
    top_.cost = cost;
    top_.value = x;
    top_.type = VALUE;
}


const_reference 
top() const {
    assert(top_.type != EMPTY);
    return top_.type == VALUE ? top_.value : *top_.iterator;
}

bool
is_top_empty() const {
    return top_.type == EMPTY;
}


void
release(
void
fixed_cost_accumulator::accumulate(const_reference x) {
    
    if (!is_top_empty()) {
        combine_impl::combine(top(), *it, state_.front());
        clear_top();
    }
    else {
        lod_iterator it      = state_.begin();
        lod_iterator prev_it = it;
        //compute_value_type result = rit->value; 
        //compute_traits::move(rit->value, result); 
        
        for (++it; it != state_.end(); prev_it = it, ++it) {
            it (it->cost == prev_it->cost) {
                //T result = combine(prev_rit->second, rit->second);
                compute_traits::inplace_combine(it->value, prev_it->value);
                it->cost += prev_it->cost;
                state_.erase(prev_it);   // release(prev_it)
                break;
            }
        }
        state_.push_front(cost_value_pair(1, x);
    }
}


const aggregate_type&
fixed_cost_accumulator::result() {

    lod_iterator it = state_.begin();
    result = it->value; 
   
    for (++rit; rit != state_.rend(); ++rit) {
        result = combine(rit->value, result);     
    }
    

    if (!is_top_empty()) {
        lod_iterator result = allocate_item();
        move_free_to_use(result);
        compute_traits::assign(top(), result.value);
        clear_top();
    }
    lod_iterator it      = state_.begin();
    lod_iterator prev_it = it;
    //compute_value_type result = rit->value; 
    //compute_traits::move(rit->value, result); 
    
    for (++it; it != state.end(); prev_it = it, ++it) {
        //result = combine(rit->value, result);
        compute_traits::inplace_combine(it->value, prev_it->value);
        it->cost += prev_it->cost;
        // or
        // compute_traits::combine(it->value, prev_it->value, result);
        // compute_traits::move(result, it->value);
        // release(prev_it);
        
    }
    
    aggregate_type tmp[2];
    current = 0;
    tmp[current]  = state_.front().value;
    cost = tmp[current].cost;
    release(state_.begin());
    revert_current = !state_.empty();
    for (it = state_.begin(); it != state_.end(); ++it) {
        int other = 1 - current;
        compute_traits::clear(tmp[other]);
        compute_traits::combine(it->value, tmp[current].value, tmp[other]);
        release(it);
        cost += it->cost;
        current = other;
    }
    // current now refers to the "other" element if the for loop has exec-ed.
    if (revert_current) {
        current = 1 - current;
    }
    assert(state_.empty());
    allocate();    
    compute_traits::assign(tmp[current], state_.front().value);
    //r = result;
    return it->value;
}


void
fixed_cost_accumulate(
    InIt        first,
    InIt        last,
    T           init,
    BinaryOp    combine
) {
    typedef std::pair<size_t, T> cost_value_pair;
    std::list<cost_value_pair> state;
    typedef std::list<cost_value_pair>::iterator forward_it;
    typedef std::list<cost_value_pair>::reverse_iterator reverse_it;
    
    state.push_back(std::make_pair(1, init);
    for (InIt it = first; it != last; ++it) {
        state.push_back(std::make_pair(1, *it));
        reverse_it rit = state.rbegin();
        reverse_rit prev_rit = rit;
        ++rit;
        for (; rit != state.rend(); ++rit) {
            it (rit->first == prev_rit->first) {
                T result = combine(prev_rit->second, rit->second);
                size_t cost = rit->first + prev_rit->first;
                *rit = std::make_pair(cost, result);
                state.erase(prev_rit);
                break;
            }
        }
    }

    reverse_it rit = state.rbegin();
    T result = *rit;
    ++rit;
    
    for (; rit != state.rend(); ++rit) {
        result = combine(rit->second, result);
    }
    
    return result;
    
}

void variable_cost_accumulator();
void unordered_var_cost_accumulator();
