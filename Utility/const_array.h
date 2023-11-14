#pragma once

#include <limits>
#include <Utility/iterator.h>

#include <assert.h>

namespace oyrke { namespace basic {

// remove max_size & capacity --> remove limits.h
// use reference_wrapper, should value_typee be T or reference_wrapper<T>
// pass more args to random_index_iterator
template <class T>
class const_array {
public:
    typedef const_array                         self_type;

    typedef T                                   value_type;
    typedef size_t                              size_type;
    typedef ptrdiff_t                           difference_type;
    typedef const value_type&                   const_reference;
    typedef const_reference                     reference;
    // need something for pointer...
    typedef void                                pointer;
    typedef random_index_iterator<const_array, value_type, void>  const_iterator;
    typedef const_iterator                      iterator;

    // no need (except type safety) to specialize reverse iterators.
    // a const is a const is a const ...
    typedef const_iterator                      const_reverse_iterator;
    typedef const_iterator                      reverse_iterator;

    //
    // Creators
    //
    const_array(size_type sz = 0)
        : size_(sz), value_(value_type()) { }
    const_array(size_type sz, const value_type& x)
        : size_(sz), value_(x) { }
    const_array(const self_type& rhs)
        : size_(rhs.size_), value_(rhs.value_) { }

    self_type& operator=(const self_type& rhs) {
        size_ = rhs.size_;
        value_ = rhs.value_;
        return *this;
    }

    //
    // iterator access
    //
    iterator begin() const          { return iterator(const_cast<self_type*>(this), 0); }
    iterator end() const            { return iterator(const_cast<self_type*>(this), size_); }
    reverse_iterator rbegin() const { return begin(); }
    reverse_iterator rend() const   { return end(); }

    //
    //
    // queries
    //
    size_type size() const        { return size_; }
    bool empty() const            { return size_ == 0; }
    size_type max_size() const    { return numeric_limits<size_type>::max(); }
    size_type capacity() const    { return max_size(); }


    //
    // element access
    //
    const_reference front() const { assert(size_ > 0); return value_; }
    const_reference back() const  { assert(size_ > 0); return value_; }
    const_reference operator[](size_type n) const {
        assert(n >= 0 && n < size_);
        return value_;
    }


    const_reference reference_at(size_type n) const {
        assert(n >= 0 && n < size_);
        return value_;
    }

    //
    // commands
    //
    void resize(size_type n_size) { size_ = size; }
    void set_value(const_reference x) { value_ = x; }
    const_reference value() const    { return value_; }


    bool operator==(const self_type& rhs) const {
        return size_ == 0 && rhs.size_ == 0
            || size_ == rhs.size_ && value_ == rhs.value_;
    }

    bool operator<(const self_type& rhs) const {
        bool lt = false;

        if (size_ == 0) {
            lt = rhs.size_ > 0;
        }
        else {
            lt = rhs.size_ != 0 ? value_ < rhs.value_ : false;
        }
        return lt;
    }


private:
    size_type  size_;
    value_type value_;
};

}}