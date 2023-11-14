#pragma once

namespace oyrke {

    template <typename Owner, typename InVT = typename Owner::value_type, typename OutVT = InVT>
    class random_index_iterator {
        typedef random_index_iterator self_type;
    public:
        typedef std::random_access_iterator_tag     iterator_category;
        typedef InVT                                value_type;
        typedef ::ptrdiff_t                         difference_type;
        typedef value_type*                         pointer;
        typedef value_type&                         reference;

        random_index_iterator() : current_(0), owner_(0) {}
        random_index_iterator(Owner *owner, difference_type pos) : owner_(owner), current_(pos) {}
        random_index_iterator(const self_type& rhs) {
            current_ = rhs.current_;
            owner_   = rhs.owner_;
        }

        self_type& operator=(const self_type& rhs) {
            current_ = rhs.current_;
            owner_   = rhs.owner_;
            return *this;
        }

        self_type& operator+=(difference_type rhs) {
            current_ += rhs;
            return *this;
        }

        self_type& operator-=(difference_type rhs) {
            current_ -= rhs;
            return *this;
        }

        
        self_type& operator++() {
            ++current_;
            return *this;
        }

        self_type  operator++(int) {
            self_type tmp(*this);
            ++(*this);
            return tmp;
        }

        self_type& operator--() {
            --current_;
            return *this;
        }


        self_type  operator--(int) {
            self_type tmp(*this);
            --(*this);
            return tmp;
        }

        bool operator==(const self_type& rhs) const {
            return current_ == rhs.current_;
        }

        bool operator!=(const self_type& rhs) const {
            return current_ != rhs.current_;
        }

        bool operator< (const self_type& rhs) const {
            return current_ < rhs.current_;
        }

        bool operator<=(const self_type& rhs) const {
            return current_ <= rhs.current_;
        }

        bool operator> (const self_type& rhs) const {
            return current_ > rhs.current_;
        }

        bool operator>=(const self_type& rhs) const {
            return current_ >= rhs.current_;
        }

        difference_type current() const { return current_; }
        Owner*          owner()  const  { return owner_; }

        struct ref_proxy {
            Owner *container_;
            difference_type pos_;

            ref_proxy(Owner *o, difference_type pos) : container_(o), pos_(pos) {}

            ref_proxy& operator=(const OutVT& rhs)     { container_->put_at(pos_, rhs);  return *this;}
            ref_proxy& operator=(const ref_proxy& rhs) { container_->put_at(pos_, rhs.container_->get_at(rhs.pos_)); return *this;}
            operator InVT() const					   { return container_->get_at(pos_); }
        };

        struct deref_proxy {
            Owner *container_;
            difference_type pos_;

            deref_proxy(Owner *o, difference_type pos) : container_(o), pos_(pos) {}
            operator InVT() const { return container_->reference_at(pos_); }
        };

        //typedef ref_proxy dereference_type; // if OutV==void, deref type = InVT or const InVT&
        using dereference_type = std::conditional_t<std::is_void_v<OutVT>, InVT, ref_proxy>;
        //typedef typename boost::mpl::if_<typename boost::is_void<OutVT>, InVT, ref_proxy>::type  dereference_type;

        
        dereference_type operator*() const {
          using helper_type = std::conditional_t<std::is_void_v<OutVT>, deref_proxy, ref_proxy>;
            //typedef typename boost::mpl::if_<typename boost::is_void<OutVT>, deref_proxy, ref_proxy>::type  helper_type;
            return helper_type(owner_, current_);
        }
        /*
        // does not work with VS7.1
        typename boost::enable_if<boost::is_void<OutVT>, dereference_type>::type operator*() const {
            return owner_->get_at(current_);
        }

        typename boost::disable_if<boost::is_void<OutVT>, dereference_type>::type operator*() const {
            return ref_proxy(owner_, current_);
        }
        */


    private:
        difference_type    current_;
        Owner             *owner_;
    };


    template <typename Owner, typename InVT, typename OutVT>
    typename random_index_iterator<Owner, InVT, OutVT>::difference_type
    operator-(const random_index_iterator<Owner, InVT, OutVT>& lhs, 
              const random_index_iterator<Owner, InVT, OutVT>& rhs) {
        return lhs.current() - rhs.current();
    }

    template <typename Owner, typename InVT, typename OutVT>
    random_index_iterator<Owner, InVT, OutVT>
    operator-(const random_index_iterator<Owner, InVT, OutVT>& lhs, 
              typename random_index_iterator<Owner, InVT, OutVT>::difference_type rhs) {
        return random_index_iterator<Owner, InVT, OutVT>(lhs.owner(), lhs.current() - rhs);
    }

    template <typename Owner, typename InVT, typename OutVT>
    random_index_iterator<Owner, InVT, OutVT>
    operator+(const random_index_iterator<Owner, InVT, OutVT>& lhs, 
              typename random_index_iterator<Owner, InVT, OutVT>::difference_type rhs) {
        return random_index_iterator<Owner, InVT, OutVT>(lhs.owner(), lhs.current() + rhs);
    }

    template <typename Owner, typename InVT, typename OutVT>
    random_index_iterator<Owner, InVT, OutVT>
    operator+(typename random_index_iterator<Owner, InVT, OutVT>::difference_type lhs, 
              const random_index_iterator<Owner, InVT, OutVT> & rhs) {
        return random_index_iterator<Owner, InVT, OutVT>(rhs.owner(), rhs.current() + lhs);
    }

#if 0
    template <typename InUnaryFunc, typename OutUnaryFunc>
    class function_random_iterator {
        typedef function_random_iterator self_type;
        typedef typename InUnaryFunc::argument_type     put_value_type;
        typedef typename OutUnaryFunc::result_type      get_value_type;

    public:
        typedef std::random_access_iterator_tag     iterator_category;
        typedef typename InUnaryFunc::argument_type value_type;  // TODO use boost type traits
        typedef ::ptrdiff_t                         difference_type;
        typedef value_type*                         pointer;
        typedef value_type&                         reference;

        explicit function_random_iterator(difference_type pos, const InUnaryFunc &in_func, const OutUnaryFunc& out_func)
            : current_(pos), in_(in_func), out_(out_func) {}

        function_random_iterator(const self_type& rhs);
        self_type& operator=(const self_type& rhs) {
            current_  = rhs.current_;
            in_func_  = rhs.in_func_;
            out_func_ = rhs.out_func_;
        }

        struct ref_proxy {

            operator=(const put_value_type& rhs) {
                in_func_(current_, rhs);
                
            }

            operator=(const ref_proxy& rhs) {
                in_func_(current_, rhs.out_func_(rhs.current_));
            }

            operator get_value_type() const {
                return out_func_(current_);
            }
        };
        

        ref_proxy operator*();

        self_type& operator++() {
            ++current_;
            return *this;
        }

        self_type  operator++(int) {
            self_type tmp(*this);
            ++(*this);
            return tmp;
        }

        self_type& operator--() {
            --current_;
            return *this;
        }


        self_type  operator--(int) {
            self_type tmp(*this);
            --(*this);
            return tmp;
        }

        bool operator==(const self_type& rhs) const {
            return current_ == rhs.current_;
        }

        bool operator!=(const self_type& rhs) const {
            return current_ != rhs.current_;
        }

        bool operator< (const self_type& rhs) const {
            return current_ < rhs.current_;
        }

        bool operator<=(const self_type& rhs) const {
            return current_ <= rhs.current_;
        }

        bool operator> (const self_type& rhs) const {
            return current_ > rhs.current_;
        }

        bool operator>=(const self_type& rhs) const {
            return current_ >= rhs.current_;
        }


    private:
        difference_type     current_;
        InUnaryFunc         in_func_;
        OutUnaryFunc        out_func_;
    };

    template <typename InUnaryFunc, typename OutUnaryFunc>
    typename function_random_iterator<InFunc,OutFunc>::difference_type
    operator-(const function_random_iterator<InFunc,OutFunc>& lhs, const function_random_iterator<InFunc,OutFunc>& rhs) {
    }

    template <typename InUnaryFunc, typename OutUnaryFunc>
    function_random_iterator<InFunc,OutFunc>
    operator-(const function_random_iterator<InFunc,OutFunc>& lhs, typename function_random_iterator<InFunc,OutFunc>::difference_type rhs) {
    }

    template <typename InUnaryFunc, typename OutUnaryFunc>
    function_random_iterator<InFunc,OutFunc>
    operator+(const function_random_iterator<InFunc,OutFunc>& lhs, typename function_random_iterator<InFunc,OutFunc>::difference_type rhs) {
    }

    template <typename InUnaryFunc, typename OutUnaryFunc>
    function_random_iterator<InFunc,OutFunc>
    operator+(typename function_random_iterator<InFunc,OutFunc>::difference_type lhs, const function_random_iterator<InFunc,OutFunc>& rhs) {
    }
#endif
}