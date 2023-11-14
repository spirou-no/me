#include <Utility/algorithm.h>

#include <vector>
#include <list>

// make compile tests with as strict iterators as possible.  E.g. use forward or input where spec-ed.
// for random_xxx algorithms, try combinations (in, fwd), (fwd, out), (in,out)->error

// code snippets to make sure code compiles
namespace {

using namespace oyrke::algorithm;

#if 0
template <typename T>
T random(T) { return T(1); }
#endif

struct user_def {};
bool operator<(int, user_def) { return false; }
bool operator<(user_def, int) { return false; }

template <typename It>
class input_iterator {
public:
    typedef typename std::iterator_traits<It>::value_type       value_type;
    typedef typename std::iterator_traits<It>::difference_type  difference_type;
    typedef difference_type    distance_type;
    typedef typename std::iterator_traits<It>::value_type       reference;
    typedef typename std::iterator_traits<It>::pointer          pointer;
    typedef std::input_iterator_tag                             iterator_category;
    
    input_iterator(It it) : it_(it) {}
    input_iterator(const input_iterator& rhs) : it_(rhs.it_) {}
    input_iterator& operator=(const input_iterator& rhs) { it_ = rhs.it_; return *this; }
    value_type operator*() const { return *it_; }
    input_iterator& operator++() { ++it_; return *this; }
    input_iterator  operator++(int) { input_iterator tmp(*this); ++(*this); return tmp; }
    bool operator==(const input_iterator& rhs) const { return it_ == rhs.it_; }
    bool operator!=(const input_iterator& rhs) const { return it_ != rhs.it_; }
    
private:
    It it_;    
};

template <typename It>
class forward_iterator {
public:
    typedef typename std::iterator_traits<It>::value_type       value_type;
    typedef typename std::iterator_traits<It>::difference_type  difference_type;
    typedef typename difference_type                            distance_type;
    typedef typename std::iterator_traits<It>::reference        reference;
    typedef typename std::iterator_traits<It>::pointer          pointer;
    typedef std::forward_iterator_tag                           iterator_category;
    
    forward_iterator(It it) : it_(it) {}
    forward_iterator(const forward_iterator& rhs) : it_(rhs.it_) {}
    forward_iterator& operator=(const forward_iterator& rhs) { it_ = rhs.it_; return *this; }
    reference operator*() const { return *it_; }
    forward_iterator& operator++() { ++it_; return *this; }
    forward_iterator  operator++(int) { forward_iterator tmp(*this); ++(*this); return tmp; }
    bool operator==(const forward_iterator& rhs) const { return it_ == rhs.it_; }
    bool operator!=(const forward_iterator& rhs) const { return it_ != rhs.it_; }
    
private:
    It it_;    
};


void
test_compile() {
  std::vector<int> v;
  std::list<int>   l;
  std::vector<__int64>  v64;
  std::vector<user_def> vud;

  typedef std::vector<int>::iterator      vector_it;
  typedef std::vector<__int64>::iterator  vector_it64;
  typedef std::vector<user_def>::iterator vector_itud;
  typedef std::list<int>::iterator        list_it;

  std::pair<vector_it, list_it> match = first_set_intersection(v.begin(), v.end(), l.begin(), l.end());
  std::pair<vector_itud, list_it> match2 = first_set_intersection(vud.begin(), vud.end(), l.begin(), l.end());
  for_each_if(v.begin(), v.end(), [](auto) {}, [](auto x) { return std::less<int>()(x, 2); });
  std::copy_if(v.begin(), v.end(), l.begin(), [](auto x) { return x < 2; } );
  auto less2 = [](auto x) { return x < 2; };
    copy_backward_if(v.begin(), v.end(), l.end(), less2); 
    while_do(v.begin(), v.end(), [](auto) {}, less2);
    transform_if(v.begin(), v.end(), l.begin(), std::negate<int>(), less2); 
    transform_if(v.begin(), v.end(), l.begin(), l.begin(), std::plus<int>(), std::less<int>()); 
    
    vector_it vit = min_element_if(v.begin(), v.end(), less2); 
    vit = max_element_if(v.begin(), v.end(), less2); 

    bool yes = contains(v.begin(), v.end(), 5);
    yes = contains_if(v.begin(), v.end(), less2);
    vit = lower_bound_hinted(v.begin(), v.end(), 5, v.begin());
    vit = lower_bound_hinted(v.begin(), v.end(), 5, v.begin(), v.begin()+3);
    vit = upper_bound_hinted(v.begin(), v.end(), 5, v.begin());
    vit = upper_bound_hinted(v.begin(), v.end(), 5, v.begin(), v.begin()+3);
    
    std::pair<vector_it, vector_it> eqr = equal_range_hinted(v.begin(), v.end(), 5, v.begin());
    eqr = equal_range_hinted(v.begin(), v.end(), 5, v.begin(), v.begin()+3);

    yes = binary_search_hinted(v.begin(), v.end(), 5, v.begin());
    yes = binary_search_hinted(v.begin(), v.end(), 5, v.begin(), v.begin()+3);


    vit = random_1(v.begin(), v.end());
    vit = random_1_if(v.begin(), v.end(), less2); 
    random_n(v.begin(), v.end(), 3, l.begin());
    random_n_of_m(v.begin(), 5, 3, l.begin());
    random_n_if(v.begin(), v.end(), 3, l.begin(), less2);
    
    /* TODO
    o binary_find
    o seeded_binary_find
    */
    
    
}


void
test_compile_stress() {
    std::vector<int> v;
    std::list<int>   l;
    std::vector<__int64>  v64;
    std::vector<user_def> vud;
    
    typedef std::vector<int>::iterator      vector_it;
    typedef std::vector<__int64>::iterator  vector_it64;
    typedef std::vector<user_def>::iterator vector_itud;
    typedef std::list<int>::iterator        list_it;
    
    typedef std::vector<int>::iterator      random_it_t;
    typedef std::list<int>::iterator        bidir_it_t;
    typedef std::back_insert_iterator<std::vector<int> > out_it_t;
    typedef input_iterator<random_it_t>     in_it_t;
    typedef forward_iterator<random_it_t>   fwd_it_t;
    
    random_it_t rnd1 = v.begin();
    random_it_t rnd2 = v.end();
    bidir_it_t bi1   = l.begin();
    bidir_it_t bi2   = l.end();
    fwd_it_t   fwd1  = fwd_it_t(v.begin());
    fwd_it_t   fwd2  = fwd_it_t(v.end());
    in_it_t    in1   = in_it_t(v.begin());
    in_it_t    in2   = in_it_t(v.end());
    out_it_t   out   = out_it_t(v);
    std::negate<int> unop;
    auto noop = [](auto) {};
    auto unpred = [](auto x) { return x < 2; };
    
    std::pair<vector_it, list_it> match = first_set_intersection(v.begin(), v.end(), l.begin(), l.end());
    
    // try intersection comparing different types which overload op< (i.e. ok),
    std::pair<vector_itud, list_it> match2 = first_set_intersection(vud.begin(), vud.end(), l.begin(), l.end());
    
    first_set_intersection(in1, in2, fwd1, fwd2);
    
    for_each_if(in1, in2, noop, unpred); 
    std::copy_if(in1, in2, out, unpred);
    copy_backward_if(bi1, bi2, bi1, unpred);
    while_do(in1, in2, noop, unpred);
    transform_if(in1, in2, out, unop, unpred);
    transform_if(in1, in2, in1, out, std::plus<int>(), std::less<int>()); 
    
    min_element_if(in1, in2, unpred);
    max_element_if(in1, in2, unpred);

    bool yes = contains(in1, in2, 5);
    yes = contains_if(in1, in2, unpred);
    
    lower_bound_hinted(bi1, bi2, 5, bi1);
    lower_bound_hinted(bi1, bi2, 5, bi1, bi1);
    upper_bound_hinted(bi1, bi2, 5, bi1);
    upper_bound_hinted(bi1, bi2, 5, bi1, bi1);
    
    equal_range_hinted(bi1, bi2, 5, bi1);
    equal_range_hinted(bi1, bi2, 5, bi1, bi1);
    
    binary_search_hinted(bi1, bi2, 5, bi1);
    binary_search_hinted(bi1, bi2, 5, bi1, bi1);


    random_1(in1, in2);   // impl will iterate through range
    random_1(fwd1, fwd2);  // impl will compute distance fwd1..fwd2, then pick a random index
    random_1_if(in1, in2, unpred); 
    random_n(in1, in2, 3, fwd1);
    random_n(fwd1, fwd2, 3, out);
    //random_n(in1, in2, 3, out);  // error, can not advance(out, ...)
    //random_n(in1, in2, 3, in1);  // error, can not write to in1
    random_n_of_m(in1, 5, 3, out);
    random_n_if(in1, in2, 3, fwd1, unpred); 
    
    /* TODO
    o binary_find
    o seeded_binary_find
    */
    
    
}

}
