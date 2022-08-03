#ifndef LOOP_UTILS_HH__
#define LOOP_UTILS_HH__

#include<stdlib.h>

namespace advect_1d { namespace utils {

template< std::size_t i_min, std::size_t i_max, typename F> 
static inline __attribute__((always_inline)) 
void loop_simple(F&& func) {
    for( std::size_t ii=i_min; ii<i_max; ++ii) 
        func(ii) ;
}

template< std::size_t i_min, std::size_t i_max, typename F> 
static inline __attribute__((always_inline)) 
void loop_inverted(F&& func) {
    for( std::size_t ii=i_max; ii>=i_min; --ii) 
        func(ii) ;
}

}
}

#endif