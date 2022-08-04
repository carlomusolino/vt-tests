#ifndef LOOP_UTILS_HH__
#define LOOP_UTILS_HH__

#include<stdlib.h>

namespace advect_1d { namespace utils {

template<typename F>
static inline __attribute__((always_inline)) 
void avt_simd_loop(std::size_t lo_, std::size_t hi_,
                    F&& func ) 
{
    #pragma omp simd
    for(std::size_t ii=lo_; ii<hi_; ++ii)
        func(ii) ; 

}

}
}

#endif