#ifndef MESH_DATA_HH__
#define MESH_DATA_HH__

#include <vt/transport.h>

#define DBG_ASSERT_BNDCHECK

namespace advect_1d {

struct mesh_data_t { 
    private: 
    std::size_t npoints ;

    struct vt_deallocator {
        void operator() (void* buf) const { vt::thePool()->dealloc(buf);  }
    };

    template<typename T>
    using vta_unique_ptr = std::unique_ptr<double, vt_deallocator> ;

    vta_unique_ptr<double> data    ;
    vta_unique_ptr<double> data_old ; 
    vta_unique_ptr<double> rhs_data ;  
    public:
    
    mesh_data_t(): data(nullptr), data_old(nullptr), rhs_data(nullptr) {} ;

    mesh_data_t(std::size_t npoints__) : npoints(npoints__) {
        data = vta_unique_ptr<double>( (double*) vt::thePool()->alloc(sizeof(double)*npoints )) ;
        data_old = vta_unique_ptr<double>( (double*) vt::thePool()->alloc(sizeof(double)*npoints )) ;
        rhs_data = vta_unique_ptr<double>( (double*) vt::thePool()->alloc(sizeof(double)*npoints )) ;
    }

    double& operator() (std::size_t idx) { 
        #ifdef DBG_ASSERT_BNDCHECK
        assert(idx < npoints) ;
        #endif 
        return data.get()[idx];     
    }

    double old(std::size_t idx) { 
        #ifdef DBG_ASSERT_BNDCHECK
        assert(idx < npoints) ;
        #endif 
        return data_old.get()[idx]; 
    }

    double& rhs(std::size_t idx) { 
        #ifdef DBG_ASSERT_BNDCHECK
        assert(idx < npoints) ;
        #endif 
        return rhs_data.get()[idx]; 
    }

    void rotate_timelevels() {
        memcpy(data_old.get(), data.get(), sizeof(double)*npoints ) ;
    }
};
}

#endif