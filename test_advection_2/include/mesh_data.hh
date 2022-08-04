#ifndef MESH_DATA_HH__
#define MESH_DATA_HH__

#include <vt/transport.h>
#include <array>
#include "mesh_array.hh"

#define DBG_ASSERT_BNDCHECK

namespace advect_1d {

template<std::size_t ntls_, bool have_rhs> 
struct mesh_data_t { 
    private: 
    std::size_t npoints ;
    std::array<mesh_array_t<double>, ntls_> data ;
    mesh_array_t<double> rhs_data    ;

    public:
    /**
     * @brief Default ctor
     * data objects are initialised to nullptr's 
     */
    mesh_data_t(): 
    npoints(0), data(), 
    rhs_data() 
    { } ;

    mesh_data_t(std::size_t npoints__) : 
    npoints(npoints__)
    {
        for(std::size_t ii=0; ii<ntls_; ++ii)
            data[ii] = mesh_array_t<double>(npoints ) ;

        if( have_rhs )
            rhs_data = mesh_array_t<double>(npoints) ;
    }

    double& operator() (std::size_t idx) { 
        #ifdef DBG_ASSERT_BNDCHECK
        assert(ntls_ > 0);
        assert(idx < npoints) ;
        #endif 
        return data[0](idx) ;     
    }

    mesh_array_t<double>& operator[](std::size_t tl_idx_) { 
        #ifdef DBG_ASSERT_BNDCHECK
        assert(tl_idx_ < ntls_ ) ;
        #endif 
        return data[tl_idx_] ; 
    }   

    void resize(const std::size_t & new_size__ ) {
        npoints = new_size__;

        for( int ii=0; ii<ntls_; ++ii)
            data[ii].resize(npoints) ;
        
        if( have_rhs )
            rhs_data.resize(npoints) ;

    }

    double& rhs(std::size_t idx) { 
        static_assert(have_rhs,
        "Trying to access rhs of non-evolved variable!\n") ;
        #ifdef DBG_ASSERT_BNDCHECK
        assert(idx < npoints) ;
        #endif 
        return rhs_data(idx); 
    }

    std::size_t size() { return npoints ; }

    std::size_t get_n_timelevels() { return data.size() ; }

    void rotate_timelevels() {
        // this is implemented in a 
        // clever way under the hood
        if( ntls_ == 0 )
            return ;
        
        for( int ii=0; ii< data.size(); ii++ )
            data[ii+1] = data[ii] ;
    }

    mesh_data_t<ntls_, have_rhs>& operator=(const mesh_data_t<ntls_, have_rhs>& rhs_ ) {
        if( this == &rhs_ ){
            return *this ;
        }

        data = rhs_.data ;

        if( have_rhs )
            rhs_data = rhs_.rhs_data ;
        
        return *this ;
    }

    template< typename Serializer >
    void serialize( Serializer& s ) {
        s | have_rhs ;
        s | npoints ;
        s | data ;
        if( have_rhs )
            s | rhs_data ;
    }
};
template<std::size_t ntls_ >
using evolved_field_t = mesh_data_t<ntls_, true> ;
template<std::size_t ntls_ >
using aux_field_t = mesh_data_t<ntls_, false> ;

}

#endif