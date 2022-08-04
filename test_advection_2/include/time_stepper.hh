#ifndef AVT_TIME_STEPPER_HH__
#define AVT_TIME_STEPPER_HH__

#include "mesh.hh"
#include "mesh_data.hh"
#include "mesh_array.hh"
#include "loop_utils.hh"
namespace advect_1d {

/** 
 * @brief base class for timestepper
 */ 
template<std::size_t scratch_size_> 
class time_stepper_t {
    protected: 
    // delta t 
    double dt ;
    // number of substeps 
    std::size_t n_substeps_ ;
    // is the timestep done?
    bool done_ ;
    // scratch space for old rhs's
    std::vector<mesh_array_t> sratch_ ;;
    // number of evolved variables,
    // same across simulation
    // probably known at compile time (?)
    static std::size_t num_evolved ;
    
    public:

    virtual time_stepper_t(
        std::size_t n_substeps__ ,
        double const& dt__ ):
        dt(dt__),
        n_substeps_(n_substeps__),
        done(false) 
        {
            scratch_.resize(scratch_size * num_evolved) ;
        }

    // Pure virtual method for doing a substep 
    virtual void do_substep(evolved_field_t& ) = 0 ; 

    bool is_done() { return done ; }

    virtual ~time_stepper_t() ;
} ;

class rk2_stepper_t: public time_stepper_t<0> {
    // counter for current substep 
    std::size_t substep_counter_ ;
    public:
    rk2_stepper_t() : time_stepper_t<0>(2), substep_counter_(0)
    { } 

    void do_substep( evolved_field_t& var ){
        switch(substep_counter_) {
        case 0: 
            substep_counter_++ ;
            utils::avt_simd_loop(0, var.size(),
            [&] ( const std::size_t & idx )
            {
                var(idx) += 0.5 * dt * var.rhs(idx) ;
            } 
            ) ;
            break ;
        case 1:
            utils::avt_simd_loop(0, var.size(),
            [&] ( const std::size_t & idx )
            {
                var(idx) += var[1](idx) +  dt * var.rhs(idx) ;
            }
            ) ;
            this->done_ = true ;
            substep_counter_=0 ;
            break;
        default :
            break;
        }

    }


} ;


}

#endif 