#ifndef MESH_HH__
#define MESH_HH__

#include <vt/transport.h>
#include "node.hh"
#include "recon_utils.hh"
#include "mesh_data.hh"

#define A_VEL__ 1

namespace advect_1d {

struct mesh_t : vt::Collection<mesh_t, vt::Index1D> {
    private:
    static constexpr const std::size_t output_every = 10 ;

    std::size_t npoints_, nobj_, nghost_=2; 
    std::size_t iter_;
    double time_ ;
    std::size_t maxiter_{10000}; 
    double max_time_{10000}; 
    double CFL_{0.4} ;

    std::size_t n_neighbors_{2} ;
    std::size_t n_recv_{0} ;
    std::size_t substep_count_{0} ;

    node_proxy_t obj_proxy_ ;

    std::vector<double> x ;
    double dx ;
    double dt ;

    mesh_data_t u ;

    template<typename F> 
    inline __attribute__((always_inline)) 
    void loop_interior(F&& func) {
        for( std::size_t ii=nghost_; ii<npoints_-nghost_; ++ii) 
            func(ii) ;
    }

    template<typename F> 
    inline __attribute__((always_inline)) 
    void loop_all(F&& func) {
        for( std::size_t ii=0; ii<npoints_; ++ii) 
            func(ii) ;
    }

    template<int stagger_l, int stagger_r, typename F> 
    inline __attribute__((always_inline)) 
    void loop_staggered(F&& func) {
        for( std::size_t ii=nghost_+stagger_l; ii<npoints_-nghost_ + stagger_r; ++ii) 
            func(ii) ;
    }

    template<int stagger_l, int stagger_r, typename F> 
    inline __attribute__((always_inline)) 
    void loop_inverted(F&& func) {
        for( std::size_t ii=npoints_-nghost_ + stagger_r; ii>=nghost_+stagger_l; --ii) 
            func(ii) ;
    }

    public: 

    explicit mesh_t() :
    iter_(0), npoints_(1), nobj_(1),
    time_(0), maxiter_(1000), max_time_(10000),
    CFL_(0.4), x(), dx(1.), dt(1.), substep_count_(0),
    n_neighbors_(2), n_recv_(0)
    { } 

    using void_msg_t = vt::CollectionMessage<mesh_t> ;

    struct mesh_msg_t: vt::CollectionMessage<mesh_t> {
        std::size_t npoints, nobj, maxiter ;
        double max_time, CFL, dx ;

        node_proxy_t parent_proxy ;

        mesh_msg_t() = default ; 

        mesh_msg_t( std::size_t const np, std::size_t const no, std::size_t const m_i,
        double const m_t, double const CFL__, double const dx__, node_proxy_t proxy__) :
        npoints(np), nobj(no), maxiter(m_i), max_time(m_t), 
        CFL(CFL__), dx(dx__), parent_proxy(proxy__)
        { }
    }; 

    struct red_msg_t : vt::collective::ReduceTMsg<double> {
        red_msg_t() = default ; 
        explicit red_msg_t( double in_val ) : ReduceTMsg<double>(in_val) {} 
    } ;

    void output_done_cback_(red_msg_t* msg) {
        double umax = msg->getConstVal() ;

        // This is threadsafe, trust me!
        std::string fname("u_max.dat") ;
        std::ofstream file ;
        file.open(fname, std::ios_base::app) ;
        file << iter_ << '\t' << time_ << '\t' << umax << '\n';
        file.close() ;
        
        auto const it_max_reached = iter_ > maxiter_ ;
        auto const t_max_reached = time_ > max_time_ ;

        if( it_max_reached or t_max_reached ){
            auto const output = it_max_reached ?
            "Maximum number of iterations reached, terminating\n\n" :
            "Maximum time reached, terminating\n\n" ;
            obj_proxy_.broadcast<node_t::work_finished_msg, &node_t::work_finished_cback_>() ;
        } else {
            fmt::print(">> it {} t {} u_max {}\n", iter_,time_, umax) ;
        }
    }

    // Initialisation methods
    // The first one assumes all the members that need to come from the RTS
    // e.g. npoints, nobjects etc. are already set.
    // The second one takes a message as input and sets these values:
    // it will be invoked by the RTS. 
    void init_() { 
        const auto my_idx = getIndex().x() ;
        
        x.assign(npoints_, 0.) ;
        u = mesh_data_t(npoints_) ;

        double xL =  ( (my_idx) * (npoints_-2*nghost_) ) * dx ; // CHECK
        for( std::size_t ii=nghost_; ii<npoints_-nghost_; ++ii){
            x[ii] = xL + ii*dx ;
        }
        for( int ii=0; ii<nghost_ ; ii++){
            x[nghost_-ii-1] = x[nghost_] - ii*dx ;
            x[npoints_-nghost_ +ii] = x[npoints_-nghost_] + ii * dx ;
        }

        // Initial data: Gaussian Pulse 
        for( std::size_t ii=nghost_; ii<npoints_-nghost_; ++ii){
            u(ii) = std::exp( - (x[ii]-0.5) * (x[ii]-0.5) / 0.1 / 0.1  ) ;
        }
        u.rotate_timelevels() ;

        dt = CFL_ * dx ;
    }

    void init(mesh_msg_t* msg) {
        nghost_ = 2;
        npoints_ = msg->npoints + 2*nghost_;
        nobj_ = msg->nobj ;
        maxiter_ = msg->maxiter ;
        max_time_ = msg->max_time ; 
        CFL_ = msg->CFL ;
        obj_proxy_ = msg->parent_proxy ;
        dx = msg->dx ;
        n_neighbors_ = 2 ;
        n_recv_ = 0; 
        init_() ;
    }

    // Compute RHS based on minmod
    // second order FV scheme
    void calc_rhs() {
        // allocate a couple of temporaries
        std::vector<double> u_r(npoints_) ;
        std::vector<double> u_l(npoints_) ; 
        std::vector<double> flux(npoints_) ; 

        loop_staggered<-1,1>(
            [&] (const std::size_t & idx) {
                auto dp = u(idx+1) - u(idx) ;
                auto dm = u(idx) - u(idx-1) ;            
                u_r[idx] = u(idx) + 0.5*mc2(dp,dm) ;
                dp = u(idx-1) - u(idx) ;
                dm = u(idx) - u(idx+1) ;
                u_l[idx] = u(idx) + 0.5*mc2(dp,dm) ; 
        }
        );

        loop_inverted<0,1>(
            [&] (const std::size_t& idx) {
                auto const tmp = u_l[idx] ;
                u_l[idx] = u_r[idx-1] ;
                u_r[idx] = tmp;
            }
        ) ;

        // Now we have the reconstructed states
        // compute the rhs 
        loop_staggered<0,1>(
            [&] (const std::size_t & idx) {
                flux[idx] = (A_VEL__ < 0) ? A_VEL__ * u_r[idx] : A_VEL__ * u_l[idx] ; 
            }
        ) ;

        double const inv_dx = 1./dx; 
        loop_interior(
            [&] (const std::size_t & idx) {
                u.rhs(idx) =  inv_dx*(flux[idx] - flux[idx+1]); 
            }
        ) ; 

    } 

    

    struct gz_msg_t : vt::CollectionMessage<mesh_t> {
        using MessageParentType = vt::CollectionMessage<mesh_t>;
        vt_msg_serialize_if_needed_by_parent_or_type2(vt::IdxBase, std::vector<double>);

        gz_msg_t() = default;

        gz_msg_t(vt::IdxBase const& in_index, std::vector<double> const& ref) :
        vt::CollectionMessage<mesh_t>(),
        from_index(in_index), val(ref)
        { }

        template <typename Serializer>
        void serialize(Serializer& s) {
            MessageParentType::serialize(s);
            s | from_index;
            s | val;
            }

        vt::IdxBase from_index = 0;
        std::vector<double> val ;
    };

    void sync_data()
    {  
        if( nobj_ == 1 ) {
            for( std::size_t ii=0; ii<nghost_; ii++){

                u(npoints_-nghost_+ii) = u(nghost_+ii);
                u(nghost_-1-ii) = u(npoints_-nghost_-1-ii);

               /*
                u[npoints_-1-ii] = u[npoints_-nghost_-1] ;
                u[nghost_-1-ii]  = u[nghost_] ;
                */
            }
            if ( substep_count_ < 2 )
                do_substep() ;
        } else {

        const vt::IdxBase my_idx = getIndex().x() ;

        auto proxy = this->getCollectionProxy() ;
        std::vector<double> gz_out(nghost_) ;
        if( my_idx == 0 ) {
            for( std::size_t ii=0; ii<nghost_; ii++){ 
                gz_out[ii] = u(nghost_+ii) ;
            }
            proxy[nobj_-1].send<gz_msg_t, &mesh_t::exchange>(
                    my_idx, gz_out 
            ) ;
        }

        if( my_idx == nobj_-1 ) {
            for( std::size_t ii=0; ii<nghost_; ii++){ 
                gz_out[ii] = u(npoints_-nghost_-1-ii) ;
            }
            proxy[0].send<gz_msg_t, &mesh_t::exchange>(
                my_idx, gz_out 
            ) ;
        }

        if( my_idx > 0 ){
            for( std::size_t ii=0; ii<nghost_; ii++){ 
                gz_out[ii] = u(nghost_+ii) ;
            }
            proxy[my_idx-1].send<gz_msg_t, &mesh_t::exchange>(
                my_idx, gz_out 
            ) ; 
        }

        if( my_idx < nobj_-1 ) {
            for( std::size_t ii=0; ii<nghost_; ii++){ 
                gz_out[ii] = u(npoints_-nghost_-1-ii) ;
            }
            proxy[my_idx+1].send<gz_msg_t, &mesh_t::exchange>(
                my_idx, gz_out 
            ) ;
        }
        }
    }

    // Do one substep of a midpoint RK2 integrator  
    void do_substep() {
        switch (substep_count_) {
        case 0 :
            substep_count_++ ;
            calc_rhs() ;
            loop_interior (
                [&] ( const std::size_t & idx) {
                    u(idx) += 0.5*dt * u.rhs(idx) ;
                }
            ) ;
            sync_data();
            break ;
        case 1:
            substep_count_++ ;
            calc_rhs() ;
            loop_interior (
                [&] ( const std::size_t & idx) {
                    u(idx) = u.old(idx) + dt * u.rhs(idx) ;
                }
            ) ;
            sync_data();
            break ;
        default:
            break ;
        }
    }

    // output 1d data as well as max ;
    void do_output(void_msg_t* msg) {

        const vt::IdxBase my_idx = getIndex().x() ;

        if( iter_ % output_every == 0) {
            auto fname = fmt::format( "u.1d.{}.dat",my_idx) ;
            std::ofstream outfile;
            outfile.open(fname, std::ios_base::app) ;
            outfile << iter_ << '\t' << time_ << '\t';
            loop_interior (
                [&] ( const std::size_t& idx){
                    outfile << u(idx) << '\t' ;
                } 
            );
            outfile << '\n';
            outfile.close() ;
        }

        double umax=u(nghost_) ;
        loop_interior (
                [&] ( const std::size_t & idx) {
                    umax = std::max(umax, u(idx)) ;
                } 
        ) ;
        auto proxy = this->getCollectionProxy() ;
        auto cb = vt::theCB()->makeSend<
                mesh_t, red_msg_t, &mesh_t::output_done_cback_
                >(proxy[0]) ;
        auto msg2 = vt::makeMessage<red_msg_t>(umax) ;
        proxy.reduce<vt::collective::MaxOp<double>>(msg2.get(), cb);
        
    }

    void exchange(gz_msg_t* msg) {

        const vt::IdxBase my_idx = getIndex().x() ;

        std::vector<double>  gz_in(nghost_) ;
        // Periodic BCs 
        if ( my_idx == 0 && msg->from_index == nobj_ - 1 ) {
            gz_in = msg->val ;
            for(std::size_t ii=0; ii<nghost_; ii++){ 
                this->u(nghost_-1-ii) = gz_in[ii] ;
            }
            n_recv_ ++ ;
        } else if ( my_idx == nobj_ -1 && msg->from_index == 0 ) {
            gz_in = msg->val ;
            for(std::size_t ii=0; ii<nghost_; ii++){ 
                this->u(npoints_-nghost_+ii) = gz_in[ii] ;
            }
            n_recv_ ++ ;
        } else if ( my_idx > msg->from_index) {
            gz_in = msg->val ;
            for(std::size_t ii=0; ii<nghost_; ii++){ 
                this->u(nghost_-1-ii) = gz_in[ii] ;
            }
            n_recv_ ++ ;
        } else if ( my_idx < msg->from_index ) {
            gz_in = msg->val ;
            for(std::size_t ii=0; ii<nghost_; ii++){ 
                this->u(npoints_-nghost_+ii) = gz_in[ii] ;
            }
            n_recv_ ++ ;
        }

        if( n_recv_ == n_neighbors_ ) {
            n_recv_=0; 
            if( substep_count_ < 2 )
                do_substep() ;
        }
    }
    // This will be called in an epoch
    // It simply sends gz data 
    void do_step(void_msg_t* msg) {
        time_ += dt ;
        iter_ ++ ;
        substep_count_ = 0 ;
        // rotate timelevels 
        u.rotate_timelevels() ;
        sync_data();
    }


} ;

}

#endif