#include <vt/transport.h>
#include <vt/runnable/invoke.h>

#include <cstdlib>
#include <cassert>
#include <iostream>

static constexpr std::size_t const d_nrow = 8;
static constexpr std::size_t const d_nobj = 4;
static constexpr double const d_tol = 1e-02 ;

struct node_t {
  bool done_ = false ;
  struct work_done_msg_t : vt::Message {} ;
  void work_done_cback_ (work_done_msg_t*) { done_ = true ; } ;
  bool is_work_done() { return done_ ; } 
} ;

using node_proxy_t = vt::objgroup::proxy::Proxy<node_t> ;

struct jacobi_1d_t : vt::Collection<jacobi_1d_t, vt::Index1D> {
private:

  std::vector<double> tcur_, told_ ;
  std::vector<double> rhs_ ;
  size_t iter_ = 0 ;
  size_t msg_rec_ =0, tot_rec_ =0 ;
  size_t nobj_ = 1;
  size_t nrows_ = 1;
  size_t maxiter_ = 10 ;

  node_proxy_t obj_proxy_ ;


public:

  explicit jacobi_1d_t () :
    tcur_(), told_(), rhs_(), iter_(0),
    msg_rec_(0), tot_rec_(0),
    nobj_(1), nrows_(1), maxiter_(10)
  {}


  using blank_msg_t = vt::CollectionMessage<jacobi_1d_t> ;

  struct jacobi_msg_t : vt::CollectionMessage<jacobi_1d_t> {
    size_t nobj, nrow ;
    size_t maxiter ;
    node_proxy_t obj_proxy ;

    jacobi_msg_t() = default ;

    jacobi_msg_t( const size_t nobj__, const size_t nrow__,
		  const size_t itermax__, const node_proxy_t myproxy__):
      nobj(nobj__), nrow(nrow__), maxiter(itermax__), obj_proxy(myproxy__)
    {}

  } ;


  struct reduce_msg_t : vt::collective::ReduceTMsg<double> {
    reduce_msg_t() = default ;
    reduce_msg_t(double in__) : ReduceTMsg<double>(in__) {} 

  } ;


  void check_complete_cback_( reduce_msg_t* msg) {

    double norm_res = msg->getConstVal() ;

    auto const iter_max_reached = iter_ > maxiter_ ;
    auto const norm_res_done= norm_res < d_tol ;

    if( iter_max_reached or norm_res_done ) {
      auto const message = iter_max_reached ?
        "\n Maximum Number of Iterations Reached. \n\n" :
        fmt::format("\n Max-Norm Residual Reduced by {} \n\n", d_tol);

      fmt::print(message) ;


      obj_proxy_.broadcast<node_t::work_done_msg_t, &node_t::work_done_cback_> () ;
      
    } else {

      fmt::print(" >> Iter {} || residual {} \n ", iter_, norm_res ) ;
    }
    
  }


  void do_iter() {

    ++iter_ ;

    tcur_[0] = told_[0] ;

    tcur_[nrows_+1] = told_[nrows_+1] ;

    for( std::size_t ii=1 ;  ii<=nrows_; ++ii ) {
      tcur_[ii] = 0.5 * (rhs_[ii] + told_[ii+1] + told_[ii-1] ) ;
    }

    std::copy( tcur_.begin(), tcur_.end(), told_.begin() ) ;

    double max_norm{0.} ;

    for( std::size_t ii=1 ; ii< tcur_.size() - 1; ++ii ) {
      max_norm = std::max( max_norm, std::fabs( tcur_[ii] ) ) ;
    }


    auto proxy = this->getCollectionProxy();

    auto cb_ = vt::theCB()->makeSend<
      jacobi_1d_t, reduce_msg_t, &jacobi_1d_t::check_complete_cback_
      >(proxy[0]) ;
    auto msg2 = vt::makeMessage<reduce_msg_t>(max_norm) ;
    proxy.reduce<vt::collective::MaxOp<double>>(msg2.get(), cb_) ;
							 
  }
  
  struct vec_msg_t : vt::CollectionMessage<jacobi_1d_t> {
    using msg_parent_t = vt::CollectionMessage<jacobi_1d_t> ;

    vt_msg_serialize_if_needed_by_parent_or_type1(vt::IdxBase) ;

    vec_msg_t() = default ;

    vec_msg_t( vt::IdxBase const& in_idx, double const& ref) :
      vt::CollectionMessage<jacobi_1d_t>() ,
      from_index(in_idx), val(ref)
    {}

    template<typename serializer_t>
    void serialize( serializer_t& s) {
      msg_parent_t::serialize(s) ;
      s | from_index ;
      s | val ;
    }

    vt::IdxBase from_index = 0;
    double val = 0. ;
    
  };

  void exchange(vec_msg_t* msg) {

    const vt::IdxBase my_idx = getIndex().x() ;

    if( my_idx > msg->from_index) {
      this->told_[0] = msg->val ;
      msg_rec_ ++ ;
    }

    if ( my_idx < msg->from_index) {
      this->told_[nrows_ + 1] = msg->val ;
      msg_rec_ ++ ;
    }

    if ( msg_rec_ == tot_rec_ ) {
      msg_rec_ = 0;
      do_iter() ;
    }
    
  }


  void do_iteration(blank_msg_t *msg) {
    
    //
    // Treat the particular case of 1 object
    // where no communication is needed.
    // Without this treatment, the code would not iterate.
    //
    
    if (nobj_ == 1) {
      do_iter();
      return;
    }
    //---------------------------------------
    
    //
    // Routine to send information to a different object
    //
    
    vt::IdxBase const myIdx = getIndex().x();

    //--- Send the values to the left
    auto proxy = this->getCollectionProxy();
    if (myIdx > 0) {
      proxy[myIdx - 1].send<vec_msg_t, &jacobi_1d_t::exchange>(
							       myIdx, told_[1]
							       );
    }
    
    //--- Send values to the right
    if (size_t(myIdx) < nobj_ - 1) {
      proxy[myIdx + 1].send<vec_msg_t, &jacobi_1d_t::exchange>(
							       myIdx, told_[nrows_]
							       );
    }
  }
  
  void init() {

    tcur_.assign(nrows_ + 2, 0. ) ;
    told_.assign(nrows_ + 2, 0. ) ;
    rhs_.assign(nrows_ + 2, 0. ) ;


    double h = 1. / ( nrows_ * nobj_ + 1. ) ;
    int nf = 3 * int( nrows_ * nobj_ + 1 ) / 4 ; 


    std::size_t const my_idx = getIndex().x() ;

    for( std::size_t ii=0; ii<tcur_.size(); ++ii)  {
      double x0 = ( nrows_ * my_idx + ii ) * h ;
      tcur_[ii] = sin( nf * M_PI * x0 * x0 ) ;
    }

    tot_rec_ = 2 ;

    if( my_idx == 0 ) {
      tcur_[0] = 0 ;
      tot_rec_ -- ;
    }

    if( my_idx == nobj_ - 1 ) {
      tcur_[nrows_ + 1 ] = 0. ;
      tot_rec_ -- ;
    }

    std::copy( tcur_.begin(), tcur_.end(), told_.begin() ) ;
      
  }


  void init( jacobi_msg_t* msg) {
    nobj_=msg->nobj ;
    nrows_ = msg->nrow;
    maxiter_ = msg->maxiter ;
    obj_proxy_ = msg->obj_proxy ;

    init() ;
  }

} ; // struct jacobi_1d_t



bool is_work_done( vt::objgroup::proxy::Proxy<node_t> const& proxy) {
  auto const this_node = vt::theContext()->getNode() ;
  return proxy[this_node].invoke<decltype(&node_t::is_work_done), &node_t::is_work_done> () ;
}



int main( int argc, char** argv) {

  size_t nobj = d_nobj ;
  size_t nrow = d_nrow ;
  size_t maxiter = 10 ;


  std::string name(argv[0]) ;
  
  vt::initialize(argc, argv);


  vt::NodeType this_node = vt::theContext()->getNode() ;
  vt::NodeType num_nodes = vt::theContext()->getNumNodes() ;

  if (argc == 1) {
    if (this_node == 0) {
      fmt::print(
        stderr, "{}: using default arguments since none provided\n", name
      );
    }
    nobj = d_nobj * num_nodes;
  } else if (argc == 2) {
    nobj = static_cast<size_t>(strtol(argv[1], nullptr, 10));
  }
  else if (argc == 3) {
    nobj = static_cast<size_t>(strtol(argv[1], nullptr, 10));
    nrow = static_cast<size_t>(strtol(argv[2], nullptr, 10));
  }
  else if (argc == 4) {
    nobj = static_cast<size_t>(strtol(argv[1], nullptr, 10));
    nrow = static_cast<size_t>(strtol(argv[2], nullptr, 10));
    maxiter = static_cast<size_t>(strtol(argv[3], nullptr, 10));
  }
  else {
    fmt::print(
      stderr, "usage: {} <num-objects> <num-rows-per-object> <maxiter>\n",
      name
    );
    return 1;
  }


  auto group_proxy = vt::theObjGroup()->makeCollective<node_t>("example_jacobi_1d") ;

  using base_idx_t = typename vt::Index1D::DenseIndexType;

  auto range = vt::Index1D(static_cast<base_idx_t>(nobj)) ;

  auto col_proxy = vt::makeCollection<jacobi_1d_t>( "example_jacobi_1d" )
    .bounds(range)
    .bulkInsert()
    .wait() ;

  vt::runInEpochCollective([col_proxy, group_proxy, nobj, nrow, maxiter] {
			     col_proxy.broadcastCollective<jacobi_1d_t::jacobi_msg_t, &jacobi_1d_t::init> (
													   nobj, nrow, maxiter, group_proxy
													   ) ;
			   }) ;

  while( !is_work_done(group_proxy) ) {
    vt::runInEpochCollective([col_proxy] {
			       col_proxy.broadcastCollective<
				 jacobi_1d_t::blank_msg_t, &jacobi_1d_t::do_iteration
				 >() ;
			     }) ;

    vt::thePhase()->nextPhaseCollective() ;

  }
    
  vt::finalize() ;
  
  return 0 ;

} ;
