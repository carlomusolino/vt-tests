#include <vt/transport.h>
#include <cmath>

static constexpr int verbose = 1 ;
static constexpr std::size_t const default_n_parts_objects = 16 ;
static constexpr std::size_t const default_n_objects = 8;
static constexpr vt::NodeType const reduce_root_node = 0;

static bool done=false ;
static const double exact_integral= M_2_PI ;


double f(double x) {
  return sin(M_PI * x ) ;
}

using reduce_msg = vt::collective::ReduceTMsg<double> ;

struct integrate_1D : vt::Collection<integrate_1D, vt::Index1D> {

private:

  size_t n_objects_ = default_n_objects ;
  size_t n_parts_per_object_ = default_n_parts_objects ;


public:

  explicit integrate_1D() :
    n_objects_(default_n_objects),
    n_parts_per_object_(default_n_parts_objects)
  { }

  struct check_result {
    void operator() (reduce_msg* msg) {
      fmt::print( " >> The integral over [0, 1] is {}\n", msg->getConstVal() ) ;
      fmt::print(
		 " >> The absolute error is {}\n",
		 std::fabs( msg->getConstVal() - exact_integral)
		 ) ;

      done = true ;
    }
  } ;

  struct init_msg: vt::CollectionMessage<integrate_1D> {
    size_t n_obj = 0;
    size_t n_interval_per_obj = 0;
    init_msg( const size_t n_obj__ , const size_t n_int) :
      n_obj(n_obj__), n_interval_per_obj(n_int)
    {}
  } ;
  
  void compute(init_msg* msg) {
    n_objects_ = msg->n_obj ;
    n_parts_per_object_ = msg->n_interval_per_obj ;

    double h = 1./(n_parts_per_object_ * n_objects_ ) ;
    double quadsum = 0.0 ;

    double a = n_parts_per_object_ * getIndex().x() * h ;

    for( int  ii=0; ii<n_parts_per_object_; ii++) {
      double x0 = a + ii * h ;
      quadsum += 0.5 * h * ( f(x0)+f(x0+h) ) ;
    }

    if ( verbose > 0 ){
      auto b = a + h * n_parts_per_object_ ; 
      fmt::print(
		 " Interval [{}, {}], on node {} & object {}, "
		 "has integral {}.\n", a, b, vt::theContext()->getNode(),
		 getIndex(), quadsum 
      );
    }


    auto proxy = this->getCollectionProxy() ;
    auto msgCB = vt::makeMessage<reduce_msg>(quadsum) ;
    auto cback = vt::theCB()->makeSend<check_result>( reduce_root_node ) ;


    proxy.reduce<vt::collective::PlusOp<double>>(msgCB.get(), cback) ;
    
  }
  
};


int main( int argc, char** argv) {

  size_t num_objs = default_n_objects ; 
  size_t n_parts = default_n_parts_objects ;

  std::string name( argv[0] ) ;

  vt::initialize(argc, argv) ;

  vt::NodeType this_node = vt::theContext()->getNode() ;
  vt::NodeType num_nodes = vt::theContext()->getNumNodes() ;

  if (argc == 1) {
    if (this_node == 0) {
      fmt::print(
        stderr, "{}: using default arguments since none provided\n", name
      );
    }
    num_objs = default_n_objects * num_nodes;
  } else {
    if (argc == 2) {
      num_objs = (size_t) strtol(argv[1], nullptr, 10);
    }
    else if (argc == 3) {
      num_objs = (size_t) strtol(argv[1], nullptr, 10);
      n_parts = (size_t) strtol(argv[2], nullptr, 10);
    }
    else {
      fmt::print(
        stderr,
        "usage: {} <num-objects> <num-interval-per-object>\n", name
      );
      return 1;
    }
  }

  using base_index_type = typename vt::Index1D::DenseIndexType ;

  auto range = vt::Index1D( static_cast<base_index_type> (num_objs) ) ;


  auto proxy = vt::makeCollection<integrate_1D> ("examples-reduce")
    .bounds(range)
    .bulkInsert()
    .wait() ;

  vt::runInEpochCollective([proxy, num_objs, n_parts]{
			     proxy.broadcastCollective<integrate_1D::init_msg, &integrate_1D::compute> (
													    num_objs, n_parts
													    );
			   }) ;

  if( vt::theContext()->getNode() == reduce_root_node ) {
    vtAssertExpr( done == true ) ;
  }

  vt::finalize() ;
  return 0 ;
}
