#include "include/node.hh"
#include "include/mesh.hh"
#include <vt/transport.h>

static constexpr std::size_t const d_no = 8 ;
static constexpr std::size_t const d_np = 8 ;
static constexpr std::size_t const d_maxiter = 50000;
static constexpr double const d_max_time = 1 ;
static constexpr double const d_CFL = 0.4 ;


int main( int argc, char** argv ){
    using namespace advect_1d ;

    size_t num_obj = d_no ;
    size_t num_points = d_np ;
    size_t max_iter = d_maxiter ;
    double max_time =  d_max_time ;
    double CFL = d_CFL ;

    double L = 1. ;
    std::string const name = argv[0];

    vt::initialize(argc, argv) ;

    vt::NodeType this_node = vt::theContext()->getNode();
    vt::NodeType num_nodes = vt::theContext()->getNumNodes();

    if (argc == 1) {
    if (this_node == 0) {
      fmt::print(
        stderr, "{}: using default arguments since none provided\n", name
      );
    }

    } else if (argc == 2) {
        num_obj = static_cast<size_t>(strtol(argv[1], nullptr, 10));
    }
    else if (argc == 3) {
        num_obj = static_cast<size_t>(strtol(argv[1], nullptr, 10));
        num_points = static_cast<size_t>(strtol(argv[2], nullptr, 10));
    }
    else if (argc == 4) {
        num_obj = static_cast<size_t>(strtol(argv[1], nullptr, 10));
        num_points = static_cast<size_t>(strtol(argv[2], nullptr, 10));
        max_iter = static_cast<size_t>(strtol(argv[3], nullptr, 10));
    }
    else {
        fmt::print(
        stderr, "usage: {} <num-objects> <num-rows-per-object> <maxiter>\n",
        name
        );
    return 1;
    }

    double const dx = L/( num_obj * num_points ) ;

    // make node collection 

    auto group_proxy = vt::theObjGroup()->makeCollective<node_t>("advect_1d") ;
 
    using base_idx_t = typename vt::Index1D::DenseIndexType ;
    auto range = vt::Index1D(static_cast<base_idx_t>(num_obj));

    auto collection_proxy = vt::makeCollection<mesh_t>("advect_1d")
    .bounds(range)
    .bulkInsert()
    .wait() ;

    vt::runInEpochCollective(
        [collection_proxy, group_proxy, num_obj, num_points, max_iter, max_time, CFL, dx] {
            collection_proxy.broadcastCollective<mesh_t::mesh_msg_t, &mesh_t::init>(
                num_points, num_obj, max_iter, max_time, CFL, dx, group_proxy
                ) ;
        }
    ) ;

    vt::runInEpochCollective(
            [collection_proxy]{
                collection_proxy.broadcastCollective<
                mesh_t::void_msg_t, &mesh_t::do_output
                > () ;
            }
        ) ;
    
    // Until termination run evolve / output 
    while( ! is_work_done(group_proxy)  ) {

        vt::runInEpochCollective(
            [collection_proxy] {
                collection_proxy.broadcastCollective<
                mesh_t::void_msg_t, &mesh_t::do_step
                > () ;
            }
        );

        vt::runInEpochCollective(
            [collection_proxy]{
                collection_proxy.broadcastCollective<
                mesh_t::void_msg_t, &mesh_t::do_output
                > () ;
            }
        ) ;

        vt::thePhase()->nextPhaseCollective() ;
    }


    vt::finalize() ;

    return 0 ;
}