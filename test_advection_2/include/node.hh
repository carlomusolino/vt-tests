#ifndef NODE_T_HH
#define NODE_T_HH

#include <vt/transport.h>

namespace advect_1d {

  struct node_t {
    bool is_finished_ = false;
    struct work_finished_msg : vt::Message {};
    // this will be called once the time of the simulation reaches max_time
    void work_finished_cback_(work_finished_msg*) { is_finished_ = true; }
    // this will be called on the root node only
    bool is_work_finished() { return is_finished_; }
  };
  using node_proxy_t = vt::objgroup::proxy::Proxy<node_t>;

  bool is_work_done( vt::objgroup::proxy::Proxy<node_t> const& this_proxy ) {
    auto const this_node = vt::theContext()->getNode() ;
    return this_proxy[this_node].invoke<
    decltype(&node_t::is_work_finished),
    &node_t::is_work_finished
    > () ;
  }

}

#endif
