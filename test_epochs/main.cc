#include <vt/transport.h>

using TestMsg = vt::Message;

vt::NodeType nextNode() {
  vt::NodeType this_node = vt::theContext()->getNode();
  vt::NodeType num_nodes = vt::theContext()->getNumNodes();
  return (this_node + 1) % num_nodes;
}

static void test_handler(TestMsg* msg) {
  static int num = 3;

  vt::NodeType this_node = vt::theContext()->getNode();

  auto epoch = vt::envelopeGetEpoch(msg->env);
  fmt::print("{}: test_handler: num={}, epoch={:x}\n", this_node, num, epoch);

  num--;
  if (num > 0) {
    auto msg_send = vt::makeMessage<TestMsg>();
    vt::theMsg()->sendMsg<TestMsg, test_handler>(nextNode(), msg_send);
  }
}

int main(int argc, char** argv) {
  vt::initialize(argc, argv);

  vt::NodeType this_node = vt::theContext()->getNode();
  vt::NodeType num_nodes = vt::theContext()->getNumNodes();

  if (num_nodes == 1) {
    return vt::rerror("requires at least 2 nodes");
  }

  auto epoch = vt::theTerm()->makeEpochCollective();

  // This action will not run until all messages originating from the
  // sends are completed
  vt::theTerm()->addAction(epoch, [=]{
    fmt::print("{}: finished epoch={:x}\n", this_node, epoch);
  });

  // Message must go out of scope before finalize
  {
    auto msg = vt::makeMessage<TestMsg>();
    vt::envelopeSetEpoch(msg->env, epoch);
    vt::theMsg()->sendMsg<TestMsg, test_handler>(nextNode(), msg);
  }

  vt::theTerm()->finishedEpoch(epoch);

  vt::finalize();

  return 0;
}
