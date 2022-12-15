#pragma once

#include <cstdint>
#include <climits>
#include <vector>

namespace road_network {

typedef uint32_t NodeID;
typedef uint32_t SubgraphID;
typedef uint32_t distance_t;

const distance_t infinity = UINT32_MAX;
const NodeID NO_NODE = 0;

struct Neighbor
{
    NodeID node;
    distance_t distance;
    Neighbor(NodeID node, distance_t distance);
};

struct Node
{
    // outgoing neighbors
    std::vector<Neighbor> out;
    // subgraph identifier
    SubgraphID subgraph_id;
    Node(SubgraphID subgraph_id);
private:
    // temporary data used by algorithms
    distance_t distance;
    NodeID flow_in, flow_out;

    friend class Graph;
};

class Graph
{
    // global graph
    static std::vector<Node> node_data;
    // subgraph info
    std::vector<NodeID> nodes;
    SubgraphID subgraph_id;
public:
    Graph(uint32_t node_count = 0);
    void add_edge(NodeID v, NodeID w, distance_t distance);
    // create subgraph
    Graph subgraph(const std::vector<NodeID> &nodes);
    // sample implementation of Dijktra over subgraph
    distance_t get_distance(NodeID u, NodeID v) const;
};

} // road_network
