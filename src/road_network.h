#pragma once

#include <cstdint>
#include <climits>
#include <vector>

namespace road_network {

typedef uint32_t NodeID;
typedef uint32_t SubgraphID;
typedef uint32_t distance_t;

const distance_t infinity = UINT32_MAX;

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

    // check if node is contained in subgraph
    bool contains(NodeID node) const;
    // run dijkstra from node v, storing distance results in node_data
    void run_dijkstra(NodeID v);
    // run BFS from node v, storing distance results in node_data
    void run_bfs(NodeID v);
    // returns distances from v to all subgraph nodes, with or without edge weights
    std::vector<distance_t> get_distances(NodeID v, bool weighted);
    // sorts nodes by difference in distance to v and w
    void diff_sort(NodeID v, NodeID w);
public:
    Graph(uint32_t node_count = 0);
    // set number of nodes; must not have been set in constructor
    void resize(uint32_t node_count);
    // insert edge from v to w
    void add_edge(NodeID v, NodeID w, distance_t distance);

    uint32_t node_count() const;
    uint32_t edge_count() const;

    // create subgraph
    Graph subgraph(const std::vector<NodeID> &nodes);
    // (re-)assign nodes to subgraph
    void assign_nodes();
    // returns distance between u and v in subgraph
    distance_t get_distance(NodeID v, NodeID w, bool weighted);
    // find node with maximal unweighted distance from given node
    NodeID get_furthest(NodeID v);
};

struct CutIndex
{
    uint64_t partition; // partition at level k is stored in k-lowest bit
    uint8_t cut_level; // level in the partition tree where vertex becomes cut-vertex (0=root, up to 63)
    uint16_t dist_index[64]; // sum of cut-sizes up to level k (indices into distances)
    std::vector<distance_t> distances; // distance to cut vertices of all levels, up to (excluding) the point where vertex becomes cut vertex
};

distance_t get_distance(const CutIndex &a, const CutIndex &b);

} // road_network
