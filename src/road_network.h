#pragma once

#include <cstdint>
#include <climits>
#include <vector>

namespace road_network {

typedef uint32_t NodeID;
typedef uint32_t SubgraphID;
typedef uint32_t distance_t;

const distance_t infinity = UINT32_MAX;

//--------------------------- CutIndex ------------------------------

struct CutIndex
{
    uint64_t partition; // partition at level k is stored in k-lowest bit
    uint8_t cut_level; // level in the partition tree where vertex becomes cut-vertex (0=root, up to 63)
    uint16_t dist_index[64]; // sum of cut-sizes up to level k (indices into distances)
    std::vector<distance_t> distances; // distance to cut vertices of all levels, up to (excluding) the point where vertex becomes cut vertex
};

distance_t get_distance(const CutIndex &a, const CutIndex &b);

//--------------------------- Graph ---------------------------------

SubgraphID next_subgraph_id(bool reset = false);

struct Neighbor
{
    NodeID node;
    distance_t distance;
    Neighbor(NodeID node, distance_t distance);
};

struct Node
{
    std::vector<Neighbor> neighbors;
    // subgraph identifier
    SubgraphID subgraph_id;
    Node(SubgraphID subgraph_id);
private:
    // temporary data used by algorithms
    distance_t distance;
    NodeID inflow, outflow;
    bool is_redundant, in_partition, in_border;

    friend class Graph;
};

struct Partition
{
    std::vector<NodeID> left, right, cut;
};

class Graph
{
    // global graph
    static std::vector<Node> node_data;
    static NodeID s,t; // virtual nodes for max-flow
    // subgraph info
    std::vector<NodeID> nodes;
    SubgraphID subgraph_id;

    // check if node is contained in subgraph
    bool contains(NodeID node) const;
    // run dijkstra from node v, storing distance results in node_data
    void run_dijkstra(NodeID v);
    // run BFS from node v, storing distance results in node_data
    void run_bfs(NodeID v);
    // run BFS from t on the residual graph, storing distance results in node_data
    void run_flow_bfs();
    // returns distances from v to all subgraph nodes, with or without edge weights
    std::vector<distance_t> get_distances(NodeID v, bool weighted);
    // sort nodes by difference in distance to v and w
    void diff_sort(NodeID v, NodeID w);
    // find minimal s-t vertex cut set
    std::vector<NodeID> min_vertex_cut();
    // insert non-redundant shortcuts between border vertices
    void add_shortcuts(const std::vector<NodeID> &cut, const std::vector<NodeID> &partition);
    // recursively decompose graph and extend cut index
    void extend_cut_index(std::vector<CutIndex> &ci, float balance, uint8_t cut_level);
public:
    // create top-level graph
    Graph(uint32_t node_count = 0);
    // create subgraph
    template<typename It>
    Graph(It begin, It end) : nodes(begin, end)
    {
        subgraph_id = next_subgraph_id(false);
        assign_nodes();
    }
    // set number of nodes in global graph; global graph must currently be empty
    void resize(uint32_t node_count);
    // insert edge from v to w into global graph
    void add_edge(NodeID v, NodeID w, distance_t distance, bool add_reverse);
    // insert node into subgraph
    void add_node(NodeID v);
    // remove set of nodes from subgraph
    void remove_nodes(const std::vector<NodeID> &node_set);

    uint32_t node_count() const;
    uint32_t edge_count() const;

    // (re-)assign nodes to subgraph
    void assign_nodes();
    // returns distance between u and v in subgraph
    distance_t get_distance(NodeID v, NodeID w, bool weighted);
    // find node with maximal unweighted distance from given node
    NodeID get_furthest(NodeID v);
    // decompose graph into connected components
    void get_connected_components(std::vector<std::vector<NodeID>> &cc);
    // partition graph into balanced subgraphs using minimal cut
    void create_partition(Partition &p, float balance);
    // insert non-redundant shortcuts between border vertices
    void add_shortcuts(const Partition &p);
    // decompose graph and construct cut index
    void create_cut_index(std::vector<CutIndex> &ci, float balance);
};

} // road_network
