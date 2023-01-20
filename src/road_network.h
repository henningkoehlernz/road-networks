#pragma once

#define NDEBUG
#define NPROFILE
#define CHECK_CONSISTENT //assert(is_consistent())
// algorithm config
//#define NO_SHORTCUTS
#ifdef NO_SHORTCUTS
    //#define CUT_BOUNDS
#endif

// use multi-threading for index construction
#define MULTI_THREAD 32 // determines threshold for multi-threading
#ifdef MULTI_THREAD
    #define MULTI_THREAD_DISTANCES 4 // number of parallel threads for label & shortcut computation
#endif

#include <cstdint>
#include <climits>
#include <vector>
#include <ostream>
#include <cassert>

namespace road_network {

typedef uint32_t NodeID;
typedef uint32_t SubgraphID;
typedef uint32_t distance_t;

const distance_t infinity = UINT32_MAX;

enum class DistanceMeasure { unweighted, weighted, mixed };

//--------------------------- CutIndex ------------------------------

struct CutIndex
{
    uint64_t partition; // partition at level k is stored in k-lowest bit
    uint8_t cut_level; // level in the partition tree where vertex becomes cut-vertex (0=root, up to 63)
    std::vector<uint16_t> dist_index; // sum of cut-sizes up to level k (indices into distances)
    std::vector<distance_t> distances; // distance to cut vertices of all levels, up to (excluding) the point where vertex becomes cut vertex
#ifdef CUT_BOUNDS
    std::vector<distance_t> cut_bounds; // min(distances[0..(index+1)*constant])
#endif
};

// compute distance between two vertices using their cut index data
distance_t get_distance(const CutIndex &a, const CutIndex &b);
// number of distance pairs used for distance computation
size_t get_hoplinks(const CutIndex &a, const CutIndex &b);
double avg_hoplinks(const std::vector<CutIndex> &ci, const std::vector<std::pair<NodeID,NodeID>> &queries);
// sums up total number of labels in index
size_t label_count(const std::vector<CutIndex> &ci);
// compute size of cut index in bytes
size_t index_size(const std::vector<CutIndex> &ci);
// average cut size, weighted by partition size
double avg_cut_size(const std::vector<CutIndex> &ci);

size_t max_cut_size(const std::vector<CutIndex> &ci);
// height of partitioning tree
size_t index_height(const std::vector<CutIndex> &ci);

std::ostream& operator<<(std::ostream& os, const CutIndex &ci);

//--------------------------- Graph ---------------------------------

SubgraphID next_subgraph_id(bool reset = false);

struct Neighbor
{
    NodeID node;
    distance_t distance;
    Neighbor(NodeID node, distance_t distance);
    bool operator<(const Neighbor &other) const;
};

struct Node
{
    std::vector<Neighbor> neighbors;
    // subgraph identifier
    SubgraphID subgraph_id;
    Node(SubgraphID subgraph_id);
private:
    // temporary data used by algorithms
    distance_t distance, outcopy_distance;
#ifdef MULTI_THREAD_DISTANCES
    distance_t distances[MULTI_THREAD_DISTANCES];
#endif
    NodeID inflow, outflow;

    friend class Graph;
};

// multi-threading requires thread-local data for s & t nodes
class MultiThreadNodeData : public std::vector<Node>
{
    thread_local static Node s_data, t_data;
public:
    Node& operator[](size_type pos);
    const Node& operator[](size_type pos) const;
};

struct Partition
{
    std::vector<NodeID> left, right, cut;
    // rates quality of partition (cutsize + balance)
    double rating() const;
    friend std::ostream& operator<<(std::ostream& os, const Partition &p);
};

struct Edge
{
    NodeID a, b;
    distance_t d;
    Edge(NodeID a, NodeID b, distance_t d);
    bool operator<(Edge other) const;
};

class Graph
{
    // global graph
#ifdef MULTI_THREAD
    static MultiThreadNodeData node_data;
    static size_t thread_threshold;
#else
    static std::vector<Node> node_data;
#endif
    static NodeID s,t; // virtual nodes for max-flow
    // subgraph info
    std::vector<NodeID> nodes;
    SubgraphID subgraph_id;

    // create subgraph
    template<typename It>
    Graph(It begin, It end) : nodes(begin, end)
    {
        subgraph_id = next_subgraph_id(false);
        assign_nodes();
        CHECK_CONSISTENT;
    }

    // (re-)assign nodes to subgraph
    void assign_nodes();
    // check if node is contained in subgraph
    bool contains(NodeID node) const;
    // insert node into subgraph
    void add_node(NodeID v);
    // remove set of nodes from subgraph
    void remove_nodes(const std::vector<NodeID> &node_set);

    // run dijkstra from node v, storing distance results in node_data
    void run_dijkstra(NodeID v);
#ifdef MULTI_THREAD_DISTANCES
    // run dijkstra from multiple nodes in parallel
    void run_dijkstra_par(const std::vector<NodeID> &vertices);
#endif
    // run BFS from node v, storing distance results in node_data
    void run_bfs(NodeID v);
    // run BFS from t on the residual graph, storing distance results in node_data
    void run_flow_bfs();

    // find node with maximal distance from given node
    std::pair<NodeID,distance_t> get_furthest(NodeID v, DistanceMeasure dm);
    // find pair of nodes with maximal distance
    Edge get_furthest_pair(DistanceMeasure dm);
    // sort nodes by difference in distance to v and w
    void diff_sort(NodeID v, NodeID w, bool precomputed);
    // find minimal s-t vertex cut set
    std::vector<NodeID> min_vertex_cut();
    // insert non-redundant shortcuts between border vertices
    void add_shortcuts(const std::vector<NodeID> &cut, const std::vector<CutIndex> &ci);
    // recursively extend cut index onto given partition, using given cut
    static void extend_on_partition(std::vector<CutIndex> &ci, double balance, uint8_t cut_level, const std::vector<NodeID> &p, const std::vector<NodeID> &cut);
    // recursively decompose graph and extend cut index
    void extend_cut_index(std::vector<CutIndex> &ci, double balance, uint8_t cut_level);

    // check if subgraph_id assignment is consistent with nodes
    bool is_consistent() const;
    // check if neighorhood relationship (with distances) is symmetrical
    bool is_undirected() const;
    // return internal node distances as vector
    std::vector<std::pair<distance_t,distance_t>> distances() const;
    // return internal flow values as vector
    std::vector<std::pair<NodeID,NodeID>> flow() const;
public:
    // create top-level graph
    Graph(size_t node_count = 0);
    Graph(size_t node_count, const std::vector<Edge> &edges);
    // set number of nodes in global graph; global graph must currently be empty
    void resize(size_t node_count);
    // insert edge from v to w into global graph
    void add_edge(NodeID v, NodeID w, distance_t distance, bool add_reverse);
    // remove edge between v and w from global graph
    void remove_edge(NodeID v, NodeID w);
    // remove isolated nodes from subgraph
    void remove_isolated();
    // undo changes made during subgraph construction
    void reset();

    size_t node_count() const;
    size_t edge_count() const;
    size_t degree(NodeID v) const;
    // number of nodes in the top-level graph
    static size_t super_node_count();
    // approximate diameter
    distance_t diameter(bool weighted);
    // returns list of all edges (one per undirected edge)
    void get_edges(std::vector<Edge> &edges) const;

    // returns distance between u and v in subgraph
    distance_t get_distance(NodeID v, NodeID w, bool weighted);
    // decompose graph into connected components
    void get_connected_components(std::vector<std::vector<NodeID>> &cc);
    // partition graph into balanced subgraphs using minimal cut
    void create_partition(Partition &p, double balance);
    // decompose graph and construct cut index
    void create_cut_index(std::vector<CutIndex> &ci, double balance);
    // returns edges that don't affect distances between nodes
    void get_redundant_edges(std::vector<Edge> &edges, const std::vector<CutIndex> &ci) const;
    void get_redundant_edges(std::vector<Edge> &edges);

    // generate random node
    NodeID random_node() const;
    // generate random pair of nodes through random walk (0 = fully random)
    std::pair<NodeID,NodeID> random_pair(size_t steps = 0) const;
    // generate batch of random node pairs, filtered into buckets by distance (as for H2H/P2H)
    void random_pairs(std::vector<std::vector<std::pair<NodeID,NodeID>>> &buckets, distance_t min_dist, size_t bucket_size, const std::vector<CutIndex> &ci);
    // randomize order of nodes and neighbors
    void randomize();

    // verify correctness of distance computed via cut index for a particular query
    bool check_cut_index(const std::vector<CutIndex> &ci, std::pair<NodeID,NodeID> query);

    friend std::ostream& operator<<(std::ostream& os, const Neighbor &n);
    friend std::ostream& operator<<(std::ostream& os, const Node &n);
    friend std::ostream& operator<<(std::ostream& os, const Graph &g);
    friend MultiThreadNodeData;
};

// print graph in DIMACS format
void print_graph(const Graph &g, std::ostream &os);
// read graph in DIMACS format
void read_graph(Graph &g, std::istream &in);

} // road_network
