#pragma once

#define NDEBUG
#define NPROFILE
#define CHECK_CONSISTENT //assert(is_consistent())
// algorithm config
//#define NO_SHORTCUTS // turns off shortcut computation, resulting in smaller indexes but slower local queries
#define CONTRACT2D // contract nodes of degree 2 for faster cut computation
#define PRUNING // enables tail-pruning, resulting in smaller indexes but increased construction time

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

const distance_t infinity = UINT32_MAX >> 1;

struct Neighbor;
class Graph;

//--------------------------- CutIndex ------------------------------

struct CutIndex
{
    uint64_t partition; // partition at level k is stored in k-lowest bit
    uint8_t cut_level; // level in the partition tree where vertex becomes cut-vertex (0=root, up to 58)
    std::vector<uint16_t> dist_index; // sum of cut-sizes up to level k (indices into distances)
    std::vector<distance_t> distances; // distance to cut vertices of all levels, up to (excluding) the point where vertex becomes cut vertex
#ifdef PRUNING
    // track number of labels that could be or are pruned
    size_t pruning_2hop, pruning_3hop, pruning_tail;
    // prune maximum tail of labels in current cut; distance labels must already containing pruning flag
    void prune_tail();
#endif

    CutIndex();
    bool is_consistent(bool partial=false) const;
    bool empty() const;
};

std::ostream& operator<<(std::ostream& os, const CutIndex &ci);

// helper functions for manipulating partition bitvectors
namespace PBV
{
    // construct partition bitvector from bit pattern and length
    uint64_t from(uint64_t bits, uint16_t length);
    // split partition bitvector into components
    uint64_t partition(uint64_t bv);
    uint16_t cut_level(uint64_t bv);
    // compute cut level of least common ancestor of given bitvectors
    uint16_t lca_level(uint64_t bv1, uint64_t bv2);
    // compute bitvector for least common ancestor of given bitvectors
    uint64_t lca(uint64_t bv1, uint64_t bv2);
    // check whether node is an ancestor of another, based on their bitvectors
    bool is_ancestor(uint64_t bv_ancestor, uint64_t bv_descendant);
}

class FlatCutIndex
{
    char* data; // stores partition bitvector, dist_index and distances
public:
    FlatCutIndex();
    FlatCutIndex(const CutIndex &ci);

    bool operator==(FlatCutIndex other) const;

    // return pointers to partition bitvector, dist_index and distances array
    uint64_t* partition_bitvector();
    const uint64_t* partition_bitvector() const;
    uint16_t* dist_index();
    const uint16_t* dist_index() const;
    distance_t* distances();
    const distance_t* distances() const;
    // split partition_bitvector into components
    uint64_t partition() const;
    uint16_t cut_level() const;

    // number of bytes allocated for index data
    size_t size() const;
    // number of labels
    size_t label_count() const;
    // number of labels at given cut level
    size_t cut_size(size_t cl) const;
    // number of labels at lowest cut level
    size_t bottom_cut_size() const;
    // returns whether index data has been allocated
    bool empty() const;

    // start of distance labels for given cut level
    const distance_t* cl_begin(size_t cl) const;
    // end of distance labels for given cut level
    const distance_t* cl_end(size_t cl) const;

    // returns labels in list-of-list format
    std::vector<std::vector<distance_t>> unflatten() const;

    friend class ContractionIndex;
};

std::ostream& operator<<(std::ostream& os, const FlatCutIndex &ci);

struct ContractionLabel
{
    FlatCutIndex cut_index;
    distance_t distance_offset; // distance to node owning the labels
    NodeID parent; // parent in tree rooted at label-owning node

    ContractionLabel();
    // index size in bytes
    size_t size() const;
};

std::ostream& operator<<(std::ostream& os, const ContractionLabel &ci);

class ContractionIndex
{
    std::vector<ContractionLabel> labels;

    static distance_t get_cut_level_distance(FlatCutIndex a, FlatCutIndex b, size_t cut_level);
    static distance_t get_distance(FlatCutIndex a, FlatCutIndex b);
    static size_t get_cut_level_hoplinks(FlatCutIndex a, FlatCutIndex b, size_t cut_level);
    static size_t get_hoplinks(FlatCutIndex a, FlatCutIndex b);
public:
    // populate from ci and closest, draining ci in the process
    ContractionIndex(std::vector<CutIndex> &ci, std::vector<Neighbor> &closest);
    // populate from binary source
    ContractionIndex(std::istream& is);
    // wrapper when not contracting
    explicit ContractionIndex(std::vector<CutIndex> &ci);
    ~ContractionIndex();

    // compute distance between v and w
    distance_t get_distance(NodeID v, NodeID w) const;
    // verify correctness of distance computed via index for a particular query
    bool check_query(std::pair<NodeID,NodeID> query, Graph &g) const;

    // check whether node had its labels pruned during contraction
    bool is_contracted(NodeID node) const;
    size_t uncontracted_count() const;
    // check whether node lies in contracted partition identified by bitvector
    bool in_partition_subgraph(NodeID node, uint64_t partition_bitvector) const;

    // compute number of hoplinks examined during distance computation
    size_t get_hoplinks(NodeID v, NodeID w) const;
    double avg_hoplinks(const std::vector<std::pair<NodeID,NodeID>> &queries) const;
    // index size in bytes
    size_t size() const;
    double avg_cut_size() const;
    size_t max_cut_size() const;
    size_t height() const;
    size_t max_label_count() const;
    size_t label_count() const;
    // number of non-empty cuts
    size_t non_empty_cuts() const;

    // generate random query
    std::pair<NodeID,NodeID> random_query() const;
    // write index in binary format
    void write(std::ostream& os) const;
    // write index in json format
    void write_json(std::ostream& os) const;
};

//--------------------------- Graph ---------------------------------

SubgraphID next_subgraph_id(bool reset = false);

struct Neighbor
{
    NodeID node;
    distance_t distance;
    Neighbor(NodeID node, distance_t distance);
    bool operator<(const Neighbor &other) const;
};

std::ostream& operator<<(std::ostream& os, const Neighbor &n);

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
#ifdef CONTRACT2D
    std::vector<uint32_t> deg2path_ids;
#endif
    NodeID inflow, outflow;
    uint16_t landmark_level;

    friend class Graph;
};

std::ostream& operator<<(std::ostream& os, const Node &n);

// multi-threading requires thread-local data for s & t nodes
class MultiThreadNodeData : public std::vector<Node>
{
    thread_local static Node s_data, t_data;
public:
    Node& operator[](size_type pos);
    const Node& operator[](size_type pos) const;
    void normalize();
};

struct Partition
{
    std::vector<NodeID> left, right, cut;
    // rates quality of partition (cutsize + balance)
    double rating() const;
};

std::ostream& operator<<(std::ostream& os, const Partition &p);

struct Edge
{
    NodeID a, b;
    distance_t d;
    Edge(NodeID a, NodeID b, distance_t d);
    bool operator<(Edge other) const;
};

// helper structure for pre-partitioning
struct DiffData
{
    NodeID node;
    distance_t dist_a, dist_b;
    int32_t diff() const;
    distance_t min() const;

    DiffData(NodeID node, distance_t dist_a, distance_t dist_b);
    // comparison function for easy sorting by diff values
    static bool cmp_diff(DiffData x, DiffData y);

    friend std::ostream& operator<<(std::ostream& os, const DiffData &dd);
};

/**
 * full graph information (edges and weights) is only stored once, as static data; graph instances describe induced subgraphs, storing only a list of nodes;
 * this approach speeds up creation of subgraphs, and saves memory, but complicates usage;
 * when traversing a subgraph by visiting neighbors (e.g. distance computation), we need to check for each neighbor whether it lies in the subgraph;
 * to speed up this check, each node stores the subgraph it currently belongs to - this information must be carefully maintained;
 * it also causes conflicts between overlapping subgraphs during parallel processing, hence such cases must be avoided
 */
class Graph
{
    // global graph
#ifdef MULTI_THREAD
    static MultiThreadNodeData node_data;
    static size_t thread_threshold; // minimum subgraph size for which processing will be split across multiple threads
#else
    static std::vector<Node> node_data;
#endif
#ifdef CONTRACT2D
    static std::vector<std::vector<NodeID>> deg2paths; // contracted paths of degree two nodes
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
    // remove set of nodes from subgraph; node_set must be sorted
    void remove_nodes(const std::vector<NodeID> &node_set);
    // return single neighbor of degree one node, or NO_NODE otherwise
    Neighbor single_neighbor(NodeID v) const;
    // return neighbors of degree two node, or pair of NO_NODE if degree > 2
    std::pair<NodeID,NodeID> pair_of_neighbors(NodeID v) const;
    // find Neighbor structure in v pointing to w
    Neighbor& get_neighbor(NodeID v, NodeID w);

    // run dijkstra from node v, storing distance results in node_data
    void run_dijkstra(NodeID v);
    // run dijkstra from node v, in subgraph excluding lower-level landmarks
    void run_dijkstra_llsub(NodeID v);
    // stores whether all shortest paths bypass other landmarks in lowest distance bit
    void run_dijkstra_ll(NodeID v);
#ifdef MULTI_THREAD_DISTANCES
    // run dijkstra from multiple nodes in parallel
    void run_dijkstra_par(const std::vector<NodeID> &vertices);
    // run dijkstra from multiple nodes in parallel, in subgraph excluding lower-level landmarks
    void run_dijkstra_llsub_par(const std::vector<NodeID> &vertices);
    // stores whether all shortest paths bypass other landmarks in lowest distance bit
    void run_dijkstra_ll_par(const std::vector<NodeID> &vertices);
#endif
    // run BFS from node v, storing distance results in node_data
    void run_bfs(NodeID v);
    // run BFS from s (forward) or t (backward) on the residual graph, storing distance results in node_data
    void run_flow_bfs_from_s();
    void run_flow_bfs_from_t();

#ifdef CONTRACT2D
    // contract paths of nodes of degree two
    void contract_deg2paths();
#endif

    // find node with maximal distance from given node
    std::pair<NodeID,distance_t> get_furthest(NodeID v, bool weighted);
    // find pair of nodes with maximal distance
    Edge get_furthest_pair(bool weighted);
    // get distances of nodes to a and b; pre-computed indicates that node_data already holds distances to a
    void get_diff_data(std::vector<DiffData> &diff, NodeID a, NodeID b, bool weighted, bool pre_computed = false);
    // find one or more minimal s-t vertex cut sets
    void min_vertex_cuts(std::vector<std::vector<NodeID>> &cuts);
    // find cut from given rough partition
    void rough_partition_to_cuts(std::vector<std::vector<NodeID>> &cuts, const Partition &p);
    // compute left/right partitions based on given cut
    void complete_partition(Partition &p);
    // insert non-redundant shortcuts between border vertices
    void add_shortcuts(const std::vector<NodeID> &cut, const std::vector<CutIndex> &ci);
    // order cut vertices in order of pruning potential (latter nodes can be pruned better)
    void sort_cut_for_pruning(std::vector<NodeID> &cut, std::vector<CutIndex> &ci);
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
    // turn progress tracking on/off
    static void show_progress(bool state);
    // number of nodes in the top-level graph
    static size_t super_node_count();

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
    // reset graph to contain all nodes in global graph
    void reset();

    size_t node_count() const;
    size_t edge_count() const;
    size_t degree(NodeID v) const;
    // approximate diameter
    distance_t diameter(bool weighted);
    // returns list of nodes
    const std::vector<NodeID>& get_nodes() const;
    // returns list of all edges (one per undirected edge)
    void get_edges(std::vector<Edge> &edges) const;

    // returns distance between u and v in subgraph
    distance_t get_distance(NodeID v, NodeID w, bool weighted);
    // decompose graph into connected components
    void get_connected_components(std::vector<std::vector<NodeID>> &cc);
    // computed rough partition with wide separator, returned in p; returns if rough partition is already a partition
    bool get_rough_partition(Partition &p, double balance, bool disconnected);
    // partition graph into balanced subgraphs using minimal cut
    void create_partition(Partition &p, double balance);
    // decompose graph and construct cut index; returns number of shortcuts used
    size_t create_cut_index(std::vector<CutIndex> &ci, double balance);
    // returns edges that don't affect distances between nodes
    void get_redundant_edges(std::vector<Edge> &edges);
    // repeatedly remove nodes of degree 1, populating closest[removed] with next node on path to closest unremoved node
    void contract(std::vector<Neighbor> &closest);

    // generate random node
    NodeID random_node() const;
    // generate random pair of nodes through random walk (0 = fully random)
    std::pair<NodeID,NodeID> random_pair(size_t steps = 0) const;
    // generate batch of random node pairs, filtered into buckets by distance (as for H2H/P2H)
    void random_pairs(std::vector<std::vector<std::pair<NodeID,NodeID>>> &buckets, distance_t min_dist, size_t bucket_size, const ContractionIndex &ci);
    // randomize order of nodes and neighbors
    void randomize();

    friend std::ostream& operator<<(std::ostream& os, const Graph &g);
    friend MultiThreadNodeData;
};

// print graph in DIMACS format
void print_graph(const Graph &g, std::ostream &os);
// read graph in DIMACS format
void read_graph(Graph &g, std::istream &in);

} // road_network
