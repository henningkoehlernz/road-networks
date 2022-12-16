#include "road_network.h"

#include <vector>
#include <queue>
#include <cassert>

using namespace std;
using namespace road_network;

//--------------------------- Graph ---------------------------------

const NodeID NO_NODE = 0;

SubgraphID next_subgraph_id(bool reset = false)
{
    static SubgraphID next_id = 0;
    if (reset)
        next_id = 0;
    return next_id++;
}

Neighbor::Neighbor(NodeID node, distance_t distance) : node(node), distance(distance)
{
}

Node::Node(SubgraphID subgraph_id) : subgraph_id(subgraph_id)
{
}

vector<Node> Graph::node_data; // definition of static member

bool Graph::contains(NodeID node) const
{
    return node_data[node].subgraph_id == subgraph_id;
}

Graph::Graph(uint32_t node_count)
{
    subgraph_id = next_subgraph_id(true);
    node_data.clear();
    resize(node_count);
}

void Graph::resize(uint32_t node_count)
{
    assert(nodes.empty());
    // node numbering starts from 1
    node_data.resize(node_count + 1, Node(subgraph_id));
    nodes.reserve(node_count);
    for (NodeID node = 1; node <= node_count; node++)
        nodes.push_back(node);
}

void Graph::add_edge(NodeID v, NodeID w, distance_t distance)
{
    assert(v < node_data.size());
    assert(w < node_data.size());
    node_data[v].out.push_back(Neighbor(w, distance));
}

uint32_t Graph::node_count() const
{
    return nodes.size();
}

uint32_t Graph::edge_count() const
{
    uint32_t ecount = 0;
    for (NodeID node : nodes)
        for (Neighbor n : node_data[node].out)
            if (contains(n.node))
                ecount++;
    return ecount;
}

Graph Graph::subgraph(const vector<NodeID> &nodes)
{
    Graph g;
    g.nodes = nodes;
    g.subgraph_id = next_subgraph_id();
    // assign nodes to subgraph (inverse index)
    for (NodeID node : nodes)
        node_data[node].subgraph_id = g.subgraph_id;
    return g;
}

// helper struct to enque nodes by distance
struct SearchNode {
    distance_t distance;
    NodeID node;
    // reversed for min-heap ordering
    bool operator<(const SearchNode &other) const { return distance > other.distance; }
    SearchNode(distance_t distance, NodeID node) : distance(distance), node(node) {}
};

distance_t Graph::get_distance(NodeID v, NodeID w) const
{
    assert(contains(v));
    assert(contains(w));
    // init distances
    for (NodeID node : nodes)
        node_data[node].distance = infinity;
    node_data[v].distance = 0;
    // init queue
    priority_queue<SearchNode> q;
    q.push(SearchNode(0, v));
    // dijkstra
    while (!q.empty()) {
        SearchNode next = q.top();
        q.pop();

        for (Neighbor n : node_data[next.node].out) {
            // filter out nodes not belonging to subgraph
            if (!contains(n.node))
                continue;
            // update distance and enque
            distance_t new_dist = next.distance + n.distance;
            if (new_dist < node_data[n.node].distance) {
                node_data[n.node].distance = new_dist;
                q.push(SearchNode(new_dist, n.node));
            }
        }
    }
    return node_data[w].distance;
}

NodeID Graph::get_furthest(NodeID v) const
{
    assert(contains(v));
    // init distances - only tracks visited or not here
    for (NodeID node : nodes)
        node_data[node].distance = infinity;
    node_data[v].distance = 0;
    // init queue
    queue<NodeID> q;
    q.push(v);
    // BFS
    NodeID last_node = v;
    while (!q.empty()) {
        last_node = q.front();
        q.pop();

        for (Neighbor n : node_data[last_node].out) {
            // filter out nodes not belonging to subgraph or already visited
            if (!contains(n.node) || node_data[n.node].distance == 0)
                continue;
            // update distance and enque
            node_data[n.node].distance = 0;
            q.push(n.node);
        }
    }
    return last_node;
}

//--------------------------- CutIndex ------------------------------

const uint8_t NON_CUT_VERTEX = 64; // for use with cut_level

// get distance when one vertex is a cut vertex for a subgraph containing both
distance_t direct_dist(const CutIndex &a, const CutIndex &b)
{
    uint16_t a_index = a.distances.size();
    uint16_t b_index = b.distances.size();
    return a_index < b_index ? b.distances[a_index]
        : a_index > b_index ? a.distances[b_index]
        : 0;
}

distance_t get_dist(const CutIndex &a, const CutIndex &b)
{
    // same leaf node, or one vertex is cut vertex
    if (a.partition == b.partition)
        return direct_dist(a, b);
    // find lowest level at which partitions differ
    uint64_t pdiff = a.partition ^ b.partition;
    // partition level used for comparison (upper bound initially)
    int pindex = min(a.cut_level, b.cut_level);
    int diff_level = __builtin_ctz(pdiff); // count trailing zeros
    // one vertex is cut vertex
    if (pindex <= diff_level)
        return direct_dist(a,b);
    pindex = diff_level;
    // compute iterator range
    const distance_t* a_end = &a.distances[0] + a.dist_index[pindex];
    const uint16_t offset = pindex ? a.dist_index[pindex - 1] : 0; // same for a and b
    const distance_t* a_ptr = &a.distances[0] + offset;
    const distance_t* b_ptr = &b.distances[0] + offset;
    // find min 2-hop distance within partition
    distance_t min_dist = infinity;
    while (a_ptr != a_end)
    {
        distance_t dist = *a_ptr + *b_ptr;
        if (dist < min_dist)
            min_dist = dist;
        a_ptr++;
        b_ptr++;
    }
    return min_dist;
}
