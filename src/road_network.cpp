#include "road_network.h"
#include "util.h"

#include <vector>
#include <queue>
#include <cassert>
#include <algorithm>

using namespace std;
using namespace road_network;

//--------------------------- Graph ---------------------------------

const NodeID NO_NODE = 0;

SubgraphID road_network::next_subgraph_id(bool reset)
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

// definition of static members
vector<Node> Graph::node_data;
NodeID Graph::s, Graph::t;

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
    // node numbering starts from 1, and we reserve two additional nodes for s & t
    node_data.resize(node_count + 3, Node(subgraph_id));
    s = node_count + 1;
    t = node_count + 2;
    nodes.reserve(node_count);
    for (NodeID node = 1; node <= node_count; node++)
        nodes.push_back(node);
}

void Graph::add_edge(NodeID v, NodeID w, distance_t distance, bool add_reverse)
{
    assert(v < node_data.size());
    assert(w < node_data.size());
    node_data[v].neighbors.push_back(Neighbor(w, distance));
    if (add_reverse)
        node_data[w].neighbors.push_back(Neighbor(v, distance));
}

void Graph::add_node(NodeID v)
{
    assert(v < node_data.size());
    nodes.push_back(v);
    node_data[v].subgraph_id = subgraph_id;
}

void Graph::remove_nodes(const vector<NodeID> &node_set)
{
    util::remove_set(nodes, node_set);
}

uint32_t Graph::node_count() const
{
    return nodes.size();
}

uint32_t Graph::edge_count() const
{
    uint32_t ecount = 0;
    for (NodeID node : nodes)
        for (Neighbor n : node_data[node].neighbors)
            if (contains(n.node))
                ecount++;
    return ecount;
}

void Graph::assign_nodes()
{
    for (NodeID node : nodes)
        node_data[node].subgraph_id = subgraph_id;
}

//--------------------------- Graph algorithms ----------------------

// helper struct to enque nodes by distance
struct SearchNode
{
    distance_t distance;
    NodeID node;
    // reversed for min-heap ordering
    bool operator<(const SearchNode &other) const { return distance > other.distance; }
    SearchNode(distance_t distance, NodeID node) : distance(distance), node(node) {}
};

void Graph::run_dijkstra(NodeID v)
{
    assert(contains(v));
    // init distances
    for (NodeID node : nodes)
        node_data[node].distance = infinity;
    node_data[v].distance = 0;
    // init queue
    priority_queue<SearchNode> q;
    q.push(SearchNode(0, v));
    // dijkstra
    while (!q.empty())
    {
        SearchNode next = q.top();
        q.pop();

        for (Neighbor n : node_data[next.node].neighbors)
        {
            // filter neighbors nodes not belonging to subgraph
            if (!contains(n.node))
                continue;
            // update distance and enque
            distance_t new_dist = next.distance + n.distance;
            if (new_dist < node_data[n.node].distance)
            {
                node_data[n.node].distance = new_dist;
                q.push(SearchNode(new_dist, n.node));
            }
        }
    }
}

void Graph::run_bfs(NodeID v)
{
    assert(contains(v));
    // init distances
    for (NodeID node : nodes)
        node_data[node].distance = infinity;
    node_data[v].distance = 0;
    // init queue
    queue<NodeID> q;
    q.push(v);
    // BFS
    while (!q.empty())
    {
        NodeID next = q.front();
        q.pop();

        distance_t new_dist = node_data[next].distance + 1;
        for (Neighbor n : node_data[next].neighbors)
        {
            // filter neighbors nodes not belonging to subgraph or already visited
            if (contains(n.node) && node_data[n.node].distance == infinity)
            {
                // update distance and enque
                node_data[n.node].distance = new_dist;
                q.push(n.node);
            }
        }
    }
}

// node in flow graph which splits nodes into incoming and outgoing copies
struct FlowNode
{
    NodeID node;
    bool outcopy; // outgoing copy of node?
    FlowNode(NodeID node, bool outcopy) : node(node), outcopy(outcopy) {}
};

void Graph::run_flow_bfs()
{
    assert(contains(s) && contains(t));
    // init distances
    for (NodeID node : nodes)
        node_data[node].distance = infinity;
    node_data[t].distance = 0;
    // init queue
    queue<FlowNode> q;
    q.push(FlowNode(t, false));
    // BFS
    while (!q.empty())
    {
        FlowNode next = q.front();
        q.pop();

        distance_t new_dist = node_data[next.node].distance + 1;
        NodeID outflow = node_data[next.node].outflow;
        // special treatment is needed for node with flow through it
        if (outflow != NO_NODE && next.outcopy)
        {
            // outflow is only valid neighbor
            if (node_data[outflow].distance == infinity)
            {
                node_data[outflow].distance = new_dist;
                q.push(FlowNode(outflow, false));
            }
        }
        // when arriving at the incoming copy of flow node, all neighbors except inflow are valid
        // inflow must have been already visited in this case, so checking all neighbors is fine
        else for (Neighbor n : node_data[next.node].neighbors)
        {
            // filter neighbors nodes not belonging to subgraph or already visited
            if (contains(n.node) && node_data[n.node].distance == infinity)
            {
                // update distance and enque
                node_data[n.node].distance = new_dist;
                q.push(FlowNode(n.node, n.node != outflow));
            }
        }
    }
}

distance_t Graph::get_distance(NodeID v, NodeID w, bool weighted)
{
    weighted ? run_dijkstra(v) : run_bfs(v);
    assert(contains(w));
    return node_data[w].distance;
}

NodeID Graph::get_furthest(NodeID v)
{
    run_bfs(v);
    NodeID furthest = v;
    for (NodeID node : nodes)
        if (node_data[node].distance > node_data[furthest].distance)
            furthest = node;
    return furthest;
}

vector<distance_t> Graph::get_distances(NodeID v, bool weighted)
{
    weighted ? run_dijkstra(v) : run_bfs(v);
    vector<distance_t> d(nodes.size());
    for (NodeID node : nodes)
        d.push_back(node_data[node].distance);
    return d;
}

void Graph::diff_sort(NodeID v, NodeID w)
{
    // compute distance difference
    uint32_t node_count = nodes.size();
    vector<pair<int32_t,NodeID>> diff(node_count);
    run_bfs(v);
    for (NodeID node : nodes)
        diff.push_back(pair(node_data[node].distance, node));
    run_bfs(w);
    for (uint32_t i = 0; i < node_count; i++)
        diff[i].first -= node_data[nodes[i]].distance;
    // sort & replace
    std::sort(diff.begin(), diff.end());
    for (uint32_t i = 0; i < node_count; i++)
        nodes[i] = diff[i].second;
}

vector<NodeID> Graph::min_vertex_cut()
{
    assert(contains(s) && contains(t));
    // set flow to empty
    for (NodeID node : nodes)
        node_data[node].inflow = node_data[node].outflow = NO_NODE;
    // find max s-t flow using Dinitz' algorithm
    while (true)
    {
        // construct BFS tree from t
        run_flow_bfs();
        distance_t s_distance = node_data[s].distance;
        if (s_distance == infinity)
            break;
        // run DFS from s along inverse BFS tree edges
        vector<NodeID> path, stack;
        // iterating over neighbors of s directly simplifies stack cleanup after new s-t path is found
        for (Neighbor sn : node_data[s].neighbors)
        {
            if (!contains(sn.node) || node_data[sn.node].distance != s_distance - 1)
                continue;
            stack.push_back(sn.node);
            while (!stack.empty())
            {
                NodeID node = stack.back();
                stack.pop_back();
                // clean up path (back tracking)
                while (node_data[path.back()].distance <= node_data[node].distance)
                    path.pop_back();
                // increase flow when s-t path is found
                if (node == t)
                {
                    assert(node_data[path.front()].inflow == NO_NODE);
                    node_data[path.front()].inflow = s;
                    for (size_t path_pos = 1; path_pos < path.size(); path_pos++)
                    {
                        NodeID from = path[path_pos - 1];
                        NodeID to = path[path_pos];
                        // we might be reverting existing flow
                        // from.inflow has already been changed => check outflow
                        if (node_data[from].outflow != NO_NODE)
                        {
                            assert(node_data[to].outflow == from);
                            node_data[to].outflow = NO_NODE;
                        }
                        else
                        {
                            node_data[from].outflow = to;
                            node_data[to].inflow = from;
                        }
                    }
                    assert(node_data[path.back()].outflow == NO_NODE);
                    node_data[path.back()].outflow = t;
                    // ensure vertices in path are not re-visited during current DFS iteration
                    for (NodeID path_node : path)
                        node_data[path_node].distance = infinity;
                    // skip to next neighbor of s
                    stack.clear();
                    path.clear();
                    break;
                }
                // continue DFS from node
                path.push_back(node);
                distance_t node_distance = node_data[node].distance;
                for (Neighbor n : node_data[s].neighbors)
                {
                    if (contains(n.node) && node_data[n.node].distance == node_distance - 1)
                        stack.push_back(n.node);
                }
            }
        }
    }
    // find min cut
}

void Graph::create_partition(Partition &p, float balance)
{
    assert(nodes.size() > 1);
    // find two extreme points
    NodeID a = nodes[0];
    NodeID b = get_furthest(a);
    a = get_furthest(b);
    // create pre-partition
    diff_sort(a,b);
    uint32_t max_left = 1 + nodes.size() * balance;
    uint32_t min_right = nodes.size() * (1 - balance);
    assert(max_left <= min_right);
    auto max_left_it = nodes.begin() + max_left;
    auto min_right_it = nodes.begin() + min_right;
    Graph left(nodes.begin(), max_left_it);
    Graph center(max_left_it, min_right_it);
    Graph right(min_right_it, nodes.end());
    // construct s-t flow graph
    center.add_node(s);
    center.add_node(t);
    // handle corner case of edges between left and right partition
    // do this first as it can eliminate other s/t neighbors
    vector<NodeID> s_neighbors, t_neighbors;
    for (NodeID node : left.nodes)
        for (Neighbor n : node_data[node].neighbors)
            if (right.contains(n.node))
            {
                s_neighbors.push_back(node);
                t_neighbors.push_back(n.node);
            }
    util::make_set(s_neighbors);
    util::make_set(t_neighbors);
    // update pre-partition
    left.remove_nodes(s_neighbors);
    for (NodeID node : s_neighbors)
        center.add_node(node);
    right.remove_nodes(t_neighbors);
    for (NodeID node : t_neighbors)
        center.add_node(node);
    // identify additional neighbors of s and t
    for (NodeID node : left.nodes)
        for (Neighbor n : node_data[node].neighbors)
            if (center.contains(n.node))
                s_neighbors.push_back(n.node);
    for (NodeID node : right.nodes)
        for (Neighbor n : node_data[node].neighbors)
            if (center.contains(n.node))
                t_neighbors.push_back(n.node);
    util::make_set(s_neighbors);
    util::make_set(t_neighbors);
    // add edges incident to s and t
    for (NodeID node : s_neighbors)
        center.add_edge(s, node, 1, true);
    for (NodeID node : t_neighbors)
        center.add_edge(t, node, 1, true);
    // find minimum cut
    p.cut = min_vertex_cut();
    // revert s-t addition
    // create partition
}

//--------------------------- CutIndex ------------------------------

const uint8_t NON_CUT_VERTEX = 64; // for use with cut_level

// get distance when one vertex is a cut vertex for a subgraph containing both
distance_t direct_distance(const CutIndex &a, const CutIndex &b)
{
    uint16_t a_index = a.distances.size();
    uint16_t b_index = b.distances.size();
    return a_index < b_index ? b.distances[a_index]
        : a_index > b_index ? a.distances[b_index]
        : 0;
}

distance_t get_distance(const CutIndex &a, const CutIndex &b)
{
    // same leaf node, or one vertex is cut vertex
    if (a.partition == b.partition)
        return direct_distance(a, b);
    // find lowest level at which partitions differ
    uint64_t pdiff = a.partition ^ b.partition;
    // partition level used for comparison (upper bound initially)
    int pindex = min(a.cut_level, b.cut_level);
    int diff_level = __builtin_ctz(pdiff); // count trailing zeros
    // one vertex is cut vertex
    if (pindex <= diff_level)
        return direct_distance(a,b);
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
