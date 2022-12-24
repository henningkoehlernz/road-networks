#include "road_network.h"
#include "util.h"

#include <vector>
#include <queue>
#include <cassert>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <bitset>

using namespace std;

#define DEBUG(X) //cerr << X << endl

namespace road_network {

void log_progress(size_t p, ostream &os = cout)
{
    static const size_t P_DIFF = 1000000L;
    static size_t progress = 0;
    size_t old_log = progress / P_DIFF;
    progress += p;
    size_t new_log = progress / P_DIFF;
    if (old_log < new_log)
    {
        for (size_t i = old_log; i < new_log; i++)
            os << '.';
        os << flush;
    }
}

//--------------------------- CutIndex ------------------------------

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
    int diff_level = __builtin_ctzll(pdiff); // count trailing zeros
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

size_t label_count(const vector<CutIndex> &ci)
{
    size_t total = 0;
    for (const CutIndex &i : ci)
        total += i.distances.size();
    return total;
}

size_t index_size(const vector<CutIndex> &ci)
{
    // vertices start from 1
    return (ci.size() - 1) * (8 + 1 + 2*64 + 4)
        + label_count(ci) * 4;
}

double avg_cut_size(const vector<CutIndex> &ci)
{
    double cut_sum = 0, labels = 0;
    for (size_t i = 1; i < ci.size(); i++)
    {
        cut_sum += ci[i].cut_level + 1;
        // adjust for label pruning
        size_t offset = ci[i].cut_level == 0 ? 0 : ci[i].dist_index[ci[i].cut_level - 1];
        labels += 2 * ci[i].distances.size() - offset + 1;
    }
    return labels / cut_sum;
}

//--------------------------- Graph ---------------------------------

const NodeID NO_NODE = 0;
const SubgraphID NO_SUBGRAPH = 0;

SubgraphID next_subgraph_id(bool reset)
{
    static SubgraphID next_id = 1;
    if (reset)
        next_id = 1;
    return next_id++;
}

Neighbor::Neighbor(NodeID node, distance_t distance) : node(node), distance(distance)
{
}

bool Neighbor::operator<(const Neighbor &other) const
{
    return node < other.node;
}

Node::Node(SubgraphID subgraph_id) : subgraph_id(subgraph_id)
{
    distance = outcopy_distance = 0;
    inflow = outflow = NO_NODE;
    is_redundant = in_partition = in_border = false;
}

Edge::Edge(NodeID a, NodeID b, distance_t d) : a(a), b(b), d(d)
{
}

// definition of static members
vector<Node> Graph::node_data;
NodeID Graph::s, Graph::t;

bool Graph::contains(NodeID node) const
{
    return node_data[node].subgraph_id == subgraph_id;
}

Graph::Graph(size_t node_count)
{
    subgraph_id = next_subgraph_id(true);
    node_data.clear();
    resize(node_count);
    CHECK_CONSISTENT;
}

Graph::Graph(size_t node_count, const vector<Edge> &edges) : Graph(node_count)
{
    for (Edge e : edges)
        add_edge(e.a, e.b, e.d, true);
}

void Graph::resize(size_t node_count)
{
    assert(nodes.empty());
    // node numbering starts from 1, and we reserve two additional nodes for s & t
    node_data.clear();
    node_data.resize(node_count + 3, Node(subgraph_id));
    s = node_count + 1;
    t = node_count + 2;
    node_data[0].subgraph_id = node_data[s].subgraph_id = node_data[t].subgraph_id = NO_SUBGRAPH;
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

void Graph::remove_edge(NodeID v, NodeID w)
{
    std::erase_if(node_data[v].neighbors, [w](const Neighbor &n) { return n.node == w; });
    std::erase_if(node_data[w].neighbors, [v](const Neighbor &n) { return n.node == v; });
}

void Graph::reset()
{
    assign_nodes();
    node_data[s].subgraph_id = NO_SUBGRAPH;
    node_data[t].subgraph_id = NO_SUBGRAPH;
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

size_t Graph::node_count() const
{
    return nodes.size();
}

size_t Graph::edge_count() const
{
    size_t ecount = 0;
    for (NodeID node : nodes)
        for (Neighbor n : node_data[node].neighbors)
            if (contains(n.node))
                ecount++;
    return ecount;
}

void Graph::get_edges(vector<Edge> &edges) const
{
    edges.clear();
    for (NodeID a : nodes)
        for (const Neighbor &n : node_data[a].neighbors)
            if (n.node > a && contains(n.node))
                edges.push_back(Edge(a, n.node, n.distance));
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
    CHECK_CONSISTENT;
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
    CHECK_CONSISTENT;
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
ostream& operator<<(ostream &os, FlowNode fn)
{
    return os << "(" << fn.node << "," << (fn.outcopy ? "T" : "F") << ")";
}

// helper function
bool update_distance(distance_t &d, distance_t d_new)
{
    if (d > d_new)
    {
        d = d_new;
        return true;
    }
    return false;
}

void Graph::run_flow_bfs()
{
    CHECK_CONSISTENT;
    assert(contains(s) && contains(t));
    // init distances
    for (NodeID node : nodes)
        node_data[node].distance = node_data[node].outcopy_distance = infinity;
    node_data[t].distance = node_data[t].outcopy_distance = 0;
    // init queue - start with neighbors of t as t requires special flow handling
    queue<FlowNode> q;
    for (Neighbor n : node_data[t].neighbors)
        if (contains(n.node) && node_data[n.node].outflow != t)
        {
            assert(node_data[n.node].outflow == NO_NODE);
            node_data[n.node].outcopy_distance = 1;
            node_data[n.node].distance = 1; // treat inner-node edges as length 0
            q.push(FlowNode(n.node, true));
        }
    // BFS
    while (!q.empty())
    {
        FlowNode fn = q.front();
        q.pop();

        distance_t fn_dist = fn.outcopy ? node_data[fn.node].outcopy_distance : node_data[fn.node].distance;
        NodeID outflow = node_data[fn.node].outflow;
        // special treatment is needed for node with flow through it
        if (outflow != NO_NODE && fn.outcopy)
        {
            // outflow is only valid neighbor
            if (update_distance(node_data[outflow].distance, fn_dist + 1))
            {
                // need to set distance for 0-distance nodes immediately
                // otherwise a longer path may set wrong distance value first
                update_distance(node_data[outflow].outcopy_distance, fn_dist + 1);
                q.push(FlowNode(outflow, false));
            }
        }
        else
        {
            // when arriving at the incoming copy of flow node, all neighbors except inflow are valid
            // inflow must have been already visited in this case, so checking all neighbors is fine
            for (Neighbor n : node_data[fn.node].neighbors)
            {
                // filter neighbors nodes not belonging to subgraph
                if (!contains(n.node))
                    continue;
                // following outflow by inverting flow requires special handling
                if (n.node == outflow)
                {
                    if (update_distance(node_data[n.node].distance, fn_dist + 1))
                    {
                        // neighbor must be a flow node
                        update_distance(node_data[n.node].outcopy_distance, fn_dist + 1);
                        q.push(FlowNode(n.node, false));
                    }
                }
                else
                {
                    if (update_distance(node_data[n.node].outcopy_distance, fn_dist + 1))
                    {
                        // neighbor may be a flow node
                        if (node_data[n.node].outflow == NO_NODE)
                            update_distance(node_data[n.node].distance, fn_dist + 1);
                        q.push(FlowNode(n.node, true));
                    }
                }
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

void Graph::diff_sort(NodeID v, NodeID w)
{
    CHECK_CONSISTENT;
    // compute distance difference
    size_t node_count = nodes.size();
    vector<pair<int32_t,NodeID>> diff;
    diff.reserve(node_count);
    run_bfs(v);
    for (NodeID node : nodes)
        diff.push_back(pair(node_data[node].distance, node));
    run_bfs(w);
    for (size_t i = 0; i < node_count; i++)
        diff[i].first -= node_data[nodes[i]].distance;
    // sort & replace
    std::sort(diff.begin(), diff.end());
    for (size_t i = 0; i < node_count; i++)
        nodes[i] = diff[i].second;
    DEBUG("diff-sorted: " << *this);
}

vector<NodeID> Graph::min_vertex_cut()
{
    DEBUG("min_vertex_cut over " << *this);
    CHECK_CONSISTENT;
    assert(contains(s) && contains(t));
    // set flow to empty
    for (NodeID node : nodes)
        node_data[node].inflow = node_data[node].outflow = NO_NODE;
#ifndef NDEBUG
    size_t last_s_distance = 1; // min s_distance is 2
#endif
    // find max s-t flow using Dinitz' algorithm
    while (true)
    {
        // construct BFS tree from t
        run_flow_bfs();
        DEBUG("BFS-tree: " << distances());
        const distance_t s_distance = node_data[s].outcopy_distance;
        if (s_distance == infinity)
            break;
        assert(s_distance > last_s_distance && (last_s_distance = s_distance));
        // run DFS from s along inverse BFS tree edges
        vector<NodeID> path;
        vector<FlowNode> stack;
        // iterating over neighbors of s directly simplifies stack cleanup after new s-t path is found
        for (Neighbor sn : node_data[s].neighbors)
        {
            if (!contains(sn.node) || node_data[sn.node].distance != s_distance - 1)
                continue;
            // ensure edge from s to neighbor exists in residual graph
            if (node_data[sn.node].inflow != NO_NODE)
            {
                assert(node_data[sn.node].inflow == s);
                continue;
            }
            stack.push_back(FlowNode(sn.node, false));
            while (!stack.empty())
            {
                FlowNode fn = stack.back();
                stack.pop_back();
                DEBUG("fn=" << fn);
                // clean up path (back tracking)
                distance_t fn_dist = fn.outcopy ? node_data[fn.node].outcopy_distance : node_data[fn.node].distance;
                // safeguard against re-visiting node during DFS (may have been enqueued before first visit)
                if (fn_dist == infinity)
                    continue;
                assert(fn_dist < s_distance && s_distance - fn_dist - 1 <= path.size());
                path.resize(s_distance - fn_dist - 1);
                // increase flow when s-t path is found
                if (fn.node == t)
                {
                    DEBUG("flow path=" << path);
                    assert(node_data[path.front()].inflow == NO_NODE);
                    node_data[path.front()].inflow = s;
                    for (size_t path_pos = 1; path_pos < path.size(); path_pos++)
                    {
                        NodeID from = path[path_pos - 1];
                        NodeID to = path[path_pos];
                        // we might be reverting existing flow
                        // from.inflow may have been changed already => check outflow
                        if (node_data[to].outflow == from)
                        {
                            node_data[to].outflow = NO_NODE;
                            if (node_data[from].inflow == to)
                                node_data[from].inflow = NO_NODE;
                        }
                        else
                        {
                            node_data[from].outflow = to;
                            node_data[to].inflow = from;
                        }
                    }
                    assert(node_data[path.back()].outflow == NO_NODE);
                    node_data[path.back()].outflow = t;
                    // skip to next neighbor of s
                    stack.clear();
                    path.clear();
                    DEBUG("new flow=" << flow());
                    break;
                }
                // ensure vertex is not re-visited during current DFS iteration
                if (fn.outcopy)
                    node_data[fn.node].outcopy_distance = infinity;
                else
                    node_data[fn.node].distance = infinity;
                // continue DFS from node
                path.push_back(fn.node);
                distance_t next_distance = fn_dist - 1;
                // when arriving at outgoing copy of a node with flow through it,
                // we are inverting outflow, so all neighbors are valid (except outflow)
                // otherwise inverting the inflow is the only possible option
                NodeID inflow = node_data[fn.node].inflow;
                if (inflow != NO_NODE && !fn.outcopy)
                {
                    if (node_data[inflow].outcopy_distance == next_distance)
                        stack.push_back(FlowNode(inflow, true));
                }
                else
                {
                    for (Neighbor n : node_data[fn.node].neighbors)
                    {
                        if (!contains(n.node))
                            continue;
                        // inflow inversion requires special handling
                        if (n.node == inflow)
                        {
                            if (node_data[inflow].outcopy_distance == next_distance)
                                stack.push_back(FlowNode(inflow, true));
                        }
                        else
                        {
                            if (node_data[n.node].distance == next_distance)
                                stack.push_back(FlowNode(n.node, false));
                        }
                    }
                }
            }
        }
    }
    // find min cut
    vector<NodeID> cut;
    // node-internal edge appears in cut iff outgoing copy is reachable from t in inverse residual graph and incoming copy is not
    // for node-external edges reachability of endpoint but unreachability of starting point is only possible if endpoint is t
    // in that case, starting point must become the cut vertex
    for (NodeID node : nodes)
    {
        NodeID outflow = node_data[node].outflow;
        // distance already stores distance from t in inverse residual graph
        if (outflow != NO_NODE)
        {
            assert(node_data[node].inflow != NO_NODE);
            if (node_data[node].outcopy_distance < infinity)
            {
                // check inner edge
                if (node_data[node].distance == infinity)
                    cut.push_back(node);
            }
            else
            {
                // check outer edge
                if (outflow == t)
                    cut.push_back(node);
            }
        }
    }
    DEBUG("cut=" << cut);
    return cut;
}

void Graph::get_connected_components(vector<vector<NodeID>> &components)
{
    CHECK_CONSISTENT;
    for (NodeID start_node : nodes)
    {
        // visited nodes are temporarily removed
        if (!contains(start_node))
            continue;
        node_data[start_node].subgraph_id = NO_SUBGRAPH;
        // create new connected component
        components.push_back(vector<NodeID>());
        vector<NodeID> &cc = components.back();
        vector<NodeID> stack;
        stack.push_back(start_node);
        while (!stack.empty())
        {
            NodeID node = stack.back();
            stack.pop_back();
            cc.push_back(node);
            for (Neighbor n : node_data[node].neighbors)
                if (contains(n.node))
                {
                    node_data[n.node].subgraph_id = NO_SUBGRAPH;
                    stack.push_back(n.node);
                }
        }
    }
    // reset subgraph IDs
    assign_nodes();
    DEBUG("components=" << components);
    assert(util::size_sum(components) == nodes.size());
}

void Graph::create_partition(Partition &p, double balance)
{
    CHECK_CONSISTENT;
    assert(nodes.size() > 1);
    // find two extreme points
    NodeID a = nodes[0];
    NodeID b = get_furthest(a);
    a = get_furthest(b);
    // create pre-partition
    diff_sort(a,b);
    // round up if possible
    size_t max_left = min(nodes.size() / 2, static_cast<size_t>(ceil(nodes.size() * balance)));
    size_t min_right = nodes.size() - max_left;
    DEBUG("max_left=" << max_left << ", min_right=" << min_right);
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
    DEBUG("moving " << s_neighbors << " and " << t_neighbors << " to center");
    left.remove_nodes(s_neighbors);
    for (NodeID node : s_neighbors)
        center.add_node(node);
    right.remove_nodes(t_neighbors);
    for (NodeID node : t_neighbors)
        center.add_node(node);
    DEBUG("pre-partition=" << left.nodes << "|" << center.nodes << "|" << right.nodes);
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
    p.cut = center.min_vertex_cut();
    // revert s-t addition
    for (NodeID node : t_neighbors)
    {
        assert(node_data[node].neighbors.back().node == t);
        node_data[node].neighbors.pop_back();
    }
    node_data[t].neighbors.clear();
    for (NodeID node : s_neighbors)
    {
        assert(node_data[node].neighbors.back().node == s);
        node_data[node].neighbors.pop_back();
    }
    node_data[s].neighbors.clear();
    // create partition
    util::make_set(p.cut);
    remove_nodes(p.cut);
    assign_nodes(); // cut vertices stay assigned to center
    vector<vector<NodeID>> components;
    get_connected_components(components);
    sort(components.begin(), components.end(), [](const vector<NodeID> &a, const vector<NodeID> &b) { return a.size() > b.size(); });
    for (const vector<NodeID> &cc : components)
        if (p.left.size() <= p.right.size())
            p.left.insert(p.left.end(), cc.begin(), cc.end());
        else
            p.right.insert(p.right.end(), cc.begin(), cc.end());
    // add cut vertices back to graph
    for (NodeID node : p.cut)
        add_node(node);
    DEBUG("partition=" << p);
    assert(p.left.size() + p.right.size() + p.cut.size() == nodes.size());
    assert(p.left.size() <= nodes.size() - max_left);
    assert(p.right.size() <= nodes.size() - max_left);
}

void Graph::add_shortcuts(const vector<NodeID> &cut, const vector<NodeID> &partition)
{
    CHECK_CONSISTENT;
    DEBUG("add_shortcuts(cut=" << cut << ", partition=" << partition << ")");
    // set flags for partition containment checks
    for (NodeID p : partition)
        node_data[p].in_partition = true;
    // compute border nodes & set flags
    vector<NodeID> border;
    for (NodeID cut_node : cut)
        for (Neighbor n : node_data[cut_node].neighbors)
            if (node_data[n.node].in_partition)
            {
                border.push_back(n.node);
                node_data[n.node].in_border = true;
            }
    util::make_set(border);
    // compute distances, redundancy flags and add shortcuts
    for (NodeID b : border)
    {
        // initialization
        for (NodeID node : nodes)
        {
            node_data[node].distance = infinity;
            node_data[node].is_redundant = true;
        }
        priority_queue<SearchNode> q;
        // redundancy check for b-edges is special, so we handle them here
        for (Neighbor bn : node_data[b].neighbors)
            if (contains(bn.node))
            {
                node_data[bn.node].distance = bn.distance;
                node_data[bn.node].is_redundant = node_data[bn.node].in_partition;
                q.push(SearchNode(bn.distance, bn.node));
            }
        // dijkstra
        while (!q.empty())
        {
            SearchNode next = q.top();
            q.pop();

            bool redundant_or_in_border = node_data[next.node].is_redundant || node_data[next.node].in_border;
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
                    node_data[n.node].is_redundant = redundant_or_in_border;
                    q.push(SearchNode(new_dist, n.node));
                }
                else if (redundant_or_in_border && new_dist == node_data[n.node].distance)
                {
                    node_data[n.node].is_redundant = true;
                }
            }
        }
        // add shortcuts
        for (NodeID x : border)
            if (x != b && !node_data[x].is_redundant)
            {
                DEBUG("shortcut: " << b << "-[" << node_data[x].distance << "]-" << x);
                add_edge(b, x, node_data[x].distance, true);
            }
    }
    // reset flags for border/partition containment checks
    for (NodeID b : border)
        node_data[b].in_border = false;
    for (NodeID p : partition)
        node_data[p].in_partition = false;
}

void Graph::extend_cut_index(std::vector<CutIndex> &ci, double balance, uint8_t cut_level)
{
    //cout << (int)cut_level << flush;
    DEBUG("extend_cut_index at level " << (int)cut_level << " on " << *this);
    CHECK_CONSISTENT;
    assert(cut_level <= 64);
    assert(node_count() > 1);
    const size_t base = ci[nodes[0]].distances.size();
    // find balanced cut
    Partition p;
    if (cut_level < 64)
        create_partition(p, balance);
    else
        p.cut = nodes;
    //cout << "[" << p.cut.size() << "/" << nodes.size() << "]" << flush;
    // compute distances from cut vertices
    for (NodeID c : p.cut)
    {
        run_dijkstra(c);
        for (NodeID node : nodes)
            ci[node].distances.push_back(node_data[node].distance);
        log_progress(nodes.size());
    }
    // truncate distances stored for cut vertices
    for (size_t c_pos = 0; c_pos < p.cut.size(); c_pos++)
        ci[p.cut[c_pos]].distances.resize(base + c_pos);
    // update dist_index
    if (cut_level < 64)
    {
        for (NodeID node : nodes)
            ci[node].dist_index[cut_level] = ci[node].distances.size();
    }
    // set cut_level
    for (NodeID c : p.cut)
        ci[c].cut_level = cut_level;
    // update partition bitstring
    for (NodeID node : p.right)
        ci[node].partition |= (static_cast<uint64_t>(1) << cut_level);

    // add_shortcuts (need to do this before creating subgraphs)
    if (p.left.size() > 1)
        add_shortcuts(p.cut, p.left);
    if (p.right.size() > 1)
        add_shortcuts(p.cut, p.right);

    // recursively extend index for partitions
    if (p.left.size() > 1)
    {
        Graph g(p.left.begin(), p.left.end());
        g.extend_cut_index(ci, balance, cut_level + 1);
    }
    else if (p.left.size() == 1)
        ci[p.left[0]].cut_level = cut_level + 1;

    if (p.right.size() > 1)
    {
        Graph g(p.right.begin(), p.right.end());
        g.extend_cut_index(ci, balance, cut_level + 1);
    }
    else if (p.right.size() == 1)
        ci[p.right[0]].cut_level = cut_level + 1;
}

void Graph::create_cut_index(std::vector<CutIndex> &ci, double balance)
{
    assert(is_undirected());
    ci.clear();
    ci.resize(node_data.size() - 2);
#ifndef NDEBUG
    // sort neighbors to make algorithms deterministic
    for (NodeID node : nodes)
        sort(node_data[node].neighbors.begin(), node_data[node].neighbors.end());
#endif
    extend_cut_index(ci, balance, 0);
}

//--------------------------- Graph debug ---------------------------

bool Graph::is_consistent() const
{
    // all nodes in subgraph have correct subgraph ID assigned
    for (NodeID node : nodes)
        if (node_data[node].subgraph_id != subgraph_id)
        {
            DEBUG("wrong subgraph ID for " << node << " in " << *this);
            return false;
        }
    // number of nodes with subgraph_id of subgraph equals number of nodes in subgraph
    size_t count = 0;
    for (NodeID node = 1; node < node_data.size(); node++)
        if (contains(node))
            count++;
    if (count != nodes.size())
    {
        DEBUG(count << "/" << nodes.size() << " nodes contained in " << *this);
        return false;
    }
    return true;
}

bool Graph::is_undirected() const
{
    for (NodeID node : nodes)
        for (Neighbor n : node_data[node].neighbors)
        {
            bool found = false;
            for (Neighbor nn : node_data[n.node].neighbors)
                if (nn.node == node && nn.distance == n.distance)
                {
                    found = true;
                    break;
                }
            if (!found)
                return false;
        }
    return true;
}

vector<pair<distance_t,distance_t>> Graph::distances() const
{
    vector<pair<distance_t,distance_t>> d;
    for (const Node &n : node_data)
        d.push_back(pair(n.distance, n.outcopy_distance));
    return d;
}

vector<pair<NodeID,NodeID>> Graph::flow() const
{
    vector<pair<NodeID,NodeID>> f;
    for (const Node &n : node_data)
        f.push_back(pair(n.inflow, n.outflow));
    return f;
}

NodeID Graph::random_node() const
{
    return nodes[rand() % nodes.size()];
}

bool Graph::check_cut_index(const vector<CutIndex> &ci, pair<NodeID,NodeID> query)
{
    distance_t d_index = road_network::get_distance(ci[query.first], ci[query.second]);
    distance_t d_dijkstra = get_distance(query.first, query.second, true);
    if (d_index != d_dijkstra)
    {
        cerr << "BUG: d_index=" << d_index << ", d_dijkstra=" << d_dijkstra << endl;
        cerr << "index[" << query.first << "]=" << ci[query.first] << endl;
        cerr << "index[" << query.second << "]=" << ci[query.second] << endl;
    }
    return d_index == d_dijkstra;
}

//--------------------------- ostream -------------------------------

// for easy distance printing
struct Dist
{
    distance_t d;
    Dist(distance_t d) : d(d) {}
};

ostream& operator<<(ostream& os, Dist distance)
{
    if (distance.d == infinity)
        return os << "inf";
    else
        return os << distance.d;
}

ostream& operator<<(ostream& os, const CutIndex &ci)
{
    vector<uint16_t> dist_index(ci.dist_index, ci.dist_index + 64);
    return os << "CI(p=" << bitset<64>(ci.partition) << ",c=" << (int)ci.cut_level
        << ",di=" << dist_index << ",d=" << ci.distances << ")";
}

ostream& operator<<(ostream& os, const Neighbor &n)
{
    if (n.distance == 1)
        return os << n.node;
    else
        return os << n.node << "@" << Dist(n.distance);
}

ostream& operator<<(ostream& os, const Node &n)
{
    return os << "N(" << n.subgraph_id << "#" << n.neighbors << ")";
}

ostream& operator<<(ostream& os, const Partition &p)
{
    return os << "P(" << p.left << "|" << p.cut << "|" << p.right << ")";
}

ostream& operator<<(ostream& os, const Graph &g)
{
    return os << "G(" << g.subgraph_id << "#" << g.nodes << " over " << g.node_data << ")";
}

} // road_network
