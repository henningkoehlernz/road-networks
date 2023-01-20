#include "road_network.h"
#include "util.h"

#include <vector>
#include <queue>
#include <cassert>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <bitset>
#include <unordered_set>
#include <thread>
#include <atomic>

using namespace std;

#define DEBUG(X) //cerr << X << endl
#define CUT_DEBUG

// agorithm config
#define DIFF_WEIGHTED
//#define DIFF_SQUARED
//#define CUT_REPEAT 3

#ifdef CUT_BOUNDS
// only store cut bound for every n-th cut vertex
const size_t cut_bound_mod = 20;
#endif

namespace road_network {

// profiling
#ifndef NPROFILE
    static atomic<double> t_partition, t_label, t_shortcut;
    #define START_TIMER util::start_timer()
    #define STOP_TIMER(var) var += util::stop_timer()
#else
    #define START_TIMER
    #define STOP_TIMER(var)
#endif

// progress of 0 resets counter
void log_progress(size_t p, ostream &os = cout)
{
    static const size_t P_DIFF = 1000000L;
    static size_t progress = 0;
    size_t old_log = progress / P_DIFF;
    if (p == 0)
    {
        // terminate progress line & reset
        if (old_log > 0)
            os << endl;
        progress = 0;
        return;
    }
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

// distance addition without overflow
distance_t safe_sum(distance_t a, distance_t b)
{
    return a == infinity || b == infinity ? infinity : a + b;
}

// offset by cut level
uint16_t get_offset(const CutIndex &ci, size_t cut_level)
{
    return cut_level ? ci.dist_index[cut_level - 1] : 0;
}

// compute distance based on given cut level
#ifdef CUT_BOUNDS

distance_t get_cut_level_distance(const CutIndex &a, const CutIndex &b, size_t cut_level)
{
    distance_t min_dist = infinity;
    const uint16_t end_index = a.dist_index[cut_level]; // same for a and b
    const distance_t* a_end = &a.distances[0] + end_index;
    size_t begin_index = end_index > cut_bound_mod ? end_index / cut_bound_mod - 1 : 0;
    do
    {
        // use forward iteration within each inner loop for improved caching
        const distance_t* a_ptr = &a.distances[0] + begin_index * cut_bound_mod;
        const distance_t* b_ptr = &b.distances[0] + begin_index * cut_bound_mod;
        const distance_t* next_a_end = a_ptr;
        while (a_ptr != a_end)
        {
            distance_t dist = safe_sum(*a_ptr, *b_ptr);
            if (dist < min_dist)
                min_dist = dist;
            a_ptr++;
            b_ptr++;
        }
        if (begin_index == 0)
            break;
        begin_index--;
        a_end = next_a_end;
    } while (a.cut_bounds[begin_index] + b.cut_bounds[begin_index] < min_dist);
    return min_dist;
}

#else // CUT_BOUNDS

distance_t get_cut_level_distance(const CutIndex &a, const CutIndex &b, size_t cut_level)
{
    const distance_t* a_end = &a.distances[0] + a.dist_index[cut_level];
#ifdef NO_SHORTCUTS
    const uint16_t offset = 0;
#else
    const uint16_t offset = get_offset(a, cut_level); // same for a and b
#endif
    const distance_t* a_ptr = &a.distances[0] + offset;
    const distance_t* b_ptr = &b.distances[0] + offset;
    // find min 2-hop distance within partition
    distance_t min_dist = infinity;
    while (a_ptr != a_end)
    {
#ifdef NO_SHORTCUTS
        distance_t dist = safe_sum(*a_ptr, *b_ptr);
#else
        distance_t dist = *a_ptr + *b_ptr;
#endif
        if (dist < min_dist)
            min_dist = dist;
        a_ptr++;
        b_ptr++;
    }
    return min_dist;
}

#endif // get_cut_level_distance

// get distance when one vertex is a cut vertex for a subgraph containing both
distance_t direct_distance(const CutIndex &a, const CutIndex &b)
{
    uint16_t a_index = a.distances.size();
    uint16_t b_index = b.distances.size();
    // same node
    if (a_index == b_index)
        return 0;
    // node with more distances values stores distance (within subgraph)
    distance_t dd = a_index < b_index ? b.distances[a_index] : a.distances[b_index];
#ifdef NO_SHORTCUTS
    int cut_level = min(a.cut_level, b.cut_level);
    if (cut_level > 0)
        dd = min(dd, get_cut_level_distance(a, b, cut_level - 1));
#endif
    return dd;
}

distance_t get_distance(const CutIndex &a, const CutIndex &b)
{
    // same leaf node, or one vertex is cut vertex
    if (a.partition == b.partition)
        return direct_distance(a, b);
    // find lowest level at which partitions differ
    int diff_level = __builtin_ctzll(a.partition ^ b.partition); // count trailing zeros
    // a or b might be cut vertex
    if (a.cut_level <= diff_level || b.cut_level <= diff_level)
        return direct_distance(a, b);
    // neither vertex lies in cut
    return get_cut_level_distance(a, b, diff_level);
}

size_t direct_hoplinks(const CutIndex &a, const CutIndex &b)
{
    uint16_t a_index = a.distances.size();
    uint16_t b_index = b.distances.size();
    // same node
    if (a_index == b_index)
        return 0;
#ifdef NO_SHORTCUTS
    int cut_level = min(a.cut_level, b.cut_level);
    return 1 + get_offset(a, cut_level);
#else
    return 1;
#endif
}

size_t get_hoplinks(const CutIndex &a, const CutIndex &b)
{
    // same leaf node, or one vertex is cut vertex
    if (a.partition == b.partition)
        return direct_hoplinks(a, b);
    // find lowest level at which partitions differ
    int diff_level = __builtin_ctzll(a.partition ^ b.partition); // count trailing zeros
    // a or b might be cut vertex
    if (a.cut_level <= diff_level || b.cut_level <= diff_level)
        return direct_hoplinks(a, b);
    // neither vertex lies in cut
#ifdef NO_SHORTCUTS
    return a.dist_index[diff_level];
#else
    return a.dist_index[diff_level] - get_offset(a, diff_level);
#endif
}

double avg_hoplinks(const std::vector<CutIndex> &ci, const vector<pair<NodeID,NodeID>> &queries)
{
    size_t hop_count = 0;
    for (pair<NodeID,NodeID> q : queries)
        hop_count += get_hoplinks(ci[q.first], ci[q.second]);
    return hop_count / (double)queries.size();
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
    size_t total = (ci.size() - 1) * (8 + 1);
    for (NodeID node = 1; node < ci.size(); node++)
    {
        const CutIndex &i = ci[node];
        // no need to account for storing size of dist_index or distances
        // these are already stored inherently
        assert(i.dist_index.size() == i.cut_level + 1u);
        assert(i.distances.size() == i.dist_index[i.cut_level]);
        total += i.distances.size() * 4 + i.dist_index.size() * 2;
#ifdef CUT_BOUNDS
        assert(i.cut_bounds.size() == i.distances.size() / cut_bound_mod);
        total += i.cut_bounds.size() * 4;
#endif
    }
    return total;
}

double avg_cut_size(const vector<CutIndex> &ci)
{
    double cut_sum = 0, labels = 0;
    for (size_t i = 1; i < ci.size(); i++)
    {
        cut_sum += ci[i].cut_level + 1;
        // adjust for label pruning
        size_t offset = get_offset(ci[i], ci[i].cut_level);
        labels += 2 * ci[i].distances.size() - offset + 1;
    }
    return labels / cut_sum;
}

size_t max_cut_size(const vector<CutIndex> &ci)
{
    size_t max_cut = 0;
    for (size_t i = 1; i < ci.size(); i++)
        max_cut = max(max_cut, 1 + ci[i].distances.size() - get_offset(ci[i], ci[i].cut_level));
    return max_cut;
}

size_t index_height(const vector<CutIndex> &ci)
{
    size_t height = 0;
    for (size_t i = 1; i < ci.size(); i++)
        if (ci[i].cut_level > height)
            height = ci[i].cut_level;
    return height;
}

//--------------------------- Graph ---------------------------------

const NodeID NO_NODE = 0;
const SubgraphID NO_SUBGRAPH = 0;

SubgraphID next_subgraph_id(bool reset)
{
    static atomic<SubgraphID> next_id = 1;
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
}

Node& MultiThreadNodeData::operator[](size_type pos)
{
    if (pos == Graph::s)
        return s_data;
    if (pos == Graph::t)
        return t_data;
    return vector::operator[](pos);
}

const Node& MultiThreadNodeData::operator[](size_type pos) const
{
    if (pos == Graph::s)
        return s_data;
    if (pos == Graph::t)
        return t_data;
    return vector::operator[](pos);
}

double Partition::rating() const
{
    size_t l = left.size(), r = right.size(), c = cut.size();
    return min(l, r) / (c * c + 1.0);
}

Edge::Edge(NodeID a, NodeID b, distance_t d) : a(a), b(b), d(d)
{
}

bool Edge::operator<(Edge other) const
{
    return a < other.a
        || (a == other.a && b < other.b)
        || (a == other.a && b == other.b && d < other.d);
}

// definition of static members
thread_local Node MultiThreadNodeData::s_data(NO_SUBGRAPH), MultiThreadNodeData::t_data(NO_SUBGRAPH);
#ifdef MULTI_THREAD
MultiThreadNodeData Graph::node_data;
size_t Graph::thread_threshold;
#else
vector<Node> Graph::node_data;
#endif
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
#ifdef MULTI_THREAD
    thread_threshold = max(node_count / MULTI_THREAD, static_cast<size_t>(1000));
#endif
}

void Graph::add_edge(NodeID v, NodeID w, distance_t distance, bool add_reverse)
{
    assert(v < node_data.size());
    assert(w < node_data.size());
    // check for existing edge
    bool exists = false;
    for (Neighbor &n : node_data[v].neighbors)
        if (n.node == w)
        {
            exists = true;
            n.distance = min(n.distance, distance);
            break;
        }
    if (!exists)
        node_data[v].neighbors.push_back(Neighbor(w, distance));
    if (add_reverse)
        add_edge(w, v, distance, false);
}

void Graph::remove_edge(NodeID v, NodeID w)
{
    std::erase_if(node_data[v].neighbors, [w](const Neighbor &n) { return n.node == w; });
    std::erase_if(node_data[w].neighbors, [v](const Neighbor &n) { return n.node == v; });
}

void Graph::remove_isolated()
{
    unordered_set<NodeID> isolated;
    for (NodeID node : nodes)
        if (degree(node) == 0)
        {
            isolated.insert(node);
            node_data[node].subgraph_id = NO_SUBGRAPH;
        }
    std::erase_if(nodes, [&isolated](NodeID node) { return isolated.contains(node); });
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
    return ecount / 2;
}

size_t Graph::degree(NodeID v) const
{
    assert(contains(v));
    size_t deg = 0;
    for (Neighbor n : node_data[v].neighbors)
        if (contains(n.node))
            deg++;
    return deg;
}

size_t Graph::super_node_count()
{
    return node_data.size() - 3;
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

#ifdef MULTI_THREAD_DISTANCES
void Graph::run_dijkstra_par(const vector<NodeID> &vertices)
{
    CHECK_CONSISTENT;
    vector<thread> threads;
    auto dijkstra = [this](NodeID v, size_t distance_id) {
        assert(contains(v));
        assert(distance_id < MULTI_THREAD_DISTANCES);
        // init distances
        for (NodeID node : nodes)
            node_data[node].distances[distance_id] = infinity;
        node_data[v].distances[distance_id] = 0;
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
                if (new_dist < node_data[n.node].distances[distance_id])
                {
                    node_data[n.node].distances[distance_id] = new_dist;
                    q.push(SearchNode(new_dist, n.node));
                }
            }
        }
    };
    for (size_t i = 0; i < vertices.size(); i++)
        threads.push_back(thread(dijkstra, vertices[i], i));
    for (size_t i = 0; i < vertices.size(); i++)
        threads[i].join();
}
#endif

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
    assert(contains(v) && contains(w));
    weighted ? run_dijkstra(v) : run_bfs(v);
    return node_data[w].distance;
}

pair<NodeID,distance_t> Graph::get_furthest(NodeID v, DistanceMeasure dm)
{
    NodeID furthest = v;
    distance_t max_dist = 0;

    if (dm == DistanceMeasure::mixed)
    {
        // get unweighted distance
        run_bfs(v);
        vector<distance_t> du;
        for (NodeID node : nodes)
            du.push_back(node_data[node].distance);
        // combine with weighted distance
        run_dijkstra(v);
        uint64_t max_d_sqr = 0;
        for (size_t i = 0; i < nodes.size(); i++)
        {
            uint64_t d_sqr = static_cast<uint64_t>(du[i]) * node_data[nodes[i]].distance;
            if (d_sqr > max_d_sqr)
            {
                furthest = nodes[i];
                max_d_sqr = d_sqr;
            }
        }
        max_dist = sqrt(max_d_sqr);
    }
    else
    {
        dm == DistanceMeasure::weighted ? run_dijkstra(v) : run_bfs(v);
        for (NodeID node : nodes)
            if (node_data[node].distance > max_dist)
            {
                furthest = node;
                max_dist = node_data[node].distance;
            }
    }
    return make_pair(furthest, max_dist);
}

Edge Graph::get_furthest_pair(DistanceMeasure dm)
{
    assert(nodes.size() > 1);
    distance_t max_dist = 0;
    NodeID start = nodes[0];
    pair<NodeID,distance_t> furthest = get_furthest(start, dm);
    while (furthest.second > max_dist)
    {
        max_dist = furthest.second;
        start = furthest.first;
        furthest = get_furthest(start, dm);
    }
    return Edge(start, furthest.first, max_dist);
}

distance_t Graph::diameter(bool weighted)
{
    if (nodes.size() < 2)
        return 0;
    return get_furthest_pair(weighted ? DistanceMeasure::weighted : DistanceMeasure::unweighted).d;
}

int64_t sqr_dist(distance_t d)
{
#ifdef DIFF_SQUARED
    if (d == infinity)
        return INT64_MAX;
    return static_cast<int64_t>(d) * static_cast<int64_t>(d);
#else
    return d;
#endif
}

void Graph::diff_sort(NodeID v, NodeID w, bool precomputed)
{
    CHECK_CONSISTENT;
    // compute distance difference
    size_t node_count = nodes.size();
    vector<pair<int64_t,NodeID>> diff;
    diff.reserve(node_count);
    if (!precomputed)
#ifdef DIFF_WEIGHTED
        run_dijkstra(v);
#else
        run_bfs(v);
#endif
    for (NodeID node : nodes)
        diff.push_back(pair(sqr_dist(node_data[node].distance), node));
    #ifdef DIFF_WEIGHTED
        run_dijkstra(w);
    #else
        run_bfs(w);
    #endif
    for (size_t i = 0; i < node_count; i++)
        diff[i].first -= sqr_dist(node_data[nodes[i]].distance);
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
#ifdef NDEBUG
    NodeID a = random_node();
#else
    NodeID a = nodes[0];
#endif
    NodeID b = get_furthest(a, DistanceMeasure::weighted).first;
    a = get_furthest(b, DistanceMeasure::weighted).first;
    // create pre-partition
#ifdef DIFF_WEIGHTED
    diff_sort(b, a, true);
#else
    diff_sort(b, a, false);
#endif
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
    // debug cases of bad cut size - for planar graphs, balanced cuts (BF=1/3) of size 2 sqrt(n) must exist
#ifdef CUT_DEBUG
    if (p.cut.size() * p.cut.size() > 4 * nodes.size())
    {
        double factor = p.cut.size() / sqrt(nodes.size());
        cout << "(cut=" << p.cut.size() << " on " << node_count() << "/" << edge_count() << ", x" << factor << ")";
        cerr << "c cut=" << p.cut.size() << " on " << node_count() << "/" << edge_count() << ", x" << factor << ", a=" << a << ", b=" << b << endl;
        print_graph(*this, cerr);
    }
#endif
}

// half-matrix index for storing half-matrix in flat vector
size_t hmi(size_t a, size_t b)
{
    assert(a != b);
    return a < b ? (b * (b - 1) >> 1) + a : (a * (a - 1) >> 1) + b;
}

void Graph::add_shortcuts(const vector<NodeID> &cut, const vector<CutIndex> &ci)
{
    CHECK_CONSISTENT;
    // compute border nodes
    vector<NodeID> border;
    for (NodeID cut_node : cut)
        for (Neighbor n : node_data[cut_node].neighbors)
            if (contains(n.node))
                border.push_back(n.node);
    util::make_set(border);
    assert(!border.empty());
    // for distance in parent graph we use distances to cut nodes, which must already be in index
    size_t cut_level = ci[cut[0]].cut_level;
    // compute distances between border nodes within subgraph and parent graph
    vector<distance_t> d_partition, d_graph;
#ifdef MULTI_THREAD_DISTANCES
    if (nodes.size() > thread_threshold)
    {
        size_t next_offset;
        for (size_t offset = 0; offset < border.size(); offset = next_offset)
        {
            next_offset = min(offset + MULTI_THREAD_DISTANCES, border.size());
            const vector<NodeID> partial_cut(border.begin() + offset, border.begin() + next_offset);
            run_dijkstra_par(partial_cut);
            for (size_t distance_id = 0; distance_id < partial_cut.size(); distance_id++)
            {
                NodeID n_i = border[distance_id + offset];
                for (size_t j = 0; j < distance_id + offset; j++)
                {
                    NodeID n_j = border[j];
                    distance_t d_ij = node_data[n_j].distances[distance_id];
                    d_partition.push_back(d_ij);
                    distance_t d_cut = get_cut_level_distance(ci[n_i], ci[n_j], cut_level);
                    d_graph.push_back(min(d_ij, d_cut));
                }
            }
        }
    }
    else
#endif
    for (size_t i = 1; i < border.size(); i++)
    {
        NodeID n_i = border[i];
        run_dijkstra(n_i);
        for (size_t j = 0; j < i; j++)
        {
            assert(d_partition.size() == hmi(i, j));
            NodeID n_j = border[j];
            distance_t d_ij = node_data[n_j].distance;
            d_partition.push_back(d_ij);
            distance_t d_cut = get_cut_level_distance(ci[n_i], ci[n_j], cut_level);
            d_graph.push_back(min(d_ij, d_cut));
        }
    }
    // find & add non-redundant shortcuts
    // separate loop as d_graph must be fully computed for redundancy check
    size_t idx_ij = 0;
    for (size_t i = 1; i < border.size(); i++)
    {
        for (size_t j = 0; j < i; j++)
        {
            assert(idx_ij == hmi(i, j));
            distance_t dg_ij = d_graph[idx_ij];
            if (d_partition[idx_ij] > dg_ij)
            {
                bool redundant = false;
                // check for redundancy due to shortest path through third border node k
                for (size_t k = 0; k < border.size(); k++)
                {
                    if (k == i || k == j)
                        continue;
                    if (d_graph[hmi(i, k)] + d_graph[hmi(k, j)] == dg_ij)
                    {
                        redundant = true;
                        break;
                    }
                }
                if (!redundant)
                {
                    DEBUG("shortcut: " << border[i] << "-[" << dg_ij << "]-" << border[j]);
                    add_edge(border[i], border[j], dg_ij, true);
                }
            }
            idx_ij++;
        }
    }
}

void Graph::extend_on_partition(vector<CutIndex> &ci, double balance, uint8_t cut_level, const vector<NodeID> &p, [[maybe_unused]] const vector<NodeID> &cut)
{
    if (p.size() > 1)
    {
        Graph g(p.begin(), p.end());
#ifndef NO_SHORTCUTS
        START_TIMER;
        g.add_shortcuts(cut, ci);
        STOP_TIMER(t_shortcut);
#endif
        g.extend_cut_index(ci, balance, cut_level + 1);
    }
    else if (p.size() == 1)
    {
        ci[p[0]].cut_level = cut_level + 1;
        ci[p[0]].dist_index.push_back(ci[p[0]].distances.size());
    }
}

void Graph::extend_cut_index(vector<CutIndex> &ci, double balance, uint8_t cut_level)
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
    {
        START_TIMER;
        create_partition(p, balance);
#ifdef CUT_REPEAT
        for (size_t i = 1; i < CUT_REPEAT; i++)
        {
            Partition p_new;
            create_partition(p_new, balance);
            if (p_new.rating() > p.rating())
                p = p_new;
        }
#endif
        STOP_TIMER(t_partition);
    }
    else
        p.cut = nodes;
    //cout << "[" << p.cut.size() << "/" << nodes.size() << "]" << flush;
    // compute distances from cut vertices
    START_TIMER;
#ifdef MULTI_THREAD_DISTANCES
    if (nodes.size() > thread_threshold)
    {
        size_t next_offset;
        for (size_t offset = 0; offset < p.cut.size(); offset = next_offset)
        {
            next_offset = min(offset + MULTI_THREAD_DISTANCES, p.cut.size());
            const vector<NodeID> partial_cut(p.cut.begin() + offset, p.cut.begin() + next_offset);
            run_dijkstra_par(partial_cut);
            for (size_t distance_id = 0; distance_id < partial_cut.size(); distance_id++)
            {
                for (NodeID node : nodes)
                    ci[node].distances.push_back(node_data[node].distances[distance_id]);
                log_progress(nodes.size());
            }
        }
    }
    else
#endif
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
    for (NodeID node : nodes)
    {
        assert(ci[node].distances.size() <= UINT16_MAX);
        assert(ci[node].dist_index.size() == cut_level);
        ci[node].dist_index.push_back(ci[node].distances.size());
    }
    // set cut_level
    for (NodeID c : p.cut)
        ci[c].cut_level = cut_level;
    // update partition bitstring
    for (NodeID node : p.right)
        ci[node].partition |= (static_cast<uint64_t>(1) << cut_level);
    STOP_TIMER(t_label);

    // add shortcuts and recurse
#ifdef MULTI_THREAD
    if (nodes.size() > thread_threshold)
    {
        std::thread t_left(extend_on_partition, std::ref(ci), balance, cut_level, std::cref(p.left), std::cref(p.cut));
        extend_on_partition(ci, balance, cut_level, p.right, p.cut);
        t_left.join();
    }
    else
#endif
    {
        extend_on_partition(ci, balance, cut_level, p.left, p.cut);
        extend_on_partition(ci, balance, cut_level, p.right, p.cut);
    }
}

void Graph::create_cut_index(std::vector<CutIndex> &ci, double balance)
{
#ifndef NPROFILE
    t_partition = t_label = t_shortcut = 0;
#endif
    assert(is_undirected());
#ifndef NDEBUG
    // sort neighbors to make algorithms deterministic
    for (NodeID node : nodes)
        sort(node_data[node].neighbors.begin(), node_data[node].neighbors.end());
#endif
    // store original neighbor counts
    vector<NodeID> original_nodes = nodes;
    vector<size_t> original_neighbors(node_data.size());
    for (NodeID node :nodes)
        original_neighbors[node] = node_data[node].neighbors.size();
    // create index
    ci.clear();
    ci.resize(node_data.size() - 2);
    extend_cut_index(ci, balance, 0);
    log_progress(0);
    // reset nodes (top-level cut vertices got removed)
    nodes = original_nodes;
    // remove shortcuts
    for (NodeID node : nodes)
        node_data[node].neighbors.resize(original_neighbors[node], Neighbor(0, 0));
#ifdef CUT_BOUNDS
    // compute cut bounds
    for (NodeID node : nodes)
    {
        distance_t bound = infinity;
        for (size_t i = 0; i < ci[node].distances.size(); i++)
        {
            if (ci[node].distances[i] < bound)
                bound = ci[node].distances[i];
            if ((i + 1) % cut_bound_mod == 0)
                ci[node].cut_bounds.push_back(bound);
        }
    }
#endif
#ifndef NPROFILE
    cerr << "partitioning took " << t_partition << "s" << endl;
    cerr << "labeling took " << t_label << "s" << endl;
    cerr << "shortcuts took " << t_shortcut << "s" << endl;
#endif
}

// returns edges that don't affect distances between nodes
void Graph::get_redundant_edges(vector<Edge> &edges, const vector<CutIndex> &ci) const
{
    CHECK_CONSISTENT;
    assert(edges.empty());
    for (NodeID node : nodes)
        for (Neighbor v : node_data[node].neighbors)
        {
            // only check edges once, and only subgraph edges
            if (v.node < node || !contains(v.node))
                continue;
            for (Neighbor w : node_data[node].neighbors)
                if (w.node != v.node && w.distance + road_network::get_distance(ci[w.node], ci[v.node]) <= v.distance)
                {
                    edges.push_back(Edge(node, v.node, v.distance));
                    break;
                }
        }
}

void Graph::get_redundant_edges(std::vector<Edge> &edges)
{
    CHECK_CONSISTENT;
    assert(edges.empty());
    // reset distances for all nodes
    for (NodeID node : nodes)
        node_data[node].distance = infinity;
    // run localized Dijkstra from each node
    vector<NodeID> visited;
    priority_queue<SearchNode> q;
    for (NodeID v : nodes)
    {
        node_data[v].distance = 0;
        visited.push_back(v);
        distance_t max_dist = 0;
        // init queue - starting from neighbors ensures that only paths of length 2+ are considered
        for (Neighbor n : node_data[v].neighbors)
            if (contains(n.node))
            {
                q.push(SearchNode(n.distance, n.node));
                if (v < n.node)
                    max_dist = max(max_dist, n.distance);
            }
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
                if (new_dist <= max_dist && new_dist < node_data[n.node].distance)
                {
                    node_data[n.node].distance = new_dist;
                    q.push(SearchNode(new_dist, n.node));
                    visited.push_back(n.node);
                }
            }
        }
        // identify redundant edges
        for (Neighbor n : node_data[v].neighbors)
            // only add redundant edges once
            if (v < n.node && contains(n.node) && node_data[n.node].distance <= n.distance)
                edges.push_back(Edge(v, n.node, n.distance));
        // cleanup
        for (NodeID w : visited)
            node_data[w].distance = infinity;
        visited.clear();
    }
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

pair<NodeID,NodeID> Graph::random_pair(size_t steps) const
{
    if (steps < 1)
        return make_pair(random_node(), random_node());
    NodeID start = random_node();
    NodeID stop = start;
    for (size_t i = 0; i < steps; i++)
    {
        NodeID n = NO_NODE;
        do
        {
            n = util::random(node_data[stop].neighbors).node;
        } while (!contains(n));
        stop = n;
    }
    return make_pair(start, stop);
}

// generate batch of random node pairs, filtered into buckets by distance (as for H2H/P2H)
void Graph::random_pairs(vector<vector<pair<NodeID,NodeID>>> &buckets, distance_t min_dist, size_t bucket_size, const vector<CutIndex> &ci)
{
    assert(buckets.size() > 0);
    const distance_t max_dist = diameter(true);
    const double x = pow(static_cast<double>(max_dist) / min_dist, 1.0 / buckets.size());
    vector<distance_t> bucket_caps;
    // don't push last cap - implied and works nicely with std::upper_bound
    for (size_t i = 1; i < buckets.size(); i++)
        bucket_caps.push_back(min_dist * pow(x, i));
    size_t todo = buckets.size();
    cout << "|";
    size_t counter = 0;
    while (todo)
    {
        // generate some queries using random walks for speedup
        pair<NodeID, NodeID> q = ++counter % 5 ? make_pair(random_node(), random_node()) : random_pair(1 + rand() % 100);
        distance_t d = road_network::get_distance(ci[q.first], ci[q.second]);
        if (d >= min_dist)
        {
            size_t bucket = upper_bound(bucket_caps.begin(), bucket_caps.end(), d) - bucket_caps.begin();
            if (buckets[bucket].size() < bucket_size)
            {
                buckets[bucket].push_back(q);
                if (buckets[bucket].size() == bucket_size)
                {
                    todo--;
                    cout << bucket << "|" << flush;
                }
            }
        }
    }
}

void Graph::randomize()
{
    random_shuffle(nodes.begin(), nodes.end());
    for (NodeID node : nodes)
        random_shuffle(node_data[node].neighbors.begin(), node_data[node].neighbors.end());
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

void print_graph(const Graph &g, ostream &os)
{
    vector<Edge> edges;
    g.get_edges(edges);
    sort(edges.begin(), edges.end());
    os << "p sp " << Graph::super_node_count() << " " << edges.size() << endl;
    for (Edge e : edges)
        os << "a " << e.a << ' ' << e.b << ' ' << e.d << endl;
}

void read_graph(Graph &g, istream &in)
{
    char line_id;
    uint32_t v, w, d;

    while (in >> line_id) {
        switch (line_id)
        {
        case 'p':
            in.ignore(3);
            in >> v;
            in.ignore(1000, '\n');
            DEBUG("resizing to " << v);
            g.resize(v);
            break;
        case 'a':
            in >> v >> w >> d;
            g.add_edge(v, w, d, true);
            break;
        default:
            in.ignore(1000, '\n');
        }
    }
    g.remove_isolated();
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
    return os << "CI(p=" << bitset<64>(ci.partition) << ",c=" << (int)ci.cut_level
        << ",di=" << ci.dist_index << ",d=" << ci.distances << ")";
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
