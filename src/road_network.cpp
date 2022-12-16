#include "road_network.h"

#include <vector>
#include <queue>
#include <cassert>

using namespace std;
using namespace road_network;

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

//--------------------------- Graph ---------------------------------

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
