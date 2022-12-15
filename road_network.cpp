#include "road_network.h"

#include <vector>
#include <queue>

using namespace std;
using namespace road_network;

SubgraphID next_subgraph_id()
{
    static SubgraphID next_id = 0;
    return next_id++;
}

vector<Node> Graph::node_data; // definition of static member

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

distance_t Graph::get_distance(NodeID u, NodeID v) const
{
    // init distances
    for (NodeID node : nodes)
        node_data[node].distance = infinity;
    node_data[u].distance = 0;
    // init queue
    priority_queue<SearchNode> queue;
    queue.push(SearchNode(0, u));
    // dijkstra
    while (!queue.empty()) {
        SearchNode next = queue.top();
        queue.pop();

        for (Neighbor n : node_data[next.node].out) {
            // filter out nodes not belonging to subgraph
            if (node_data[n.node].subgraph_id != subgraph_id)
                continue;
            // update distance and enque
            distance_t new_dist = next.distance + n.distance;
            if (new_dist < node_data[n.node].distance) {
                node_data[n.node].distance = new_dist;
                queue.push(SearchNode(new_dist, n.node));
            }
        }
    }
    return node_data[v].distance;
}
