#include "road_network.h"
#include "util.h"

#include <iostream>

using namespace std;
using namespace road_network;

Graph sample_graph()
{
    Graph g(10);
    vector<pair<NodeID, NodeID>> edges = {
        pair(1,2), pair(1,3), pair(1,4), pair(2,5), pair(3,6), pair(4,6), pair(5,7), pair(5,8), pair(6,10), pair(7,9), pair(8,9), pair(9,10)
    };
    for (pair<NodeID, NodeID> e : edges)
        g.add_edge(e.first, e.second, 1, true);
    return g;
}

int main (int argc, char *argv[])
{
    cout << "Testing sample graph ..." << endl;
    Graph g = sample_graph();
    cout << g << endl;
    vector<CutIndex> ci;
    g.create_cut_index(ci, 0.25);

    return 0;
}
