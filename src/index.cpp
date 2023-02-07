#include "road_network.h"

#include <iostream>

using namespace std;
using namespace road_network;

int main()
{
    // read graph
    Graph g;
    read_graph(g, std::cin);
    vector<Edge> redundant_edges;
    g.get_redundant_edges(redundant_edges);
    for (Edge e : redundant_edges)
        g.remove_edge(e.a, e.b);
    vector<Neighbor> closest;
    g.contract(closest);
#ifdef NDEBUG
    srand(time(nullptr));
    g.randomize();
#endif
    // construct index
    vector<CutIndex> ci;
    g.create_cut_index(ci, 0.2);
    ContractionIndex con_index(ci, closest);
    // write index
    con_index.write(std::cout);
    return 0;
}
