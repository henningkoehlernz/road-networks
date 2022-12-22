#include "road_network.h"

#include <iostream>
#include <fstream>
#include <string>
#include <chrono>

using namespace std;
using namespace road_network;

#define DEBUG(X) //cerr << X << endl

Graph read_graph(istream &in)
{
    Graph g;
    char line_id;
    uint32_t v, w, d;

    while (in >> line_id) {
        DEBUG("line_id=" << line_id);
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
            DEBUG("adding edge " << v << "->" << w);
            g.add_edge(v, w, d, false);
            break;
        default:
            in.ignore(1000, '\n');
        }
    }

    return g;
}

int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        cout << "syntax: cut <filename> ... <filename>" << endl;
        return 0;
    }

    for (int f = 1; f < argc; f++)
    {
        const char* filename = argv[f];
        cout << "reading graph from " << filename << endl;
        fstream fs(filename);
        Graph g = read_graph(fs);
        cout << "read " << g.node_count() << " vertices and " << g.edge_count() << " edges" << endl;
        DEBUG(g << endl);
        vector<CutIndex> ci;
        auto t_start = chrono::high_resolution_clock::now();
        g.create_cut_index(ci, 0.25);
        auto t_stop = chrono::high_resolution_clock::now();
        long dur_ms = chrono::duration_cast<chrono::milliseconds>(t_stop - t_start).count();
        cout << "created cut index of size " << index_size(ci) / (1024.0*1024.0)
            << " MB in " << dur_ms / 1000.0 << "s" << endl;
    }

    return 0;
}
