#include "road_network.h"
#include "util.h"

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
        cout << "syntax: " << argv[0] << " [balance] <filename> ... <filename>" << endl;
        return 0;
    }
    // check for balance parameter
    double balance = atof(argv[1]);
    size_t file_start = 2;
    if (balance == 0.0)
    {
        balance = 0.25;
        file_start = 1;
    }

    for (int f = file_start; f < argc; f++)
    {
        const char* filename = argv[f];
        cout << "reading graph from " << filename << endl;
        fstream fs(filename);
        Graph g = read_graph(fs);
        cout << "read " << g.node_count() << " vertices and " << g.edge_count() << " edges" << endl;
        DEBUG(g << endl);
        vector<CutIndex> ci;
        auto t_start = chrono::high_resolution_clock::now();
        g.create_cut_index(ci, balance);
        auto t_stop = chrono::high_resolution_clock::now();
        long dur_ms = chrono::duration_cast<chrono::milliseconds>(t_stop - t_start).count();
        cout << "created cut index of size " << index_size(ci) / (1024.0*1024.0)
            << " MB in " << dur_ms / 1000.0 << "s" << endl;
        // test query speed
        g.reset(); // needed for distance testing
        vector<pair<NodeID,NodeID>> queries;
        for (size_t i = 0; i < 1000000; i++)
            queries.push_back(pair(g.random_node(), g.random_node()));
        t_start = chrono::high_resolution_clock::now();
        for (pair<NodeID,NodeID> q : queries)
            get_distance(ci[q.first], ci[q.second]);
        t_stop = chrono::high_resolution_clock::now();
        dur_ms = chrono::duration_cast<chrono::milliseconds>(t_stop - t_start).count();
        cout << "ran " << queries.size() << " queries in " << dur_ms / 1000.0 << "s" << endl;
        // test correctness of distance results
        // Dijkstra is slow => reduce number of queries to check
        util::make_set(queries);
        if (queries.size() > 1000)
            queries.resize(1000);
        t_start = chrono::high_resolution_clock::now();
        for (pair<NodeID,NodeID> q : queries)
            if (!g.check_cut_index(ci, q))
                return 0;
        t_stop = chrono::high_resolution_clock::now();
        dur_ms = chrono::duration_cast<chrono::milliseconds>(t_stop - t_start).count();
        cout << "verified " << queries.size() << " queries in " << dur_ms / 1000.0 << "s" << endl;
    }

    return 0;
}
