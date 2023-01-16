#include "road_network.h"
#include "util.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <thread>
#include <time.h>
#include <sys/resource.h>

using namespace std;
using namespace road_network;

#define DEBUG(X) //cerr << X << endl
// disable expensive query timing
//#define NQUERY
#define REMOVE_REDUNDANT
// repeat index generation & queries for each graph
#define REPEATS 10

const size_t nr_queries = 1000000;
const size_t nr_query_tests = 10;

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
            g.add_edge(v, w, d, true);
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
    int file_start = 2;
    if (balance == 0.0)
    {
        balance = 0.2;
        file_start = 1;
    }

#ifdef NO_SHORTCUTS
    cout << "shortcuts disabled" << endl;
#else
    cout << "shortcuts enabled" << endl;
#endif

#ifdef MULTI_THREAD
    cout << "multi-threading enabled" << endl;
    cout << "threads supported by hardware: " << thread::hardware_concurrency() << endl;
#else
    cout << "multi-threading disabled" << endl;
#endif

#ifdef NDEBUG
    srand(time(nullptr));
#endif
    for (int f = file_start; f < argc; f++)
    {
#ifdef REPEATS
        for (size_t i = 0; i < REPEATS; i++)
        {
#endif
        const char* filename = argv[f];
        cout << endl << "reading graph from " << filename << endl;
        fstream fs(filename);
        Graph g = read_graph(fs);
        fs.close();
        cout << "read " << g.node_count() << " vertices and " << g.edge_count() << " edges" << flush;
        cout << " (diameter=" << g.diameter(false) << ")" << endl;
        DEBUG(g << endl);
#ifdef REMOVE_REDUNDANT
        vector<Edge> redundant_edges;
        util::start_timer();
        g.get_redundant_edges(redundant_edges);
        for (Edge e : redundant_edges)
            g.remove_edge(e.a, e.b);
        cout << "removed " << redundant_edges.size() << " redundant edges in " << util::stop_timer() << "s" << endl;
#endif
#ifdef NDEBUG
        g.randomize();
#endif
        // construct index
        vector<CutIndex> ci;
        util::start_timer();
        g.create_cut_index(ci, balance);
        double duration = util::stop_timer();
        cout << "created cut index of size " << index_size(ci) / (1024*1024) << " MB in " << duration << "s" << endl;
        cout << "#labels=" << label_count(ci) << ", avg/max cut size=" << setprecision(3) << avg_cut_size(ci) << "/" << max_cut_size(ci) << ", height=" << index_height(ci) << endl;
        g.reset(); // needed for distance testing
#ifndef REMOVE_REDUNDANT
        g.get_redundant_edges(redundant_edges, ci);
        cout << "redundant edges: " << redundant_edges.size() << endl;
#endif

        // show memory consumption
        rusage usage;
        if (getrusage(RUSAGE_SELF, &usage) != -1)
            cout << "maximum memory used: " << usage.ru_maxrss / 1024 << " MB" << endl;

        // test query speed
        vector<pair<NodeID,NodeID>> queries;
        for (size_t i = 0; i < nr_queries; i++)
            queries.push_back(g.random_pair());
        util::start_timer();
        for (pair<NodeID,NodeID> q : queries)
            get_distance(ci[q.first], ci[q.second]);
        duration = util::stop_timer();

        cout << "ran " << queries.size() << " random queries in " << duration << "s (hoplinks=" << avg_hoplinks(ci, queries) << ")" << endl;
#ifndef NQUERY
        // test correctness of distance results
        // Dijkstra is slow => reduce number of queries to check
        util::make_set(queries);
        if (queries.size() > nr_query_tests)
            queries.resize(nr_query_tests);
        util::start_timer();
        for (pair<NodeID,NodeID> q : queries)
            if (!g.check_cut_index(ci, q))
                return 0;
        duration = util::stop_timer();
        cout << "verified " << queries.size() << " queries in " << duration << "s" << endl;
        // test query speed by distance, as for H2H / P2H
        cout << "generating queries by distance: " << flush;
        vector<vector<pair<NodeID,NodeID>>> query_buckets(10);
        util::start_timer();
        g.random_pairs(query_buckets, 1000, 10000, ci);
        cout << " in " << util::stop_timer() << "s" << endl;
        for (size_t bucket = 0; bucket < query_buckets.size(); bucket++)
        {
            util::start_timer();
            for (pair<NodeID,NodeID> q : query_buckets[bucket])
                get_distance(ci[q.first], ci[q.second]);
            duration = util::stop_timer();
            cout << "ran " << query_buckets[bucket].size() << " queries (bucket " << bucket << ") in " << duration << "s (hoplinks=" << avg_hoplinks(ci, query_buckets[bucket]) << ")" << endl;
        }
#endif
#ifdef REPEATS
        }
#endif
    }
    return 0;
}
