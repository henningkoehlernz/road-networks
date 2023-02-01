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

// disable expensive query timing
#define REMOVE_REDUNDANT
#define CONTRACT

const size_t repeats = 1;
const size_t nr_queries = 1000000;
const size_t nr_query_tests = 10;
const size_t nr_buckets = 10;
const size_t bucket_size = 10000;
const distance_t bucket_min = 1000;

const size_t MB = 1024 * 1024;

struct ResultData
{
    size_t label_count;
    size_t index_size;
    size_t index_height;
    double index_time;
    double avg_cut_size;
    size_t max_cut_size;
    size_t pruning_2hop;
    size_t pruning_3hop;
    size_t pruning_tail;
    double random_query_time;
    double random_hoplinks;
    vector<double> bucket_query_times;
    vector<double> bucket_hoplinks;
};

#ifdef PRUNING
size_t get_2hop_pruning(const vector<CutIndex> &ci)
{
    size_t total = 0;
    for (NodeID node = 1; node < ci.size(); node++)
        total += ci[node].pruning_2hop;
    return total;
}

size_t get_3hop_pruning(const vector<CutIndex> &ci)
{
    size_t total = 0;
    for (NodeID node = 1; node < ci.size(); node++)
        total += ci[node].pruning_3hop;
    return total;
}

size_t get_tail_pruning(const vector<CutIndex> &ci)
{
    size_t total = 0;
    for (NodeID node = 1; node < ci.size(); node++)
        total += ci[node].pruning_tail;
    return total;
}
#endif

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
        const char* filename = argv[f];
        vector<ResultData> results;
        for (size_t i = 0; i < repeats; i++)
        {
            cout << endl << "reading graph from " << filename << endl;
            fstream fs(filename);
            Graph g;
            read_graph(g, fs);
            fs.close();
            cout << "read " << g.node_count() << " vertices and " << g.edge_count() << " edges" << flush;
            cout << " (diameter=" << g.diameter(false) << ")" << endl;
            // check for redundant edges
            vector<Edge> redundant_edges;
            util::start_timer();
            g.get_redundant_edges(redundant_edges);
#ifdef REMOVE_REDUNDANT
            for (Edge e : redundant_edges)
                g.remove_edge(e.a, e.b);
            cout << "removed " << redundant_edges.size() << " redundant edges in " << util::stop_timer() << "s" << endl;
#else
            cout << "found " << redundant_edges.size() << " redundant edges in " << util::stop_timer() << "s" << endl;
#endif
#ifdef CONTRACT
            size_t old_size = g.node_count();
            vector<Neighbor> closest;
            g.contract(closest);
            cout << "contracted to " << g.node_count() << " vertices (" << g.node_count() * 100 / max<size_t>(1, old_size) << "%) and " << g.edge_count() << " edges" << endl;
#endif
#ifdef NDEBUG
            g.randomize();
#endif
            ResultData result = {};
            // construct index
            vector<CutIndex> ci;
            util::start_timer();
            size_t shortcuts = g.create_cut_index(ci, balance);
#ifdef PRUNING
            result.pruning_2hop = get_2hop_pruning(ci);
            result.pruning_3hop = get_3hop_pruning(ci);
            result.pruning_tail = get_tail_pruning(ci);
#endif
#ifdef CONTRACT
            ContractionIndex con_index(ci, closest);
#else
            ContractionIndex con_index(ci);
#endif
            result.index_time = util::stop_timer();
            result.index_size = con_index.size() / MB;
            result.label_count = con_index.label_count();
            result.index_height = con_index.height();
            result.avg_cut_size = con_index.avg_cut_size();
            result.max_cut_size = con_index.max_cut_size();
            cout << "created index of size " << result.index_size << " MB in " << result.index_time << "s using " << shortcuts << " shortcuts" << endl;
            cout << "#labels=" << result.label_count << ", avg/max cut size=" << setprecision(3) << result.avg_cut_size << "/" << result.max_cut_size << ", height=" << result.index_height << endl;
#ifdef PRUNING
            size_t unpruned_labels = max<size_t>(1, result.label_count + result.pruning_tail);
            cout << "3-HOP pruning could remove " << result.pruning_3hop << " labels (" << result.pruning_3hop * 100 / unpruned_labels << "%)" << endl;
            cout << "2-HOP pruning could remove " << result.pruning_2hop << " labels (" << result.pruning_2hop * 100 / unpruned_labels << "%)" << endl;
            cout << "tail pruning *has* removed " << result.pruning_tail << " labels (" << result.pruning_tail * 100 / unpruned_labels << "%)" << endl;
#endif
            g.reset(); // needed for distance testing

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
                con_index.get_distance(q.first, q.second);
            result.random_query_time = util::stop_timer();
            result.random_hoplinks = con_index.avg_hoplinks(queries);
            cout << "ran " << queries.size() << " random queries in " << result.random_query_time << "s (hoplinks=" << result.random_hoplinks << ")" << endl;

            // test correctness of distance results
            // Dijkstra is slow => reduce number of queries to check
            util::make_set(queries);
            if (queries.size() > nr_query_tests)
                queries.resize(nr_query_tests);
            util::start_timer();
            for (pair<NodeID,NodeID> q : queries)
                if (!con_index.check_query(q, g))
                    return 0;
            cout << "verified " << queries.size() << " queries in " << util::stop_timer() << "s" << endl;

            // test query speed by distance, as for H2H / P2H
            cout << "generating queries by distance: " << flush;
            vector<vector<pair<NodeID,NodeID>>> query_buckets(nr_buckets);
            util::start_timer();
            g.random_pairs(query_buckets, bucket_min, bucket_size, con_index);
            cout << " in " << util::stop_timer() << "s" << endl;
            for (size_t bucket = 0; bucket < query_buckets.size(); bucket++)
            {
                util::start_timer();
                for (pair<NodeID,NodeID> q : query_buckets[bucket])
                    con_index.get_distance(q.first, q.second);
                result.bucket_query_times.push_back(util::stop_timer());
                result.bucket_hoplinks.push_back(con_index.avg_hoplinks(query_buckets[bucket]));
                cout << "ran " << query_buckets[bucket].size() << " queries (bucket " << bucket << ") in " << result.bucket_query_times.back() << "s (hoplinks=" << result.bucket_hoplinks.back() << ")" << endl;
            }
            results.push_back(result);
        }
        if (repeats > 1)
        {
            cout << endl << "Summary for " << filename << ":" << endl << setprecision(5);
            cout << "Index size (MB): " << util::summarize(results, [](ResultData r) -> double { return r.index_size; }) << endl;
            cout << "Index time (s): " << util::summarize(results, [](ResultData r) -> double { return r.index_time; }) << endl;
            cout << "Index height: " << util::summarize(results, [](ResultData r) -> double { return r.index_height; }) << endl;
            cout << "Avg cut size: " << util::summarize(results, [](ResultData r) -> double { return r.avg_cut_size; }) << endl;
            cout << "Max cut size: " << util::summarize(results, [](ResultData r) -> double { return r.max_cut_size; }) << endl;
            cout << "Query time (s): " << util::summarize(results, [](ResultData r) -> double { return r.random_query_time; }) << endl;
            cout << "Avg Hoplinks: " << util::summarize(results, [](ResultData r) -> double { return r.random_hoplinks; }) << endl;
            for (size_t bucket = 0; bucket < nr_buckets; bucket++)
            {
                cout << "Bucket " << bucket << ": time = " << util::summarize(results, [bucket](ResultData r) -> double { return r.bucket_query_times[bucket]; }) * (nr_queries / bucket_size);
                cout << ", hoplinks = " << util::summarize(results, [bucket](ResultData r) -> double { return r.bucket_hoplinks[bucket]; }) << endl;
            }
        }
    }
    return 0;
}
