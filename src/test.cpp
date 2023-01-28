#include "road_network.h"
#include "util.h"

#include <iostream>
#include <cstdlib>
#include <csignal>

using namespace std;
using namespace road_network;

double balance = 0.25;

class GridEncoder
{
    size_t x_dim;
public:
    GridEncoder(size_t x_dim) : x_dim(x_dim) {}
    NodeID encode(size_t x, size_t y) { return 1 + x + x_dim * y; }
};

// random number in range (inclusive)
int rr(int lower, int upper)
{
    return lower + rand() % (1 + upper - lower);
}

Graph random_grid_graph(size_t x_dim, size_t y_dim)
{
    Graph g(x_dim * y_dim);
    GridEncoder e(x_dim);
    // create full grid
    for (size_t x = 0; x < x_dim; x++)
        for (size_t y = 0; y < y_dim; y++)
        {
            if (x > 0)
                g.add_edge(e.encode(x-1, y), e.encode(x,y), rr(2,3), true);
            if (y > 0)
                g.add_edge(e.encode(x, y-1), e.encode(x,y), rr(2,3), true);
            if (x > 0 && y > 0)
            {
                if (rand() % 2)
                    g.add_edge(e.encode(x-1, y-1), e.encode(x,y), rr(3,4), true);
                else
                    g.add_edge(e.encode(x-1, y), e.encode(x,y-1), rr(3,4), true);
            }
        }
    // sparsify
    vector<Edge> edges;
    g.get_edges(edges);
    for (Edge e : edges)
        if (rand() % 2)
        {
            g.remove_edge(e.a, e.b);
            // ensure graph remains connected
            if (g.get_distance(e.a, e.b, false) == infinity)
                g.add_edge(e.a, e.b, e.d, true);
        }
    return g;
}

void test_crash(Graph &g)
{
    vector<CutIndex> ci;
    g.create_cut_index(ci, balance);
    ContractionIndex con_index(ci);
    con_index.size();
}

void test_index(Graph &g)
{
    vector<CutIndex> ci;
    g.create_cut_index(ci, balance);
    ContractionIndex con_index(ci);
    g.reset();
    vector<pair<NodeID,NodeID>> queries;
    for (size_t i = 0; i < 100; i++)
        queries.push_back(g.random_pair());
    util::make_set(queries);
    for (pair<NodeID,NodeID> q : queries)
        assert(con_index.check_query(q, g));
}

static Graph current_graph;
static vector<Edge> current_edges;
void set_current_graph(const Graph &g)
{
    current_graph = g;
    current_graph.get_edges(current_edges);
}

void abort_handler(int signal)
{
    if (signal != SIGABRT)
    {
        cerr << "called abort_handler with signal=" << signal;
        return;
    }
    cerr << "Aborted on " << current_graph << endl;
    Graph original_graph(current_graph.node_count(), current_edges);
    print_graph(original_graph, cerr);
}

template<typename F>
void run_test(F f, int repeats)
{
    for (size_t x_dim = 2; x_dim < 10; x_dim++)
        for (size_t y_dim = 2; y_dim <= x_dim; y_dim++) {
            for (int i = 0; i < repeats; i++)
            {
                Graph g = random_grid_graph(x_dim, y_dim);
                set_current_graph(g);
                try {
                    f(g);
                } catch (...) {
                    abort();
                }
            }
            cout << "." << flush;
        }
    cout << endl;
}

int main(int argc, char *argv[])
{
#ifdef NDEBUG
    cout << "NDEBUG defined" << endl;
    return 0;
#endif
    int repeats = argc > 1 ? atoi(argv[1]) : 100;
    if (argc > 2)
        balance = atof(argv[2]);
    signal(SIGABRT, abort_handler);
    cout << "Running crash tests " << flush;
    run_test(test_crash, repeats);
    cout << "Running distance tests " << flush;
    run_test(test_index, repeats);
    return 0;
}
