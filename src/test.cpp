#include "road_network.h"
#include "util.h"

#include <iostream>
#include <cstdlib>

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
    g.create_cut_index(ci, 0.25);
}

template<typename F>
void run_test(F f)
{
    Graph g = sample_graph();
    f(g);
    for (size_t x_dim = 2; x_dim < 10; x_dim++)
        for (size_t y_dim = 2; y_dim <= x_dim; y_dim++)
            for (int i = 0; i < 100; i++)
            {
                g = random_grid_graph(x_dim, y_dim);
                f(g);
            }
}

int main()
{
    cout << "Running crash tests ..." << endl;
    run_test(test_crash);
    return 0;
}
