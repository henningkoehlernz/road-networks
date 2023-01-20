#include "road_network.h"
#include "util.h"

#include <iostream>
#include <fstream>
#include <time.h>

using namespace std;
using namespace road_network;

int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        cout << "syntax: " << argv[0] << " [balance] <filename>" << endl;
        return 0;
    }
    // check for balance parameter
    double balance = atof(argv[1]);
    const char *filename = argv[1];
    if (balance == 0.0)
        balance = 0.2;
    else
        filename = argv[2];

    srand(time(nullptr));
    cout << endl << "reading graph from " << filename << endl;
    fstream fs(filename);
    Graph g;
    read_graph(g, fs);
    fs.close();
    cout << "read " << g.node_count() << " vertices and " << g.edge_count() << " edges" << flush;
    cout << " (diameter=" << g.diameter(false) << ")" << endl;
    g.randomize();
    util::start_timer();
    Partition p;
    g.create_partition(p, balance);
    cout << "found cut of size " << p.cut.size() << " in " << util::stop_timer() << "s" << endl;

    cerr << "cut: " << p.cut << endl;
    cerr << "left: " << p.left << endl;
    cerr << "right: " << p.right << endl;

    return 0;
}
