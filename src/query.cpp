#include "road_network.h"
#include "util.h"

#include <iostream>

using namespace std;
using namespace road_network;

const size_t nr_queries = 1000000;
const size_t MB = 1024 * 1024;

int main()
{
#ifdef NDEBUG
    srand(time(nullptr));
#endif
    // read index
    util::start_timer();
    ContractionIndex con_index(std::cin);
    double read_index_time = util::stop_timer();
    cout << "read index in " << read_index_time << "s (" << con_index.size() / MB << " MB)" << endl;
    // test query speed
    vector<pair<NodeID,NodeID>> queries;
    for (size_t i = 0; i < nr_queries; i++)
        queries.push_back(con_index.random_query());
    util::start_timer();
    for (pair<NodeID,NodeID> q : queries)
        con_index.get_distance(q.first, q.second);
    double random_query_time = util::stop_timer();
    double random_hoplinks = con_index.avg_hoplinks(queries);
    cout << "ran " << queries.size() << " random queries in " << random_query_time << "s (hoplinks=" << random_hoplinks << ")" << endl;
    return 0;
}
