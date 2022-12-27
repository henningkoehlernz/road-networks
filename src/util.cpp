#include "util.h"

#include <chrono>

using namespace std;

namespace util {

// stack of times of start_timer calls
static vector<chrono::time_point<std::chrono::high_resolution_clock>> start_times;

void start_timer()
{
    start_times.push_back(chrono::high_resolution_clock::now());
}

long stop_timer()
{
    auto stop_time = chrono::high_resolution_clock::now();
    long diff_ms = chrono::duration_cast<chrono::milliseconds>(stop_time - start_times.back()).count();
    start_times.pop_back();
    return diff_ms;
}

}
