#pragma once

#include <vector>
#include <algorithm>

namespace util {

// sort vector and remove duplicate elements
template<typename T>
void make_set(std::vector<T> &v)
{
    size_t v_size = v.size();
    if (v_size == 0)
        return;
    std::sort(v.begin(), v.end());
    size_t last_distinct = 0;
    for (size_t next = 1; next < v_size; next++)
        if (v[next] != v[last_distinct])
        {
            last_distinct++;
            std::swap(v[next], v[last_distinct]);
        }
    v.resize(last_distinct + 1);
}

// remove elements in set from v
// set must be sorted
template<typename T>
void remove_set(std::vector<T> &v, const std::vector<T> set)
{
    if (v.empty() || set.empty())
        return;
    std::erase_if(v, [&set](T value) { return std::binary_search(set.cbegin(), set.cend(), value); });
}

}
