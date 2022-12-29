#pragma once

#include <vector>
#include <algorithm>
#include <ostream>

namespace util {

// start new time measurement
void start_timer();
// returns time in seconds since last unconsumed start_timer call and consumes it
double stop_timer();

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

// compute total number of elements in vector of collections
template<typename T>
size_t size_sum(const std::vector<T> &v)
{
    size_t sum = 0;
    for (const T &x : v)
        sum += x.size();
    return sum;
}

// extract size values of vector of collection
template<typename T>
std::vector<size_t> sizes(const std::vector<T> &v)
{
    std::vector<size_t> s;
    for (const T &x : v)
        s.push_back(x.size());
    return s;
}

template<typename T>
T random(const std::vector<T> &v)
{
    assert(v.size() > 0);
    return v[rand() % v.size()];
}

} // util

namespace std {

template <typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T> &v)
{
    if (v.empty())
        return os << "[]";
    os << "[0:" << v[0];
    for (size_t i = 1; i < v.size(); i++)
        os << ',' << i << ":" << v[i];
    return os << ']';
}

template <typename A, typename B>
std::ostream& operator<<(std::ostream& os, const std::pair<A,B> &p)
{
    return os << "(" << p.first << "," << p.second << ")";
}

} // std
