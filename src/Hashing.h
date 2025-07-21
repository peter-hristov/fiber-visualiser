#pragma once

#include <set>
#include <array>
#include <utility>
#include <functional>

template <typename T>
struct MyHash
{
    std::size_t operator()(const T &value) const
    {
        return std::hash<T>{}(value);
    }
};

template <>
struct MyHash<std::pair<int, int>>
{
    std::size_t operator()(const std::pair<int, int> &p) const
    {
        std::size_t h1 = std::hash<int>{}(p.first);
        std::size_t h2 = std::hash<int>{}(p.second);
        return h1 ^ (h2 << 1); // Combine the two hashes
    }
};

template <std::size_t N>
struct MyHash<std::array<int, N>>
{
    std::size_t operator()(const std::array<int, N>& arr) const
    {
        std::size_t h = 0;
        for (int x : arr)
        {
            std::size_t hx = std::hash<int>{}(x);
            h ^= hx + 0x9e3779b9 + (h << 6) + (h >> 2);
        }
        return h;
    }
};

template <>
struct MyHash<std::set<int>>
{
    std::size_t operator()(const std::set<int> &s) const
    {
        std::size_t hash_value = 0; // Start with a base hash value

        for (int num : s)
        {
            hash_value ^= std::hash<int>{}(num) + 0x9e3779b9 + (hash_value << 6) + (hash_value >> 2);
        }

        return hash_value;
        return 1;
    }
};

