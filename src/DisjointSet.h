#pragma once

#include <cassert>
#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <unordered_map>


// Generic template for hashing (fallback case)
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
        size_t h1 = std::hash<int>{}(p.first);
        size_t h2 = std::hash<int>{}(p.second);
        return h1 ^ (h2 << 1); // Combine the two hashes
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

template<typename DataType>
class DisjointSet {

    public:

        std::vector<int> parent;
        std::vector<int> rank;

        // Map a triangle to an ID (element ID)
        //std::map<std::set<int>, int> data;
        
        std::unordered_map<DataType, int, MyHash<DataType>> data;
        //std::map<DataType, int> data;


        DisjointSet() 
        {

        }

        DisjointSet(const std::set<DataType> &preimageGraph) 
        {
            initialize(preimageGraph);
        }

        bool isEmpty()
        {
            return this->parent.size() == 0;
        }

        void clear()
        {
            this->parent = std::vector<int>();
            this->rank = std::vector<int>();
            this->data  = std::unordered_map<DataType, int, MyHash<DataType>>();
            //this->data  = std::map<DataType, int>();
        }

        // Make sure everyone is poiting to the root
        void update()
        {
            for (int i = 0 ; i < parent.size() ; i++)
            {
                this->find(i);
            }
        }

        std::vector<int> getUniqueRoots()
        {
            std::set<int> uniqueRoots;

            for(int i  = 0 ; i < parent.size() ; i++)
            {
                uniqueRoots.insert(find(i));
            }

            return std::vector<int>(uniqueRoots.begin(), uniqueRoots.end());
        }

        void initialize(const std::set<DataType> &preimageGraph) 
        {
            int n = preimageGraph.size();

            rank.resize(n, 0);

            parent.resize(n);
            for (int i = 0; i < n; ++i) 
            {
                parent[i] = i;
            }

            // Map each triangle to an ID in the disjoint set
            int counter = 0;
            for(const DataType triangle: preimageGraph)
            {
                data[triangle] = counter++;
            }

        }

        // Equals the number of distinct roots
        int countConnectedComponents()
        {
            std::set<int> roots;

            for(int i  = 0 ; i < parent.size() ; i++)
            {
                roots.insert(find(i));
            }

            return roots.size();
        }

        void addElements(const DataType triangle1)
        {
            const int newElementId = parent.size();

            this->parent.push_back(newElementId);
            this->rank.push_back(0);
            this->data[triangle1] = newElementId;
        }

        // Interface for triangles
        int findTriangle(DataType triangle)
        {
            assert(data.contains(triangle));
            return find(data[triangle]);
        }

        void union_setsTriangle(const DataType triangle1, const DataType triangle2) 
        {
            assert(data.contains(triangle1));
            assert(data.contains(triangle2));
            union_sets(data[triangle1], data[triangle2]);

        }

        bool connectedTriangle(const DataType triangle1, const DataType triangle2)
        {
            assert(data.contains(triangle1));
            assert(data.contains(triangle2));
            return connected(data[triangle1], data[triangle2]);
        }

        // Find with path compression
        int find(int x) {
            if (parent[x] != x) {
                // Path compression
                parent[x] = find(parent[x]);  
            }
            return parent[x];
        }

        // Union by rank
        void union_sets(int x, int y) 
        {
            int rootX = find(x);
            int rootY = find(y);

            if (rootX != rootY) 
            {
                if (rank[rootX] > rank[rootY]) 
                {
                    parent[rootY] = rootX;
                } 
                else if (rank[rootX] < rank[rootY]) 
                {
                    parent[rootX] = rootY;
                } 
                else 
                {
                    parent[rootY] = rootX;
                    rank[rootX]++;
                }
            }
        }

        // Check if two nodes are connected
        bool connected(int x, int y) {
            return find(x) == find(y);
        }
};
