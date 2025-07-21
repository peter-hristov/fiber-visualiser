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

template<typename ElementType>
class DisjointSet {

    public:

        std::vector<int> parent;
        std::vector<int> rank;

        // Map a triangle to an ID (element ID)
        //std::map<std::set<int>, int> data;
        
        // Map from the data value to it's index in the data structure
        std::unordered_map<ElementType, int, MyHash<ElementType>> data;
        //std::map<ElementType, int> data;


        DisjointSet() 
        {

        }

        DisjointSet(const std::set<ElementType> &preimageGraph) 
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
            this->data  = std::unordered_map<ElementType, int, MyHash<ElementType>>();
            //this->data  = std::map<ElementType, int>();
        }

        // Make sure everyone is poiting to the root
        void finalise()
        {
            for (int i = 0 ; i < parent.size() ; i++)
            {
                this->findIndex(i);
            }
        }

        std::vector<std::pair<ElementType, int>> getUniqueRepresentativesAndRoots()
        {
            std::set<int> uniqueRoots;
            std::vector<std::pair<ElementType, int>> representativesAndRoots;

            for (const auto &[key, ufId] : this->data)
            {
                const int root = this->findIndex(ufId);

                if (false == uniqueRoots.contains(root))
                {
                    representativesAndRoots.push_back({key, root});
                    uniqueRoots.insert(root);
                }
            }

            return representativesAndRoots;
        }

        std::vector<int> getUniqueRepresentatives()
        {

            std::set<int> uniqueRoots;
            std::vector<ElementType> representatives;

            for (const auto &[key, value] : this->data)
            {
                const int root = this->findTriangle(key);

                if (false == uniqueRoots.contains(root))
                {
                    representatives.push_back(key);
                    uniqueRoots.insert(root);
                }
            }

            //for(int i  = 0 ; i < parent.size() ; i++)
            //{
                //if (uniqueRoots.contains())
                //uniqueRoots.insert(find(i));
            //}

            return representatives;
        }

        std::vector<int> getUniqueRoots()
        {
            std::set<int> uniqueRoots;

            for(int i  = 0 ; i < parent.size() ; i++)
            {
                uniqueRoots.insert(findIndex(i));
            }

            return std::vector<int>(uniqueRoots.begin(), uniqueRoots.end());
        }

        void initialize(const std::set<ElementType> &preimageGraph) 
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
            for(const ElementType triangle: preimageGraph)
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
                roots.insert(findIndex(i));
            }

            return roots.size();
        }

        void addElement(const ElementType &triangle1)
        {
            const int newElementId = parent.size();

            this->parent.push_back(newElementId);
            this->rank.push_back(0);
            this->data[triangle1] = newElementId;
        }

        // Interface for triangles
        int findElement(const ElementType &triangle)
        {
            assert(data.contains(triangle));
            return findIndex(data[triangle]);
        }

        void unionElements(const ElementType &triangle1, const ElementType &triangle2) 
        {
            assert(data.contains(triangle1));
            assert(data.contains(triangle2));
            unionIndices(data[triangle1], data[triangle2]);

        }

        bool connectedElements(const ElementType &triangle1, const ElementType &triangle2)
        {
            assert(data.contains(triangle1));
            assert(data.contains(triangle2));
            return connected(data[triangle1], data[triangle2]);
        }

        // Find with path compression
        int findIndex(const int &x) {
            if (parent[x] != x) {
                // Path compression
                parent[x] = findIndex(parent[x]);  
            }
            return parent[x];
        }

        // Union by rank
        void unionIndices(const int &x, const int &y) 
        {
            const int rootX = findIndex(x);
            const int rootY = findIndex(y);

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
        bool connectedIndex(const int &x, const int &y) {
            return findIndex(x) == findIndex(y);
        }
};
