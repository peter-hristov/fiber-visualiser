#pragma once

#include <iostream>
#include <vector>
#include <set>
#include <map>

class DisjointSet {

    private:
        std::vector<int> parent;
        std::vector<int> rank;

        // Map a triangle to an ID
        std::map<std::set<int>, int> data;

    public:
        DisjointSet() 
        {

        }

        DisjointSet(const std::set<std::set<int>> &preimageGraph) 
        {
            initialize(preimageGraph);
        }

        void initialize(const std::set<std::set<int>> &preimageGraph) 
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
            for(const std::set<int> triangle: preimageGraph)
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


        // Interface for triangles
        int findTriangle(std::set<int> triangle)
        {
            return find(data[triangle]);
        }

        void union_setsTriangle(const std::set<int> triangle1, const std::set<int> triangle2) 
        {
            union_sets(data[triangle1], data[triangle2]);

        }

        bool connectedTriangle(const std::set<int> triangle1, const std::set<int> triangle2)
        {
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

