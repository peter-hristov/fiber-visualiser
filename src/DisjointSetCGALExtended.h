#pragma once

#include <CGAL/Union_find.h>

template<typename DataType>
class DisjointSetCGALExtended {
private:
    CGAL::Union_find<int> uf;
public:
        std::map<DataType, int> data;

        DisjointSetCGALExtended() 
        {

        }

        DisjointSetCGALExtended(const std::set<DataType> &preimageGraph) 
        {
            initialize(preimageGraph);
        }

        // Make sure everyone is poiting to the root
        void update()
        {
            //for (int i = 0 ; i < this->data.size() ; i++)
            for (auto &[key, value] : this->data)
            {
                this->uf.find(value);
            }
        }

        std::vector<int> getUniqueRoots()
        {

            std::set<int> uniqueRoots;

            for (auto &[key, value] : this->data)
            {
                uniqueRoots.insert(this->uf.find(value));
            }

            return std::vector<int>(uniqueRoots.begin(), uniqueRoots.end());
        }

        void initialize(const std::set<DataType> &preimageGraph) 
        {
            int n = preimageGraph.size();

            // Initialize the union find
            this->uf(n);

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
            return uf.number_of_sets();
        }

        // Interface for triangles
        int findTriangle(DataType triangle)
        {
            assert(data.contains(triangle));
            //return find(data[triangle]);
            return this->uf.find(data[triangle]);
        }

        void union_setsTriangle(const DataType triangle1, const DataType triangle2) 
        {
            assert(data.contains(triangle1));
            assert(data.contains(triangle2));
            this->uf.unify_sets(data[triangle1], data[triangle2]);
        }

        bool connectedTriangle(const DataType triangle1, const DataType triangle2)
        {
            assert(data.contains(triangle1));
            assert(data.contains(triangle2));
            return this->uf.same_set(data[triangle1], data[triangle2]);
        }

        int find(int x) 
        {
            return this->uf.find(static_cast<int>(x));
        }
};
