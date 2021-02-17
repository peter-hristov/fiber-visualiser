#pragma once

#include <assert.h>
#include <iomanip>
#include <iostream>
#include <utility>
#include <vector>

// Vertex data structure used by RTF
class DisjointSet
{
    // <value, parent>

  public:
    // @TODO
    // Add an enum for the merge style

    std::vector<std::pair<int, int>> nodes;

    //
    // The labels are used to indicate the order of a connected connected
    // component in terms of the number of elements it has
    //
    std::vector<int> labels;

    DisjointSet(int size)
    {
        nodes = std::vector<std::pair<int, int>>(size, std::pair<int, int>(1, -1));
        labels = std::vector<int>(size, -1);

        for (int i = 0; i < nodes.size(); i++) {
            nodes[i].second = i;
        }
    }

    void erase() { nodes.erase(nodes.begin(), nodes.end()); }

    void add() { this->nodes.push_back(std::make_pair(1, this->nodes.size())); }

    int find(int i)
    {
        assert(i >= 0);
        assert(i < nodes.size());

        if (nodes[i].second != i) {
            nodes[i].second = find(nodes[i].second);
        }

        return nodes[i].second;

        //
        // Iterative Path Compression
        //

        // Find root
        // auto root = nodes[i].second;

        // while(nodes[root].second != root)
        //{
        // root = nodes[root].second;
        //}

        //// Compress Path
        // while (nodes[i].second != i)
        //{
        // int parent = nodes[i].second;
        // nodes[i].second = root;

        // i = parent;
        //}

        // return root;
    }

    int merge(int i, int j, bool mergeByIndex)
    {
        assert(i >= 0);
        assert(i < nodes.size());

        assert(j >= 0);
        assert(j < nodes.size());

        int rootI = find(i);
        int rootJ = find(j);

        if (rootI == rootJ) {
            return 0;
        }

        if (true == mergeByIndex) {
            if (nodes[rootI].first < nodes[rootJ].first) {
                nodes[rootJ].second = rootI;
                return -1;
            } else {
                nodes[rootI].second = rootJ;
                return 1;
            }
        }
        // Merge by size
        else {
            if (nodes[rootI].first > nodes[rootJ].first) {
                nodes[rootJ].second = rootI;
                nodes[rootI].first += nodes[rootJ].first;
                return -1;
            } else
            // else if (nodes[rootI].first < nodes[rootJ].first)
            {
                nodes[rootI].second = rootJ;
                nodes[rootJ].first += nodes[rootI].first;
                return 1;
            }
        }

        // Something went horribly wrong
        assert(false);
    }

    void print()
    {
        std::cout << "\n\nIndices";

        int w = 4;

        for (auto t = 0; t < this->nodes.size(); t++) {
            std::cout << std::setw(w) << t;
        }

        std::cout << std::endl;

        std::cout << "Values ";
        for (auto t : this->nodes) {
            std::cout << std::setw(w) << t.first;
        }

        std::cout << std::endl;

        std::cout << "Parents";
        for (auto t : this->nodes) {
            std::cout << std::setw(w) << t.second;
        }

        std::cout << std::endl;
    }
};
