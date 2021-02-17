#pragma once

#include <vector>

class MergeTree
{
  public:
    MergeTree();
    void setValues(std::vector<std::vector<std::vector<float>>>);
    void computeJoinTree();
    void computeSplitTree();
    void simplifyMergeTree(float);
    std::vector<std::vector<std::vector<int>>> computeVisitedForIsovalue(float);

    //
    // Returns the persistence pairs that include the isovalue
    //
    std::vector<int> getActiveComponents(float);

    void computeMergeTree();

    std::vector<std::tuple<int, int, int>> getAdjacent(int, int, int, int, int, int);

    // Determines whether you're computing the spit tree (1) or join tree (-1).
    float isovalueMult = 1.0;

    // Just for convenience
    int xdim, ydim, zdim;
    // 3D field the merge tree is computed on
    std::vector<std::vector<std::vector<float>>> data;
    // Parent list in he emrge tree
    std::vector<int> mergeTree;
    // Arranged by <lowerIndex, higherIndex, persistence>
    std::vector<std::tuple<int, int, int>> persistencePairs;
    std::vector<std::vector<std::vector<int>>> vis;
    std::vector<std::vector<std::vector<int>>> simplifiedVisited;
};
