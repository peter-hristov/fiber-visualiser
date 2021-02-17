#include "./MergeTree.h"
#include "./DisjointSet.hpp"
#include "./Utility.h"

#include "./external/tqdm.h"

#include <algorithm>

using namespace std;

MergeTree::MergeTree()
{
    this->xdim = 0;
    this->ydim = 0;
    this->zdim = 0;
}

void
MergeTree::setValues(vector<vector<vector<float>>> _data)
{
    this->data = _data;

    this->xdim = data.size();
    this->ydim = data[0].size();
    this->zdim = data[0][0].size();
}

void
MergeTree::computeJoinTree()
{
    //
    // Flip Values
    //
    for (int i = 0; i < xdim; i++) {
        for (int j = 0; j < ydim; j++) {
            for (int k = 0; k < zdim; k++) {
                this->data[i][j][k] *= -1;
            }
        }
    }

    this->isovalueMult = -1.0;
    this->computeMergeTree();
}

void
MergeTree::computeSplitTree()
{
    this->computeMergeTree();
}

void
MergeTree::computeMergeTree()
{
    this->vis = vector<vector<vector<int>>>(xdim, vector<vector<int>>(ydim, vector<int>(zdim, -1)));

    // startTimer();
    // printf("\n\n------------------------------------ Linearise Array. \n");
    vector<pair<int, float>> sortedVertices = vector<pair<int, float>>(xdim * ydim * zdim);

    for (int i = 0; i < xdim; i++) {
        for (int j = 0; j < ydim; j++) {
            for (int k = 0; k < zdim; k++) {
                int index = Utility::trippleToIndex(make_tuple(i, j, k), xdim, ydim, zdim);
                sortedVertices[index] = make_pair(index, data[i][j][k]);
            }
        }
    }
    // printf("------------------------------------ Linearising done - %3ld.%06ld
    // seconds.\n", endTimer().first, endTimer().second);

    // printf("\n\n------------------------------------ Sort Array. \n");
    sort(sortedVertices.begin(), sortedVertices.end(), [](pair<int, float> a, pair<int, float> b) {
        return a.second < b.second;
    });
    // printf("------------------------------------ Sorting Done - %3ld.%06ld
    // seconds.\n", endTimer().first, endTimer().second);

    //
    // Data structure to hold connected components of the merge tree
    //
    DisjointSet disjointSet(xdim * ydim * zdim);

    //
    // List of parents
    //
    this->mergeTree = vector<int>(xdim * ydim * zdim, -1);

    //
    // Persistent Homology Pairs or equivalently Branches of the Merge Tree
    //
    this->persistencePairs = vector<tuple<int, int, int>>(0);

    //
    // The vertex with the lowest index in its respective connected component
    //
    vector<int> lowestVertex(xdim * ydim * zdim);

    for (int i = 0; i < lowestVertex.size(); i++) {
        lowestVertex[i] = i;
        disjointSet.nodes[sortedVertices[i].first].first = i;
    }

    //
    // Indices instead of values to avoid numerical problems when comparing
    //
    vector<vector<vector<int>>> sortedIndices =
      vector<vector<vector<int>>>(xdim, vector<vector<int>>(ydim, vector<int>(zdim, -1)));

    for (int t = 0; t < sortedVertices.size(); t++) {
        int i, j, k;
        tie(i, j, k) = Utility::indexToTripple(sortedVertices[t].first, xdim, ydim, zdim);
        sortedIndices[i][j][k] = t;
    }

    tqdm bar;

    // printf("\n\n------------------------------------ Computing Merge Tree.
    // \n");
    for (int t = 0; t < sortedVertices.size(); t++) {
        bar.progress(t, sortedVertices.size());

        int currentPosition = sortedVertices[t].first;

        int i, j, k;
        tie(i, j, k) = Utility::indexToTripple(currentPosition, xdim, ydim, zdim);

        for (tuple<int, int, int> n : getAdjacent(i, j, k, xdim, ydim, zdim)) {
            int x, y, z;
            tie(x, y, z) = n;

            int nPosition = Utility::trippleToIndex(n, xdim, ydim, zdim);
            if (sortedIndices[i][j][k] > sortedIndices[x][y][z] &&
                disjointSet.find(currentPosition) != disjointSet.find(nPosition)) {
                mergeTree[lowestVertex[disjointSet.find(nPosition)]] = currentPosition;

                // Save the roots before merging to use as persistence pairs
                int rootCurrentPosition = disjointSet.find(currentPosition);
                int rootNPosition = disjointSet.find(nPosition);

                int merged = disjointSet.merge(currentPosition, nPosition, true);

                if (-1 == merged) {
                    int persistence = disjointSet.nodes[currentPosition].first - disjointSet.nodes[rootNPosition].first;

                    if (0 != persistence) {
                        persistencePairs.push_back(make_tuple(rootNPosition, currentPosition, persistence));
                    }
                } else if (1 == merged) {
                    int persistence =
                      disjointSet.nodes[currentPosition].first - disjointSet.nodes[rootCurrentPosition].first;

                    if (0 != persistence) {
                        persistencePairs.push_back(make_tuple(rootCurrentPosition, currentPosition, persistence));
                    }
                }
                // Something must have gone horribly wrong if we're in the else branch
                else {
                    assert(false);
                }
            }
        }

        vis[i][j][k] = disjointSet.find(currentPosition);
        lowestVertex[disjointSet.find(currentPosition)] = currentPosition;
    }
    bar.finish();

    // Push the Master Branch
    persistencePairs.push_back(
      make_tuple(sortedVertices[0].first, sortedVertices[xdim * ydim * zdim - 1].first, xdim * ydim * zdim - 1));

    //
    // Use colors that correspond to persistence
    //
    std::sort(persistencePairs.begin(), persistencePairs.end(), [](tuple<int, int, int> a, tuple<int, int, int> b) {
        return get<2>(a) > get<2>(b);
    });

    //
    // Compute the 3D visited array that maps vertices to the index of their
    // respective branch
    //
    for (int i = 0; i < persistencePairs.size(); i++) {
        int a, b, c;
        tie(a, b, c) = persistencePairs[i];

        int current = a;
        while (b != current) {
            int x, y, z;
            tie(x, y, z) = Utility::indexToTripple(current, xdim, ydim, zdim);

            this->vis[x][y][z] = i;
            current = mergeTree[current];
        }
    }

    this->simplifiedVisited = vis;
}

//
// Remove all branches with persistence less than this threshold.
// We iterate over all branches, and then for the branches bellow the threshold,
// we set their 3D vis index to -1
//
void
MergeTree::simplifyMergeTree(float threshold)
{
    this->simplifiedVisited = vis;

    for (int i = 0; i < persistencePairs.size(); i++) {
        int a, b, c;
        tie(a, b, c) = persistencePairs[i];

        if (c < threshold) {
            //
            // Set the index of the vertices in that branch to -1
            //
            int i, j, k;
            int current = a;

            while (b != current) {
                tie(i, j, k) = Utility::indexToTripple(current, xdim, ydim, zdim);
                this->simplifiedVisited[i][j][k] = -1;
                current = mergeTree[current];
            }
        }
    }
}

//
// Compute the neighbourhood of (i, j, k)
//
vector<tuple<int, int, int>>
MergeTree::getAdjacent(int i, int j, int k, int xdim, int ydim, int zdim)
{
    // Top
    //
    // /-----/-----/
    // /  2  /  1  /
    // /-----/-----/
    // /  3  /  4  /
    // /-----/-----/
    //
    // Bottom
    //
    // /-----/-----/
    // /  6  /  5  /
    // /-----/-----/
    // /  7  /  8  /
    // /-----/-----/

    // All possible adjacent vertices
    vector<tuple<int, int, int>> adjacent;

    // Square 1
    adjacent.push_back(make_tuple(i, j, k + 1));
    adjacent.push_back(make_tuple(i, j + 1, k));
    adjacent.push_back(make_tuple(i, j + 1, k + 1));
    adjacent.push_back(make_tuple(i + 1, j, k));
    adjacent.push_back(make_tuple(i + 1, j, k + 1));
    adjacent.push_back(make_tuple(i + 1, j + 1, k));
    adjacent.push_back(make_tuple(i + 1, j + 1, k + 1));

    // Square 2
    adjacent.push_back(make_tuple(i - 1, j, k));

    // Square 3
    adjacent.push_back(make_tuple(i, j - 1, k));
    adjacent.push_back(make_tuple(i - 1, j - 1, k));

    // Square 4, None here

    // Square 5
    adjacent.push_back(make_tuple(i, j, k - 1));

    // Square 6
    adjacent.push_back(make_tuple(i - 1, j, k - 1));

    // Square 7
    adjacent.push_back(make_tuple(i, j - 1, k - 1));
    adjacent.push_back(make_tuple(i - 1, j - 1, k - 1));

    // Square 8, None here

    vector<tuple<int, int, int>> neighbours;

    //
    // Filter ones that are outside the grid
    //
    for (auto neighbour : adjacent) {
        int x, y, z;
        tie(x, y, z) = neighbour;

        if (0 <= x && x < xdim && 0 <= y && y < ydim && 0 <= z && z < zdim) {
            neighbours.push_back(neighbour);
        }
    }

    return neighbours;
}

//
// The visited array needs to be adjusted for the current isovaule that is being
// displayed. This is because large components need to absord their child
// branches in order to get a proper color for the isosurface. The way this is
// done is to set the index of the vertex in a branch based on the root of the
// connected component of a sub/super level set (of the merge tree) at the
// isovalue.
//
// 1. Find all leaves
// 2. Go up until you hit the root of the sub/super level set of the tree
// 3. Remember the index of the root and set it to all vertices up that chain
// you went up on
//
vector<vector<vector<int>>>
MergeTree::computeVisitedForIsovalue(float isovalue)
{
    // Correct isovalue for join tree
    isovalue *= isovalueMult;

    vector<vector<vector<int>>> visited = simplifiedVisited;

    //
    // Find the leaves as the nodes which not the parent of any other
    //
    vector<bool> isParent(xdim * ydim * zdim);
    for (int i = 0; i < mergeTree.size(); i++) {
        // Exclude the Root as a leaf
        if (mergeTree[i] != -1) {
            isParent[mergeTree[i]] = true;
        }
    }

    //
    // Save leaves in array
    //
    vector<int> leaves;
    for (int i = 0; i < isParent.size(); i++) {
        if (false == isParent[i]) {
            leaves.push_back(i);
        }
    }

    //
    // Vector to cache whether we have already visited a node so we can stop tree
    // climbing early Reduces complexity to linear.
    //
    vector<bool> endpointFound(xdim * ydim * zdim, false);

    //
    // Go up the tree from every leaf to find and replase their visited value
    //
    for (int i = 0; i < leaves.size(); i++) {
        int index = leaves[i];

        int x, y, z;
        tie(x, y, z) = Utility::indexToTripple(index, xdim, ydim, zdim);

        int endpointValue = visited[x][y][z];

        //
        // Go up the tree to find the root
        //
        while (mergeTree[index] != -1 && data[x][y][z] <= isovalue) {
            index = mergeTree[index];

            if (endpointFound[index]) {
                break;
            }

            tie(x, y, z) = Utility::indexToTripple(index, xdim, ydim, zdim);
            // cout << index << " (" << visited[x][y][z] << ", " <<
            // data->vals[x][y][z] << ") ,";

            endpointFound[index] = true;
        }

        //
        // Save the internal (to this precedure) index of the root
        //
        int lastIndex = index;

        //
        // Save the visited ID of the root
        //
        tie(x, y, z) = Utility::indexToTripple(lastIndex, xdim, ydim, zdim);
        endpointValue = visited[x][y][z];

        //
        // Reset current index
        //
        index = leaves[i];

        //
        // Go up the parents of the leaf to set their index
        //
        tie(x, y, z) = Utility::indexToTripple(index, xdim, ydim, zdim);
        visited[x][y][z] = endpointValue;

        while (index != lastIndex) {
            index = mergeTree[index];
            tie(x, y, z) = Utility::indexToTripple(index, xdim, ydim, zdim);
            visited[x][y][z] = endpointValue;
        }
    }

    return visited;
}

vector<int>
MergeTree::getActiveComponents(float isovalue)
{
    vector<int> components;

    for (int i = 0; i < this->persistencePairs.size(); i++) {
        int x, y, z;
        auto pair = persistencePairs[i];

        tie(x, y, z) = Utility::indexToTripple(get<0>(pair), xdim, ydim, zdim);
        float minVal = data[x][y][z];

        tie(x, y, z) = Utility::indexToTripple(get<1>(pair), xdim, ydim, zdim);
        float maxVal = data[x][y][z];

        if (minVal <= isovalue && isovalue <= maxVal) {
            components.push_back(i);
        }
    }

    return components;
}
