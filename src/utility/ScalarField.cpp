#include "./ScalarField.h"

#include "./Utility.h"
#include "src/utility/MergeTree.h"
#include <algorithm>
#include <omp.h>

using namespace std;

void
tv9k::utility::ScalarField::computeMinMaxValues()
{
    this->minValues = std::vector<GLfloat>(this->values.size());
    this->maxValues = std::vector<GLfloat>(this->values.size());

    // Compute the min/max for every time step
    for (int t = 0; t < this->values.size(); t++) {
        minValues[t] = this->values[t][0][0][0];
        maxValues[t] = this->values[t][0][0][0];

        for (int i = 0; i < this->values[0].size(); i++) {
            for (int j = 0; j < this->values[0][0].size(); j++) {
                for (int k = 0; k < this->values[0][0][0].size(); k++) {
                    minValues[t] = std::min(this->minValues[t], this->values[t][i][j][k]);
                    maxValues[t] = std::max(this->maxValues[t], this->values[t][i][j][k]);
                }
            }
        }
    }

    // Compute the min/max for all time steps
    this->min = *min_element(this->minValues.begin(), this->minValues.end());
    this->max = *max_element(this->maxValues.begin(), this->maxValues.end());

    this->currentIsovalue = this->min;
}

void
tv9k::utility::ScalarField::computeJoinTrees()
{
    // Don't do anything if they're already computed
    if (this->joinTrees.size() == this->values.size()) {
        return;
    }

    // Allocate array ot empty trees
    this->joinTrees.resize(this->tDim());

#pragma omp parallel
    {
#pragma omp for
        for (int t = 0; t < this->values.size(); t++) {
            this->joinTrees[t].setValues(this->values[t]);
            this->joinTrees[t].computeJoinTree();
        }
    }
}

void
tv9k::utility::ScalarField::computeSplitTrees()
{
    // Don't do anything if they're already computed
    if (this->splitTrees.size() == this->values.size()) {
        return;
    }

    // Allocate array ot empty trees
    splitTrees.resize(this->tDim());

#pragma omp parallel
    {
#pragma omp for
        for (int t = 0; t < this->values.size(); t++) {
            splitTrees[t].setValues(this->values[t]);
            splitTrees[t].computeSplitTree();
        }
    }
}

// @TODO Deprecated
void
tv9k::utility::ScalarField::clipDomain(int xMin, int xMax, int yMin, int yMax, int zMin, int zMax, int tMin, int tMax)
{
    // Set default values if we have 0s
    // tMax = tMax == 0 ? tDim() - 1 : tMax;
    // xMax = xMax == 0 ? xDim() - 1 : xMax;
    // yMax = yMax == 0 ? yDim() - 1 : yMax;
    // zMax = zMax == 0 ? zDim() - 1 : zMax;

    //// Make sure all dimensions are set correctly
    // assert(tMin >= 0 && tMax < tDim() && tMin <= tMax);
    // assert(xMin >= 0 && xMax < xDim() && xMin <= xMax);
    // assert(yMin >= 0 && yMax < yDim() && yMin <= yMax);
    // assert(zMin >= 0 && zMax < zDim() && zMin <= zMax);

    // this->values.erase(values.begin() + tMax + 1, values.end());
    // this->values.erase(values.begin(), values.begin() + tMin);
    // for (int t = 0; t < tDim(); t++) {
    // this->values[t].erase(values[t].begin() + xMax + 1, values[t].end());
    // this->values[t].erase(values[t].begin(), values[t].begin() + xMin);
    // for (int i = 0; i < xDim(); i++) {
    // this->values[t][i].erase(values[t][i].begin() + yMax + 1, values[t][i].end());
    // this->values[t][i].erase(values[t][i].begin(), values[t][i].begin() + yMin);
    // for (int j = 0; j < yDim(); j++) {
    // this->values[t][i][j].erase(values[t][i][j].begin() + zMax + 1, values[t][i][j].end());
    // this->values[t][i][j].erase(values[t][i][j].begin(), values[t][i][j].begin() + zMin);
    //}
    //}
    //}
}

void
tv9k::utility::ScalarField::downsampleDomain(int downsampleX, int downsampleY, int downsampleZ)
{
    // Nothing to do in the identity downsample
    if (1 == downsampleX && 1 == downsampleY && 1 == downsampleZ) {
        return;
    }

    for (int t = 0; t < tDim(); t += 1) {
        auto newData = std::vector<std::vector<std::vector<GLfloat>>>(
          xDim() / downsampleX + 1,
          std::vector<std::vector<GLfloat>>(yDim() / downsampleY + 1,
                                            std::vector<GLfloat>(zDim() / downsampleZ + 1, 0)));

        for (int i = 0; i < xDim(); i += downsampleX) {
            for (int j = 0; j < yDim(); j += downsampleY) {
                for (int k = 0; k < zDim(); k += downsampleZ) {

                    float sum = 0;
                    int count = 0;

                    for (int x = 0; x < downsampleX && i + x < xDim(); x++) {
                        for (int y = 0; y < downsampleY && j + y < yDim(); y++) {
                            for (int z = 0; z < downsampleZ && k + z < zDim(); z++) {
                                sum += this->values[t][i + x][j + y][k + z];
                                count++;
                            }
                        }
                    }

                    newData[i / downsampleX][j / downsampleY][k / downsampleZ] = sum / count;
                }
            }
        }
        this->values[t] = newData;
    }
}
