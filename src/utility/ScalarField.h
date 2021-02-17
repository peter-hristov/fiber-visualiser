#pragma once

#include <cassert>
#include <memory>
#include <qpoint.h>
#include <qvector.h>
#include <string>
#include <vector>
#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif

#include "./MergeTree.h"

namespace tv9k {
namespace utility {
class ScalarField
{
  public:
    ScalarField() {}
    std::string name;
    std::vector<std::vector<std::vector<std::vector<GLfloat>>>> values;

    std::vector<MergeTree> joinTrees;
    std::vector<MergeTree> splitTrees;

    GLfloat min, max;
    std::vector<GLfloat> minValues, maxValues;

    GLfloat currentIsovalue;

    void computeMinMaxValues();
    void computeJoinTrees();
    void computeSplitTrees();

    bool hasJoinTrees() { return this->joinTrees.size() == this->values.size(); }

    bool hasSplitTrees() { return this->splitTrees.size() == this->values.size(); }

    std::string longName = "";
    std::string units = "";

    void clipDomain(int xMin, int xMax, int yMin, int yMax, int zMin, int zMax, int tMin, int tMax);

    // Split domain into cubes of size _downsampleX, int _downsampleY, int
    // _downsampleZ and average their values into a point in the new smaller grid.
    // This really worked best if the domain's dimensions are divisible by these
    // values.
    void downsampleDomain(int _downsampleX, int _downsampleY, int _downsampleZ);

    // For convenience
    int tDim() { return this->values.size(); }
    int xDim()
    {
        assert(this->values.size() > 0);
        return this->values[0].size();
    }
    int yDim()
    {
        assert(this->values[0].size() > 0);
        return this->values[0][0].size();
    }
    int zDim()
    {
        assert(this->values[0][0].size() > 0);
        return this->values[0][0][0].size();
    }
};
}
}
