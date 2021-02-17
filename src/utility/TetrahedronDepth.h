#pragma once

#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <vector>

namespace TetrahedronDepth {
// The only funcition that is meant to be called from outside
float
getDepth(std::vector<std::vector<float>>, std::vector<std::vector<float>>, bool);

// Helper functions
float dotProduct(std::vector<float>, std::vector<float>);
float getLength(std::vector<float>);
std::vector<float> get2DIntersection(std::vector<float>, std::vector<float>, std::vector<float>, std::vector<float>);
void printVector(std::vector<float>);
std::vector<float> getVectorFromPoints(std::vector<float>, std::vector<float>);
std::vector<float> get2DNormal(std::vector<float>);
std::vector<int> findCrossing(std::vector<std::vector<float>>);
bool isInsideTriangle(std::vector<float>, std::vector<float>, std::vector<float>, std::vector<float>);
std::vector<float>
carthesianToBarycentric(std::vector<std::vector<float>>, std::vector<int>, int);

// Related to testing
std::vector<float>
getRandomPoint(int, float, float);
void
runPermutationsTest();
}
