#pragma once

#include <QColor>
#include <algorithm>
#include <sys/time.h>
#include <vector>

#include "DisjointSet.hpp"

namespace Utility {
//
// Getting Color for the 3D render and the 2D plot
//
std::vector<float>
getColor(int, int, float);
std::vector<int>
getColor255(int, float);

std::vector<int>
getColor255(int, float);

QColor
getColorQt(int);

int
getColorId(int r, int g, int b);

// extern is used to that he linker can defer the loopup to the source file
extern std::vector<QColor> cols;
// Fiber surface
extern std::vector<QColor> colsFS;
// Combined surface
extern std::vector<QColor> colsCS;

//
// Computing Connected Components for Volume Estimation
//
int
trippleToIndex(std::tuple<int, int, int> current, int xdim, int ydim, int zdim);
std::tuple<int, int, int>
indexToTripple(int current, int xdim, int ydim, int zdim);
bool sortFunction(std::tuple<int, int, int>, std::tuple<int, int, int>);

// Compute connected components and return a 3D array with indices from 0 to inf
// of the index of the connected components, ordered by size and a 3D array with
// the sizes of the connected components made for exporting the last output is
// the total number of connected components
std::tuple<std::vector<std::vector<std::vector<int>>>, std::vector<std::vector<std::vector<int>>>, int>
computeVolumes(const float, const std::vector<std::vector<std::vector<float>>>& dat, int, int, int);

//
// Stopwatch for measuring performance
//
// @TODO Why does this need to be inline?
extern struct timeval startTime;
extern struct timeval endTime;

void
startTimer();
std::pair<long, long>
endTimer();
}
