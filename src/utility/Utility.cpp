#include "Utility.h"

using namespace std;

namespace Utility {
struct timeval startTime;
struct timeval endTime;

QColor
getColorQt(int index)
{
    index = index % 95;
    return cols[index];
}

std::vector<int>
getColor255(int index, float alpha)
{
    index = index % 95;

    return { cols[index].red(), cols[index].green(), cols[index].blue(), cols[index].blue() };
}

int
getColorId(int r, int g, int b)
{
    for (int i = 0; i < cols.size(); i++) {
        if (cols[i].red() == r && cols[i].green() == g && cols[i].blue() == b) {
            return i;
        }
    }

    return -1;
}

std::vector<float>
getColor(int index, int colorType, float alpha)
{
    index = index % 95;

    float red;
    float green;
    float blue;

    // Isosurface Colours
    if (colorType == 0)
    {
        red = (float)(cols[index].red()) / 255.0;
        green = (float)(cols[index].green()) / 255.0;
        blue = (float)(cols[index].blue()) / 255.0;
    }
    // Fibre surface
    else if (colorType == 1)
    {
        red = (float)(colsFS[index].red()) / 255.0;
        green = (float)(colsFS[index].green()) / 255.0;
        blue = (float)(colsFS[index].blue()) / 255.0;
    }
    else if (colorType == 2)
    {
        red = (float)(colsCS[index].red()) / 255.0;
        green = (float)(colsCS[index].green()) / 255.0;
        blue = (float)(colsCS[index].blue()) / 255.0;
    }
    else
    {
        assert(false);
    }

    std::vector<float> color = { red, green, blue, alpha };

    return color;
}

void
startTimer()
{
    gettimeofday(&startTime, NULL);
}

std::pair<long, long>
endTimer()
{
    gettimeofday(&endTime, NULL);

    long nSeconds = endTime.tv_sec - startTime.tv_sec;
    long nMSeconds = endTime.tv_usec - startTime.tv_usec;

    if (nMSeconds < 0) { // negative microsecond delta
        nSeconds--;
        nMSeconds += 1000000;
    } // negative microsecond delta

    return std::make_pair(nSeconds, nMSeconds);
}

bool
sortFunction(std::tuple<int, int, int> a, std::tuple<int, int, int> b)
{
    return std::get<0>(a) > std::get<0>(b);
}

int
trippleToIndex(std::tuple<int, int, int> current, int xdim, int ydim, int zdim)
{
    int i = std::get<0>(current);
    int j = std::get<1>(current);
    int k = std::get<2>(current);

    return j * xdim * zdim + (i * zdim) + k;
}

std::tuple<int, int, int>
indexToTripple(int current, int xdim, int ydim, int zdim)
{
    int j = current / (xdim * zdim);

    int remainder = current % (xdim * zdim);

    int i = remainder / zdim;

    int k = remainder % zdim;

    return std::make_tuple(i, j, k);
}

std::tuple<std::vector<std::vector<std::vector<int>>>, std::vector<std::vector<std::vector<int>>>, int>
computeVolumes(const float isovalue,
               const std::vector<std::vector<std::vector<float>>>& dat,
               int xdim,
               int ydim,
               int zdim)
{
    std::vector<std::vector<std::vector<int>>> vis =
      std::vector<std::vector<std::vector<int>>>(xdim, std::vector<std::vector<int>>(ydim, std::vector<int>(zdim, -1)));
    std::vector<std::vector<std::vector<int>>> vol =
      std::vector<std::vector<std::vector<int>>>(xdim, std::vector<std::vector<int>>(ydim, std::vector<int>(zdim, -1)));

    // startTimer();
    // printf("\n\n------------------------------------ Initialising Disjoit Set.
    // \n");

    DisjointSet disjointSet(xdim * ydim * zdim);

    // printf("------------------------------------ Done - %3ld.%06ld seconds.\n",
    // endTimer().first, endTimer().second);

    // startTimer();
    // printf("------------------------------------ Computing Connected
    // Components. \n");
    for (int j = 0; j < ydim; j++) {
        for (int i = 0; i < xdim; i++) {
            for (int k = 0; k < zdim; k++) {
                if (dat[i][j][k] < isovalue) {
                    auto c = trippleToIndex(std::make_tuple(i, j, k), xdim, ydim, zdim);

                    if (i + 1 < xdim && dat[i + 1][j][k] < isovalue) {
                        disjointSet.merge(c, trippleToIndex(std::make_tuple(i + 1, j, k), xdim, ydim, zdim), false);
                    }
                    if (j + 1 < ydim && dat[i][j + 1][k] < isovalue) {
                        disjointSet.merge(c, trippleToIndex(std::make_tuple(i, j + 1, k), xdim, ydim, zdim), false);
                    }
                    if (k + 1 < zdim && dat[i][j][k + 1] < isovalue) {
                        disjointSet.merge(c, trippleToIndex(std::make_tuple(i, j, k + 1), xdim, ydim, zdim), false);
                    }
                }
            }
        }
    }
    // printf("------------------------------------ Done - %3ld.%06ld seconds.\n",
    // endTimer().first, endTimer().second);

    std::vector<std::tuple<int, int, int>> cc;

    // startTimer();
    // printf("------------------------------------ Computing Map for the
    // Connected Components\n");

    // Collect the roots of the connected compoentns in the cc vector
    for (int i = 0; i < disjointSet.nodes.size(); i++) {
        // If it's the root of a connected component
        if (disjointSet.nodes[i].second == i) {
            auto current = indexToTripple(disjointSet.nodes[i].second, xdim, ydim, zdim);

            int x = std::get<0>(current);
            int y = std::get<1>(current);
            int z = std::get<2>(current);

            // If it's actually bigger than the isovalue
            if (dat[x][y][z] < isovalue) {
                cc.push_back(std::make_tuple(disjointSet.nodes[i].first, disjointSet.nodes[i].second, cc.size()));
            }
        }
    }
    // printf("------------------------------------ Done - %3ld.%06ld seconds.\n",
    // endTimer().first, endTimer().second);

    // Sort the connected components by size
    std::sort(cc.begin(), cc.end(), sortFunction);

    // Labels the roots with the index of their connected component
    for (int i = 0; i < cc.size(); i++) {
        // std::cout << "Connected component " << i << " with root " <<
        // std::get<1>(cc[i]) << " and " << std::get<0>(cc[i]) << " elements." <<
        // std::endl;
        disjointSet.labels[std::get<1>(cc[i])] = i;
    }

    // startTimer();
    // printf("------------------------------------ Computing Visited Array\n");

    // Transfer the labels to each element in the 3D array
    for (int j = 0; j < ydim; j++) {
        for (int i = 0; i < xdim; i++) {
            for (int k = 0; k < zdim; k++) {
                // Obtain the 1D index of the current element
                auto currentIndex = trippleToIndex(std::make_tuple(i, j, k), xdim, ydim, zdim);

                // Get the root of the connected component of the current vertex
                auto root = disjointSet.find(currentIndex);

                // Get the index of this connected compoennt
                auto value = disjointSet.labels[root];

                // Set that index
                vis[i][j][k] = value;

                // @Feature to get volume
                int volume = disjointSet.nodes[disjointSet.find(currentIndex)].first;

                vol[i][j][k] = volume;

                // std::cout << "Volume is : " << volume << std::endl;
            }
        }
    }
    // printf("------------------------------------ Done - %3ld.%06ld
    // seconds.\n\n\n\n", endTimer().first, endTimer().second);

    return std::make_tuple(vis, vol, cc.size());
}

vector<QColor> cols = {
    QColor(255, 0, 0, 255),     QColor(0, 255, 0, 255),     QColor(0, 0, 255, 255),     QColor(255, 255, 0, 255),
    QColor(255, 0, 255, 255),   QColor(0, 255, 255, 255),   QColor(89, 200, 0, 255),    QColor(80., 80., 255),
    QColor(120.0, 80.0, 255),   QColor(140.0, 100.0, 255),  QColor(120.0, 120.0, 255),  QColor(80.0, 120.0, 255),
    QColor(60.0, 100.0, 255),   QColor(80.0, 80.0, 255),    QColor(184, 51, 103, 255),  QColor(109, 32, 205, 255),
    QColor(82, 214, 129, 255),  QColor(124, 172, 132, 255), QColor(215, 182, 8, 255),   QColor(117, 107, 230, 255),
    QColor(171, 251, 21, 255),  QColor(131, 50, 203, 255),  QColor(4, 191, 39, 255),    QColor(48, 53, 70, 255),
    QColor(83, 198, 31, 255),   QColor(243, 254, 104, 255), QColor(29, 242, 255, 255),  QColor(51, 231, 214, 255),
    QColor(144, 199, 23, 255),  QColor(180, 126, 77, 255),  QColor(176, 123, 172, 255), QColor(240, 211, 226, 255),
    QColor(221, 196, 185, 255), QColor(175, 130, 0, 255),   QColor(85, 47, 200, 255),   QColor(239, 20, 23, 255),
    QColor(150, 94, 219, 255),  QColor(34, 39, 250, 255),   QColor(5, 181, 162, 255),   QColor(170, 240, 209, 255),
    QColor(231, 119, 95, 255),  QColor(31, 239, 50, 255),   QColor(131, 199, 142, 255), QColor(33, 164, 32, 255),
    QColor(23, 236, 113, 255),  QColor(217, 40, 9, 255),    QColor(61, 22, 32, 255),    QColor(223, 250, 25, 255),
    QColor(38, 40, 101, 255),   QColor(51, 96, 224, 255),   QColor(14, 43, 51, 255),    QColor(56, 166, 213, 255),
    QColor(197, 53, 164, 255),  QColor(116, 145, 25, 255),  QColor(19, 245, 194, 255),  QColor(104, 174, 217, 255),
    QColor(179, 117, 164, 255), QColor(199, 170, 237, 255), QColor(10, 113, 117, 255),  QColor(177, 142, 44, 255),
    QColor(158, 42, 249, 255),  QColor(12, 226, 105, 255),  QColor(130, 147, 120, 255), QColor(239, 176, 219, 255),
    QColor(152, 38, 242, 255),  QColor(83, 126, 79, 255),   QColor(107, 119, 70, 255),  QColor(47, 160, 71, 255),
    QColor(132, 113, 51, 255),  QColor(141, 133, 58, 255),  QColor(83, 8, 21, 255),     QColor(169, 138, 36, 255),
    QColor(46, 174, 179, 255),  QColor(199, 9, 163, 255),   QColor(15, 30, 28, 255),    QColor(75, 151, 159, 255),
    QColor(207, 136, 199, 255), QColor(175, 188, 116, 255), QColor(165, 5, 136, 255),   QColor(174, 243, 232, 255),
    QColor(134, 233, 236, 255), QColor(34, 210, 119, 255),  QColor(153, 243, 212, 255), QColor(164, 123, 127, 255),
    QColor(56, 130, 205, 255),  QColor(220, 104, 162, 255), QColor(225, 167, 148, 255), QColor(125, 78, 146, 255),
    QColor(141, 165, 17, 255),  QColor(132, 80, 77, 255),   QColor(118, 123, 212, 255), QColor(14, 61, 216, 255),
    QColor(75, 93, 49, 255),    QColor(180, 243, 46, 255),  QColor(221, 11, 24, 255),   QColor(75, 27, 218, 255),
    QColor(180, 107, 129, 255), QColor(134, 220, 229, 255), QColor(112, 38, 102, 255),  QColor(69, 208, 102, 255),
};

vector<QColor> colsFS = {
    QColor(0, 255, 0, 254),     QColor(0, 255, 0, 254),     QColor(255, 0, 0, 254),     QColor(255, 255, 0, 254),
    QColor(255, 0, 255, 254),   QColor(0, 255, 255, 254),   QColor(89, 200, 0, 254),    QColor(80., 80., 254),
    QColor(120.0, 80.0, 254),   QColor(140.0, 100.0, 254),  QColor(120.0, 120.0, 254),  QColor(80.0, 120.0, 254),
    QColor(60.0, 100.0, 254),   QColor(80.0, 80.0, 254),    QColor(184, 51, 103, 254),  QColor(109, 32, 205, 254),
    QColor(82, 214, 129, 254),  QColor(124, 172, 132, 254), QColor(215, 182, 8, 254),   QColor(117, 107, 230, 254),
    QColor(171, 251, 21, 254),  QColor(131, 50, 203, 254),  QColor(4, 191, 39, 254),    QColor(48, 53, 70, 254),
    QColor(83, 198, 31, 254),   QColor(243, 254, 104, 254), QColor(29, 242, 255, 254),  QColor(51, 231, 214, 254),
    QColor(144, 199, 23, 254),  QColor(180, 126, 77, 254),  QColor(176, 123, 172, 254), QColor(240, 211, 226, 254),
    QColor(221, 196, 185, 254), QColor(175, 130, 0, 254),   QColor(85, 47, 200, 254),   QColor(239, 20, 23, 254),
    QColor(150, 94, 219, 254),  QColor(34, 39, 250, 254),   QColor(5, 181, 162, 254),   QColor(170, 240, 209, 254),
    QColor(231, 119, 95, 254),  QColor(31, 239, 50, 254),   QColor(131, 199, 142, 254), QColor(33, 164, 32, 254),
    QColor(23, 236, 113, 254),  QColor(217, 40, 9, 254),    QColor(61, 22, 32, 254),    QColor(223, 250, 25, 254),
    QColor(38, 40, 101, 254),   QColor(51, 96, 224, 254),   QColor(14, 43, 51, 254),    QColor(56, 166, 213, 254),
    QColor(197, 53, 164, 254),  QColor(116, 145, 25, 254),  QColor(19, 245, 194, 254),  QColor(104, 174, 217, 254),
    QColor(179, 117, 164, 254), QColor(199, 170, 237, 254), QColor(10, 113, 117, 254),  QColor(177, 142, 44, 254),
    QColor(158, 42, 249, 254),  QColor(12, 226, 105, 254),  QColor(130, 147, 120, 254), QColor(239, 176, 219, 254),
    QColor(152, 38, 242, 254),  QColor(83, 126, 79, 254),   QColor(107, 119, 70, 254),  QColor(47, 160, 71, 254),
    QColor(132, 113, 51, 254),  QColor(141, 133, 58, 254),  QColor(83, 8, 21, 254),     QColor(169, 138, 36, 254),
    QColor(46, 174, 179, 254),  QColor(199, 9, 163, 254),   QColor(15, 30, 28, 254),    QColor(75, 151, 159, 254),
    QColor(207, 136, 199, 254), QColor(175, 188, 116, 254), QColor(165, 5, 136, 254),   QColor(174, 243, 232, 254),
    QColor(134, 233, 236, 254), QColor(34, 210, 119, 254),  QColor(153, 243, 212, 254), QColor(164, 123, 127, 254),
    QColor(56, 130, 205, 254),  QColor(220, 104, 162, 254), QColor(225, 167, 148, 254), QColor(125, 78, 146, 254),
    QColor(141, 165, 17, 254),  QColor(132, 80, 77, 254),   QColor(118, 123, 212, 254), QColor(14, 61, 216, 254),
    QColor(75, 93, 49, 254),    QColor(180, 243, 46, 254),  QColor(221, 11, 24, 254),   QColor(75, 27, 218, 254),
    QColor(180, 107, 129, 254), QColor(134, 220, 229, 254), QColor(112, 38, 102, 254),  QColor(69, 208, 102, 254),
};

vector<QColor> colsCS = {
    QColor(0, 0, 255, 254),     QColor(0, 255, 0, 254),     QColor(255, 0, 0, 254),     QColor(255, 255, 0, 254),
    QColor(255, 0, 255, 254),   QColor(0, 255, 255, 254),   QColor(89, 200, 0, 254),    QColor(80., 80., 254),
    QColor(120.0, 80.0, 254),   QColor(140.0, 100.0, 254),  QColor(120.0, 120.0, 254),  QColor(80.0, 120.0, 254),
    QColor(60.0, 100.0, 254),   QColor(80.0, 80.0, 254),    QColor(184, 51, 103, 254),  QColor(109, 32, 205, 254),
    QColor(82, 214, 129, 254),  QColor(124, 172, 132, 254), QColor(215, 182, 8, 254),   QColor(117, 107, 230, 254),
    QColor(171, 251, 21, 254),  QColor(131, 50, 203, 254),  QColor(4, 191, 39, 254),    QColor(48, 53, 70, 254),
    QColor(83, 198, 31, 254),   QColor(243, 254, 104, 254), QColor(29, 242, 255, 254),  QColor(51, 231, 214, 254),
    QColor(144, 199, 23, 254),  QColor(180, 126, 77, 254),  QColor(176, 123, 172, 254), QColor(240, 211, 226, 254),
    QColor(221, 196, 185, 254), QColor(175, 130, 0, 254),   QColor(85, 47, 200, 254),   QColor(239, 20, 23, 254),
    QColor(150, 94, 219, 254),  QColor(34, 39, 250, 254),   QColor(5, 181, 162, 254),   QColor(170, 240, 209, 254),
    QColor(231, 119, 95, 254),  QColor(31, 239, 50, 254),   QColor(131, 199, 142, 254), QColor(33, 164, 32, 254),
    QColor(23, 236, 113, 254),  QColor(217, 40, 9, 254),    QColor(61, 22, 32, 254),    QColor(223, 250, 25, 254),
    QColor(38, 40, 101, 254),   QColor(51, 96, 224, 254),   QColor(14, 43, 51, 254),    QColor(56, 166, 213, 254),
    QColor(197, 53, 164, 254),  QColor(116, 145, 25, 254),  QColor(19, 245, 194, 254),  QColor(104, 174, 217, 254),
    QColor(179, 117, 164, 254), QColor(199, 170, 237, 254), QColor(10, 113, 117, 254),  QColor(177, 142, 44, 254),
    QColor(158, 42, 249, 254),  QColor(12, 226, 105, 254),  QColor(130, 147, 120, 254), QColor(239, 176, 219, 254),
    QColor(152, 38, 242, 254),  QColor(83, 126, 79, 254),   QColor(107, 119, 70, 254),  QColor(47, 160, 71, 254),
    QColor(132, 113, 51, 254),  QColor(141, 133, 58, 254),  QColor(83, 8, 21, 254),     QColor(169, 138, 36, 254),
    QColor(46, 174, 179, 254),  QColor(199, 9, 163, 254),   QColor(15, 30, 28, 254),    QColor(75, 151, 159, 254),
    QColor(207, 136, 199, 254), QColor(175, 188, 116, 254), QColor(165, 5, 136, 254),   QColor(174, 243, 232, 254),
    QColor(134, 233, 236, 254), QColor(34, 210, 119, 254),  QColor(153, 243, 212, 254), QColor(164, 123, 127, 254),
    QColor(56, 130, 205, 254),  QColor(220, 104, 162, 254), QColor(225, 167, 148, 254), QColor(125, 78, 146, 254),
    QColor(141, 165, 17, 254),  QColor(132, 80, 77, 254),   QColor(118, 123, 212, 254), QColor(14, 61, 216, 254),
    QColor(75, 93, 49, 254),    QColor(180, 243, 46, 254),  QColor(221, 11, 24, 254),   QColor(75, 27, 218, 254),
    QColor(180, 107, 129, 254), QColor(134, 220, 229, 254), QColor(112, 38, 102, 254),  QColor(69, 208, 102, 254),
};
}
