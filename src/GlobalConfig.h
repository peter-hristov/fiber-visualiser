#pragma once

#include <algorithm>
#include <cassert>
#include <iostream>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

namespace tv9k {
class InputInformation
{
  public:
    std::string filename;
    std::string attributeNames[3];

    std::string xName{ "xt" };
    std::string yName{ "yt" };
    std::string zName{ "zt" };
    std::string tName{ "time" };

    std::string verticalLines = "";
    std::string horizontalLines = "";

    std::vector<float> verticalLineNumbers;
    std::vector<float> horizontalLineNumbers;

    int scatterplotResolution{ 400 };
    size_t histogramResolution{ 300 };

    // Default cropping is no cropping
    size_t xMin{ 0 };
    size_t xMax{ 0 };
    size_t yMin{ 0 };
    size_t yMax{ 0 };
    size_t zMin{ 0 };
    size_t zMax{ 0 };
    size_t tMin{ 0 };
    size_t tMax{ 0 };

    // Default downsample is now downsample
    int downsampleX{ 1 };
    int downsampleY{ 1 };
    int downsampleZ{ 1 };

    void explodeLinesToVector()
    {
        if (verticalLines.length() > 0) {
            std::istringstream input;
            input.str(verticalLines);

            for (std::string number; std::getline(input, number, '|');) {
                verticalLineNumbers.push_back(std::stof(number));
            }
        }
        std::cout << std::endl;
        if (horizontalLines.length() > 0) {
            std::istringstream input;
            input.str(horizontalLines);

            for (std::string number; std::getline(input, number, '|');) {
                horizontalLineNumbers.push_back(std::stof(number));
            }
        }
    }

    std::tuple<size_t, size_t, size_t, size_t> cropDimensions(const int xDim,
                                                              const int yDim,
                                                              const int zDim,
                                                              const int tDim)
    {
        // Set default values if we have 0s
        this->tMax = this->tMax == 0 ? tDim : this->tMax;
        this->xMax = this->xMax == 0 ? xDim : this->xMax;
        this->yMax = this->yMax == 0 ? yDim : this->yMax;
        this->zMax = this->zMax == 0 ? zDim : this->zMax;

        // Make sure all dimensions are set correctly
        assert(this->tMin >= 0 && this->tMax <= tDim && this->tMin <= this->tMax);
        assert(this->xMin >= 0 && this->xMax <= xDim && this->xMin <= this->xMax);
        assert(this->yMin >= 0 && this->yMax <= yDim && this->yMin <= this->yMax);
        assert(this->zMin >= 0 && this->zMax <= zDim && this->zMin <= this->zMax);

        // Return cropped values
        return { this->xMax - this->xMin, this->yMax - this->yMin, this->zMax - this->zMin, this->tMax - this->tMin };
    }
};

namespace config {}
}
