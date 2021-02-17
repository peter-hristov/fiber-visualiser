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

    int scatterplotResolution{ 400 };

    std::string verticalLines = "";
    std::string horizontalLines = "";

    std::vector<float> verticalLineNumbers;
    std::vector<float> horizontalLineNumbers;

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
};

namespace config {}
}
