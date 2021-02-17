#pragma once

#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif
#include <vector>

namespace tk9k {
namespace utility {
namespace histogram {
std::vector<unsigned int>
computeHistogram(size_t numberOfBins,
                 const GLfloat minValue,
                 const GLfloat maxValue,
                 const std::vector<std::vector<std::vector<GLfloat>>>& data)
{
    std::vector<unsigned int> histogram(numberOfBins, 0);

    for (int i = 0; i < data.size(); i++) {
        for (int j = 0; j < data[0].size(); j++) {
            for (int k = 0; k < data[0][0].size(); k++) {
                if (data[i][j][k] != 0.0) {
                    // Rescale data point
                    size_t rescaledValue = (numberOfBins / (maxValue - minValue)) * (data[i][j][k] - minValue);

                    // Clip to stay inside the histogram
                    rescaledValue = std::max(static_cast<size_t>(0), rescaledValue);
                    rescaledValue = std::min(rescaledValue, numberOfBins - 1);

                    histogram[rescaledValue]++;
                }
            }
        }
    }

    return histogram;
}
}
}
}
