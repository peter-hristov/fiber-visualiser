#include <vector>
#include <utility>
#include <string>

namespace utility
{

    //
    // Take half a tophat of the form (v0, V4, v3, V2, v1). We only use half because the other part is mirrored
    // It then arranges the points in the order around a circle.
    //
    std::pair<std::vector<float>, std::vector<float>> generateCoordinatesFromTophat(const std::vector<std::string> &tophat);
    
}
