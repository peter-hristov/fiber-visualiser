#include "./utility.h"
#include <cmath>

using namespace std;

std::pair<std::vector<float>, std::vector<float>> utility::generateCoordinatesFromTophat(const std::vector<std::string> &tophat)
{
    float step = (2 * M_PI) / (tophat.size() * 2);

    vector<float> circleOrder(tophat.size(), -1);

    auto isCapital = [](const char &c) { return (c >= 'A' && c <= 'Z'); };

    for (int i = 0 ; i < tophat.size(); i++)
    {
        char x = tophat[i].at(0);
        int id = stoi(tophat[i].substr(1));

        if (circleOrder[id] != -1)
        {
            continue;
        }

        if (isCapital(x))
        {
            circleOrder[id] = i * step + M_PI;
        }
        else
        {
            circleOrder[id] = i * step;
        }
    }


    vector<float> coordinatesF, coordinatesG;
    for (int i = 0 ; i < circleOrder.size() ; i++)
    {
        float t = circleOrder[i];
        t -= M_PI / 2;

        coordinatesF.push_back(cos(t));
        coordinatesG.push_back(-sin(t));
    }
    coordinatesF.push_back(0.0);
    coordinatesG.push_back(0.0);

    return {coordinatesF, coordinatesG};
}
