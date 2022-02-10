#pragma once

#include <vector>


// The location of a fiber on a single tet
class FaceFiberPoint{

    public:

    std::vector<float> point;
    int colour;

    // The barycentric coordinates as well as the three triangle vertices
    FaceFiberPoint(const float alpha, const float beta, const std::vector<std::vector<float>> vertices)
    {
        this->point = {0, 0, 0};

        point[0] = 
            alpha * vertices[0][0] +
            beta * vertices[1][0] +
            (1 - alpha - beta) * vertices[2][0];

        point[1] = 
            alpha * vertices[0][1] +
            beta * vertices[1][1] +
            (1 - alpha - beta) * vertices[2][1];

        point[2] = 
            alpha * vertices[0][2] +
            beta * vertices[1][2] +
            (1 - alpha - beta) * vertices[2][2];
    }
};
