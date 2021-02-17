#include "./SurfaceMesh.h"

using namespace std;

void
tv9k::utility::SurfaceMesh::clear()
{
    this->triangles.clear();
    // this->visited.clear();
}

std::pair<std::vector<GLfloat>, std::vector<GLfloat>>
tv9k::utility::SurfaceMesh::getMinMaxPointsObject(const int objectId)
{
    // Value of the extramal corners of the bouding box of the object
    vector<GLfloat> minPoint{ 0, 0, 0 };
    vector<GLfloat> maxPoint{ 0, 0, 0 };

    // Used to set a default value for the min and max point as the first found
    // point;
    bool firstFound = false;

    for (int i = 0; i < this->triangles.size(); i++) {
        // Do not consider other objects than the one with the desired id
        if (this->triangles[i].triangleId != objectId) {
            continue;
        }

        for (int j = 0; j < 3; j++) {
            vector<GLfloat> vertex = this->triangles[i].vertices[j];

            if (firstFound == false) {
                minPoint = vertex;
                maxPoint = vertex;
                firstFound = true;
            } else {
                if (minPoint[0] > vertex[0]) {
                    minPoint[0] = vertex[0];
                }
                if (minPoint[1] > vertex[1]) {
                    minPoint[1] = vertex[1];
                }
                if (minPoint[2] > vertex[2]) {
                    minPoint[2] = vertex[2];
                }

                if (maxPoint[0] < vertex[0]) {
                    maxPoint[0] = vertex[0];
                }
                if (maxPoint[1] < vertex[1]) {
                    maxPoint[1] = vertex[1];
                }
                if (maxPoint[2] < vertex[2]) {
                    maxPoint[2] = vertex[2];
                }
            }
        }
    }
    return { minPoint, maxPoint };
}

std::vector<GLfloat>
tv9k::utility::SurfaceMesh::getCenter(const int objectId)
{
    // vector<GLfloat> minPoint;
    // vector<GLfloat> maxPoint;

    // if (true == firstFound)
    //{
    //// Return the midway point on the bounding box diagonal
    // return {(minPoint[0] + maxPoint[0]) / 2, (minPoint[1] + maxPoint[1]) / 2,
    // (minPoint[2] + maxPoint[2]) / 2};
    //}
    // else
    //{
    // return {0, 0, 0};
    //}
}
