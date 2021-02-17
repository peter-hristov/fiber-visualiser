#include "./MarchingCubes.h"
#include "./Geometry.h"
#include "./SurfaceMesh.h"

#include "./mc_tables.h"

using namespace std;

void
tv9k::utility::MarchingCubes::computeTriangles(GLfloat isovalue,
                                               const vector<vector<vector<GLfloat>>>& vals,
                                               const vector<vector<vector<GLfloat>>>& vals2,
                                               const vector<vector<vector<GLfloat>>>& vals3,
                                               int mode,
                                               tv9k::utility::SurfaceMesh& mesh,
                                               const int isovalueMult,
                                               Data* data)
{
    mesh.clear();

    for (int i = 0; i < data->xdim - 1; i++) {
        for (int j = 0; j < data->ydim - 1; j++) {
            for (int k = 0; k < data->zdim - 1; k++) {
                tv9k::utility::MarchingCubes::prepareCube(
                  i, j, k, isovalue, vals, vals2, vals3, mode, mesh, isovalueMult, data);
            }
        }
    }
}

void
tv9k::utility::MarchingCubes::prepareCube(int i,
                                          int j,
                                          int k,
                                          GLfloat isovalue,
                                          const vector<vector<vector<GLfloat>>>& vals,
                                          const vector<vector<vector<GLfloat>>>& vals2,
                                          const vector<vector<vector<GLfloat>>>& vals3,
                                          int mode,
                                          tv9k::utility::SurfaceMesh& mesh,
                                          const int isovalueMult,
                                          Data* data)
{
    GLfloat values[8];

    values[0] = vals[i][j][k];
    values[1] = vals[i][j][k + 1];
    values[2] = vals[i][j + 1][k];
    values[3] = vals[i][j + 1][k + 1];
    values[4] = vals[i + 1][j][k];
    values[5] = vals[i + 1][j][k + 1];
    values[6] = vals[i + 1][j + 1][k];
    values[7] = vals[i + 1][j + 1][k + 1];

    int lookupValue = 0;

    for (int i = 0; i < 8; i++) {
        lookupValue |= (values[i] < isovalue) << i;
    }

    // 1.9 to 0.9 seconds of speedup right here
    if (lookupValue == 0) {
        return;
    }

    // Keeping or removing these coppies does not affect performance much
    GLfloat values2[8] = { vals2[i][j][k],         vals2[i][j][k + 1],        vals2[i][j + 1][k],
                           vals2[i][j + 1][k + 1], vals2[i + 1][j][k],        vals2[i + 1][j][k + 1],
                           vals2[i + 1][j + 1][k], vals2[i + 1][j + 1][k + 1] };

    GLfloat values3[8] = { vals3[i][j][k],         vals3[i][j][k + 1],        vals3[i][j + 1][k],
                           vals3[i][j + 1][k + 1], vals3[i + 1][j][k],        vals3[i + 1][j][k + 1],
                           vals3[i + 1][j + 1][k], vals3[i + 1][j + 1][k + 1] };

    // this->drawWiredCube(cubeVertices);
    tv9k::utility::MarchingCubes::processCube(
      lookupValue, values, values2, values3, isovalue, i, j, k, mode, mesh, isovalueMult, vals, data);
}

void
tv9k::utility::MarchingCubes::processCube(int lookupValue,
                                          const GLfloat values[8],
                                          const GLfloat values2[8],
                                          const GLfloat values3[8],
                                          GLfloat isovalue,
                                          int x,
                                          int y,
                                          int z,
                                          int mode,
                                          tv9k::utility::SurfaceMesh& mesh,
                                          const int isovalueMult,
                                          const vector<vector<vector<GLfloat>>>& vals,
                                          Data* data)
{
    tv9k::utility::MeshTriangle currentTriangle;

    for (int i = 1; i <= triangleTable[lookupValue][0] * 3; i++) {
        // The indices of the endpoints of the intersected edge as per the array
        // cubeVertices
        int u = edgeTable[triangleTable[lookupValue][i]][0];
        int v = edgeTable[triangleTable[lookupValue][i]][1];

        // Swap if v is smaller than u
        if (values[v] < values[u]) {
            std::swap(v, u);
        }

        // Interpolate
        GLfloat l = (values[v] - isovalue) / (values[v] - values[u]);

        vector<GLfloat> point = {
            l * tv9k::geometry::cubeVertices[u][0] + (1 - l) * tv9k::geometry::cubeVertices[v][0] + x,
            l * tv9k::geometry::cubeVertices[u][1] + (1 - l) * tv9k::geometry::cubeVertices[v][1] + y,
            l * tv9k::geometry::cubeVertices[u][2] + (1 - l) * tv9k::geometry::cubeVertices[v][2] + z,
        };

        if (mode == 0) {
            int ux = x + vertPos[u][0];
            int uy = y + vertPos[u][1];
            int uz = z + vertPos[u][2];

            // Use the bigger vertex if we are looking at superlevel sets
            if (isovalueMult == -1) {
                ux = x + vertPos[v][0];
                uy = y + vertPos[v][1];
                uz = z + vertPos[v][2];
            }

            vector<GLfloat> normal =
              tv9k::geometry::computeCentralDifferencingNormal(vector<GLfloat>{ float(x), float(y), float(z) },
                                                               u,
                                                               v,
                                                               l,
                                                               vals.size(),
                                                               vals[0].size(),
                                                               vals[0][0].size(),
                                                               vals,
                                                               isovalueMult == -1);

            currentTriangle.vertices.push_back(point);
            currentTriangle.normals.push_back(normal);

            // This is for the 2D scatterplot
            GLfloat a = l * values2[u] + (1 - l) * values2[v];
            GLfloat b = l * values3[u] + (1 - l) * values3[v];
            // scatterplotTriangles.push_back(make_tuple(a, b, visited[ux][uy][uz]));

            currentTriangle.projectedVertices.push_back({ a, b });

            // Save a triangle every 3 vertices
            if (i % 3 == 0) {
                currentTriangle.triangleId = mesh.visited[ux][uy][uz];
                currentTriangle.color = Utility::getColor(mesh.visited[ux][uy][uz], 0, 255);
                mesh.triangles.push_back(currentTriangle);
                currentTriangle.clear();
            }

        } else if (mode == 1) {
            // Integer Coordinates of the smaller endpoint on the edge
            int ux = x + vertPos[u][0];
            int uy = y + vertPos[u][1];
            int uz = z + vertPos[u][2];

            vector<GLfloat> normal =
              tv9k::geometry::computeCentralDifferencingNormal(vector<GLfloat>{ float(x), float(y), float(z) },
                                                               u,
                                                               v,
                                                               l,
                                                               vals.size(),
                                                               vals[0].size(),
                                                               vals[0][0].size(),
                                                               vals,
                                                               false);

            currentTriangle.vertices.push_back(point);
            currentTriangle.normals.push_back(normal);

            // Save a triangle every 3 vertices
            if (i % 3 == 0) {
                currentTriangle.triangleId = mesh.visited[ux][uy][uz];
                currentTriangle.color = Utility::getColor(mesh.visited[ux][uy][uz], 1, 255);
                mesh.triangles.push_back(currentTriangle);
                currentTriangle.clear();
            }
            // Combine Surface
        } else if (mode == 2) {
            // Integer Coordinates of the smaller endpoint on the edge
            int ux = x + vertPos[u][0];
            int uy = y + vertPos[u][1];
            int uz = z + vertPos[u][2];

            vector<GLfloat> normal =
              tv9k::geometry::computeCentralDifferencingNormal(vector<GLfloat>{ float(x), float(y), float(z) },
                                                               u,
                                                               v,
                                                               l,
                                                               vals.size(),
                                                               vals[0].size(),
                                                               vals[0][0].size(),
                                                               vals,
                                                               false);

            currentTriangle.vertices.push_back(point);
            currentTriangle.normals.push_back(normal);

            // Save a triangle every 3 vertices
            if (i % 3 == 0) {
                currentTriangle.triangleId = mesh.visited[ux][uy][uz];
                currentTriangle.color = Utility::getColor(mesh.visited[ux][uy][uz], 2, 255);
                mesh.triangles.push_back(currentTriangle);
                currentTriangle.clear();
            }
        } else {
            assert(false);
        }
    }
}
