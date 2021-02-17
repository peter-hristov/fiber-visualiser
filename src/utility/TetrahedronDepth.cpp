#include "TetrahedronDepth.h"

using namespace std;

float
TetrahedronDepth::dotProduct(vector<float> a, vector<float> b)
{
    float sum = 0;

    for (int i = 0; i < a.size(); i++) {
        sum += a[i] * b[i];
    }

    return sum;
}

float
TetrahedronDepth::getLength(vector<float> u)
{
    return sqrt(dotProduct(u, u));
}

vector<float>
TetrahedronDepth::get2DIntersection(vector<float> a, vector<float> b, vector<float> c, vector<float> d)
{
    float m1 = (a[1] - b[1]) / (a[0] - b[0]);
    float m2 = (c[1] - d[1]) / (c[0] - d[0]);

    float b1 = (a[1] - m1 * a[0]);
    float b2 = (c[1] - m2 * c[0]);

    float x = (b2 - b1) / (m1 - m2);
    float y = m1 * x + b1;

    return vector<float>{ x, y, 0 };
}

void
TetrahedronDepth::printVector(vector<float> a)
{
    if (a.size() == 0) {
        return;
    }

    cout << "{";

    for (int i = 0; i < a.size() - 1; i++) {
        cout << a[i] << ", ";
    }

    cout << a[a.size() - 1];

    cout << "}" << endl;
}

vector<float>
TetrahedronDepth::getVectorFromPoints(vector<float> a, vector<float> b)
{
    vector<float> c;

    for (int i = 0; i < a.size(); i++) {
        c.push_back(b[i] - a[i]);
    }

    return c;
}

vector<float>
TetrahedronDepth::get2DNormal(vector<float> a)
{
    return { -a[1], a[0] };
}

vector<int>
TetrahedronDepth::findCrossing(vector<vector<float>> points)
{
    for (int v = 1; v < 4; v++) {
        vector<float> normal = get2DNormal(getVectorFromPoints(points[0], points[v]));

        vector<int> other = { 1, 2, 3 };

        // Remove the vertices that are used for the normal (0 and v)
        other.erase(other.begin() + v - 1);

        float dot1 = dotProduct(getVectorFromPoints(points[0], points[other[0]]), normal);
        float dot2 = dotProduct(getVectorFromPoints(points[0], points[other[1]]), normal);

        // If the dot product with the two other vertices has oposite signs then
        // their product is less than 0 so BINGO
        if (dot1 * dot2 < 0) {
            return { 0, v, other[0], other[1] };
        }
    }

    // We would only reach this point if something went wrong.
    return { -1, -1, -1, -1 };
}

bool
TetrahedronDepth::isInsideTriangle(vector<float> a, vector<float> b, vector<float> c, vector<float> x)
{
    vector<float> ab = get2DNormal(getVectorFromPoints(a, b));
    vector<float> bc = get2DNormal(getVectorFromPoints(b, c));
    vector<float> ca = get2DNormal(getVectorFromPoints(c, a));

    vector<float> ax = getVectorFromPoints(a, x);
    vector<float> bx = getVectorFromPoints(b, x);
    vector<float> cx = getVectorFromPoints(c, x);

    vector<float> dotProducts = { dotProduct(ab, ax), dotProduct(bc, bx), dotProduct(ca, cx) };

    if (dotProducts[0] <= 0 && dotProducts[1] <= 0 && dotProducts[2] <= 0) {
        return true;
    }

    if (dotProducts[0] >= 0 && dotProducts[1] >= 0 && dotProducts[2] >= 0) {
        return true;
    }

    return false;
}

vector<float>
TetrahedronDepth::carthesianToBarycentric(vector<vector<float>> points, vector<int> outside, int inside)
{
    // Aliases for the points to make the following paragraphs less messy
    vector<float> A = points[outside[0]];
    vector<float> B = points[outside[1]];
    vector<float> C = points[outside[2]];
    vector<float> D = points[inside];

    // What we are going to divide by every elemnt in the following matrix
    float divideBy = (A[0] - C[0]) * (B[1] - C[1]) - (B[0] - C[0]) * (A[1] - C[1]);

    // Inverse of the Carthesian to Barycentric Matrix.
    vector<vector<float>> T(3, vector<float>(3));
    T[0][0] = (B[1] - C[1]) / divideBy;
    T[1][0] = (C[1] - A[1]) / divideBy;
    T[0][1] = (C[0] - B[0]) / divideBy;
    T[1][1] = (A[0] - C[0]) / divideBy;

    float u = T[0][0] * (D[0] - C[0]) + T[0][1] * (D[1] - C[1]);
    float v = T[1][0] * (D[0] - C[0]) + T[1][1] * (D[1] - C[1]);

    return { u, v, 1 - u - v };
}

float
TetrahedronDepth::getDepth(vector<vector<float>> points, vector<vector<float>> p, bool debug)
{
    int inside = -1;

    // Determine Case by checking if any of the points is inside the triangle
    if (isInsideTriangle(p[0], p[1], p[2], p[3])) {
        inside = 3;
    } else if (isInsideTriangle(p[0], p[1], p[3], p[2])) {
        inside = 2;
    } else if (isInsideTriangle(p[0], p[2], p[3], p[1])) {
        inside = 1;
    } else if (isInsideTriangle(p[1], p[2], p[3], p[0])) {
        inside = 0;
    }

    if (debug) {
        cout << "The point inside is " << inside << endl;
    }

    // Case 1 - One of the points is inside the triangle
    if (inside != -1) {
        // cout << "We are in CASE 1!!!" << endl;

        // Figure out the indices of the points which are outisde
        vector<int> outside;
        for (int i = 0; i < p.size(); i++) {
            if (i != inside) {
                outside.push_back(i);
            }
        }

        // Obtain the barycentric coordinates of the point inside
        vector<float> baryCoords = carthesianToBarycentric(p, outside, inside);

        // Use the same barycentric coordinates to get the point on the bottom
        vector<float> bottomPoint = { baryCoords[0] * points[outside[0]][0] + baryCoords[1] * points[outside[1]][0] +
                                        baryCoords[2] * points[outside[2]][0],
                                      baryCoords[0] * points[outside[0]][1] + baryCoords[1] * points[outside[1]][1] +
                                        baryCoords[2] * points[outside[2]][1],
                                      baryCoords[0] * points[outside[0]][2] + baryCoords[1] * points[outside[1]][2] +
                                        baryCoords[2] * points[outside[2]][2] };

        vector<float> depthVector = { points[inside][0] - bottomPoint[0],
                                      points[inside][1] - bottomPoint[1],
                                      points[inside][2] - bottomPoint[2] };

        if (debug) {
            cout << "The barycentric coordinates are : " << baryCoords[0] << " " << baryCoords[1] << " "
                 << baryCoords[2] << endl;

            cout << "The bottom point is : ";
            printVector(bottomPoint);

            cout << "The depth vector is : ";
            printVector(depthVector);
            cout << endl << "It has length : " << getLength(depthVector) << endl;
        }

        return getLength(depthVector);
    } else {

        // cout << "We are in CASE 2!!!" << endl;

        // Find which edges are crossed
        vector<int> cross = findCrossing(p);

        // Find the intersection point
        vector<float> i = get2DIntersection(p[cross[0]], p[cross[1]], p[cross[2]], p[cross[3]]);

        vector<float> e1 = getVectorFromPoints(p[cross[0]], p[cross[1]]);
        vector<float> e2 = getVectorFromPoints(p[cross[2]], p[cross[3]]);

        float ratio1 = getLength(getVectorFromPoints(p[cross[0]], i)) / getLength(e1);
        float ratio2 = getLength(getVectorFromPoints(p[cross[2]], i)) / getLength(e2);

        vector<float> E1 = getVectorFromPoints(points[cross[0]], points[cross[1]]);
        vector<float> E2 = getVectorFromPoints(points[cross[2]], points[cross[3]]);

        // Points on the edges as translted vectors
        vector<float> c1 = { ratio1 * E1[0] + points[cross[0]][0],
                             ratio1 * E1[1] + points[cross[0]][1],
                             ratio1 * E1[2] + points[cross[0]][2] };
        vector<float> c2 = { ratio2 * E2[0] + points[cross[2]][0],
                             ratio2 * E2[1] + points[cross[2]][1],
                             ratio2 * E2[2] + points[cross[2]][2] };

        vector<float> depthVector = { c2[0] - c1[0], c2[1] - c1[1], c2[2] - c1[2] };

        if (debug) {
            cout << "Here are the crossings : " << cross[0] << " " << cross[1] << " " << cross[2] << " " << cross[3]
                 << endl;
            cout << "The crossing point is ";
            printVector(i);
            cout << endl;

            cout << "e1 is ";
            printVector(e1);
            cout << "e2 is ";
            printVector(e2);
            cout << endl;

            cout << "The ratios are : " << ratio1 << " " << ratio2 << endl;

            cout << "E1 is ";
            printVector(E1);

            cout << "E2 is ";
            printVector(E2);
            cout << endl;

            cout << "c1 is ";
            printVector(c1);

            cout << "c2 is ";
            printVector(c2);
            cout << endl;

            cout << "The depth vector is : ";
            printVector(depthVector);
            cout << endl << "It has length : " << getLength(depthVector) << endl;
        }

        return getLength(depthVector);
    }
}

vector<float>
TetrahedronDepth::getRandomPoint(int dimension, float min, float max)
{
    vector<float> point;

    for (int i = 0; i < dimension; i++) {
        float a = (float)rand() / (float)RAND_MAX * (max - min) + min;

        point.push_back(a);
    }

    return point;
}

void
TetrahedronDepth::runPermutationsTest()
{
    srand((unsigned)time(0));

    vector<vector<float>> points(4);
    vector<vector<float>> p(4);

    int indices[] = { 0, 1, 2, 3 };

    for (int i = 0; i < 100; i++) {
        // cout << "Here are the points: " << endl << endl;

        for (int j = 0; j < 4; j++) {
            points[j] = getRandomPoint(3, -10000, 10000);
            p[j] = getRandomPoint(2, -10000, 10000);

            printVector(points[j]);
            printVector(p[j]);
            cout << endl;
        }

        vector<float> depths;
        do {
            vector<vector<float>> points2(4);
            vector<vector<float>> p2(4);

            points2[0] = points[indices[0]];
            points2[1] = points[indices[1]];
            points2[2] = points[indices[2]];
            points2[3] = points[indices[3]];

            p2[0] = p[indices[0]];
            p2[1] = p[indices[1]];
            p2[2] = p[indices[2]];
            p2[3] = p[indices[3]];

            depths.push_back(getDepth(points2, p2, false));
        } while (std::next_permutation(indices, indices + 4));

        for (int j = 0; j < depths.size(); j++) {
            // cout << "The depth is - " << depths[j] << endl;
        }

        for (int j = 0; j < depths.size() - 1; j++) {
            // if (depths[j] != depths[j + 1])
            // if (abs(depths[j] - depths[j + 1] >
            // std::numeric_limits<float>::epsilon()))
            if (abs(depths[j] - depths[j + 1] > 0.001)) {
                cout << endl
                     << endl
                     << setprecision(10)
                     << "---------------------------------------- DOES NOT PASS "
                        "----------------------------------------------- "
                     << depths[j] << " " << depths[j + 1] << endl
                     << endl;
            }
        }
    }
}
