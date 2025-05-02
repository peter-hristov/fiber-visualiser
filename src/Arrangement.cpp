#include "./Arrangement.h"
#include "src/CustomArrangementTypes.h"

#include <vector>

void CustomArrangement::computeArrangementCustom(Data *data)
{
    //
    // Make an array of points,
    //

    // ASSUMPTION - these are already sorted lexicographically
    std::vector<PointHandler> points(data->vertexCoordinatesF.size());

    for (int i = 0 ; i < data->vertexCoordinatesF.size() ; i++)
    {
        const double u = data->vertexCoordinatesF[i];
        const double v = data->vertexCoordinatesG[i];

        points[i] = std::make_shared<Point>(Point(u, v, i));
    }


    //
    // Make an array of segments
    //

    // Get the unique edges
    std::map<std::set<size_t>, bool> uniqueEdges;
    for (int i = 0 ; i < data->tetrahedra.size() ; i++)
    {
        // Bottom face
        uniqueEdges[std::set<size_t>({data->tetrahedra[i][0], data->tetrahedra[i][1]})] = true;
        uniqueEdges[std::set<size_t>({data->tetrahedra[i][1], data->tetrahedra[i][2]})] = true;
        uniqueEdges[std::set<size_t>({data->tetrahedra[i][2], data->tetrahedra[i][0]})] = true;

        // Connect bottom face to top vertex
        uniqueEdges[std::set<size_t>({data->tetrahedra[i][0], data->tetrahedra[i][3]})] = true;
        uniqueEdges[std::set<size_t>({data->tetrahedra[i][1], data->tetrahedra[i][3]})] = true;
        uniqueEdges[std::set<size_t>({data->tetrahedra[i][2], data->tetrahedra[i][3]})] = true;
    }

    std::vector<SegmentHandler> segments;
    for (const auto& edge : uniqueEdges) 
    {
        // Put in a vector for easy access
        std::vector<int> edgeVector(edge.first.begin(), edge.first.end());
        assert(edgeVector.size() == 2);

        // Swap if the first point is not smaller
        if (*points[edgeVector[1]] < *points[edgeVector[0]])
        {
            std::swap(edgeVector[0], edgeVector[1]);

        }

        segments.push_back(std::make_shared<Segment>(points[edgeVector[0]], points[edgeVector[1]], segments.size()));
    }

    std::cout << "Printing all points...\n";

    for (int i = 0 ; i < points.size() ; i++)
    {
        printf("Point %d (%.2f, %.2f).\n", points[i]->index, points[i]->x, points[i]->y);
    }


    std::cout << "\n\nPrinting all segments...\n";

    for (int i = 0 ; i < segments.size() ; i++)
    {
        auto pointA = segments[i]->a;
        auto pointB = segments[i]->b;
        printf("Segment %d between point %d (%.2f, %.2f) and point %d (%.2f, %.2f).\n", segments[i]->index, pointA->index, pointA->x, pointA->y, pointB->index, pointB->x, pointB->y);
    }

    data->customArrangementPoints = points;
    data->customArrangementSegments = segments;

    for (int i = 0 ; i < segments.size() ; i++)
    {
        for (int j = i + 1 ; j < segments.size() ; j++)
        {
            printf("Checking segments %d and %d.\n", i, j);
            std::optional<Point> intersectionPoint = computeIntersection(segments[i], segments[j]);


            // Add the new point of intersection to the list of all arrangement points
            if (intersectionPoint) 
            {
                std::cout << "Intersection at (" << intersectionPoint->x << ", " << intersectionPoint->y << ")\n";

                // Set the index of the point
                intersectionPoint->index = data->customArrangementPoints.size();
                data->customArrangementPoints.push_back(std::make_shared<Point>(*intersectionPoint));

            } else 
            {
                std::cout << "No intersection.\n";
            }

        }

    }
}



std::optional<CustomArrangement::Point> CustomArrangement::computeIntersection(ConstSegmentHandler &s1, ConstSegmentHandler &s2)
{
    const ConstPointHandler A = s1->a;
    const ConstPointHandler B = s1->b;
    const ConstPointHandler C = s2->a;
    const ConstPointHandler D = s2->b;

    const CustomArrangement::DataType dx1 = B->x - A->x;
    const CustomArrangement::DataType dy1 = B->y - A->y;
    const CustomArrangement::DataType dx2 = D->x - C->x;
    const CustomArrangement::DataType dy2 = D->y - C->y;

    const CustomArrangement::DataType denom = dx1 * dy2 - dy1 * dx2;

    // Lines are parallel or collinear
    if (denom == 0) 
    {
        return std::nullopt;
    }

    const CustomArrangement::DataType dx = C->x - A->x;
    const CustomArrangement::DataType dy = C->y - A->y;

    const CustomArrangement::DataType t_num = dx * dy2 - dy * dx2;
    const CustomArrangement::DataType u_num = dx * dy1 - dy * dx1;

    const CustomArrangement::DataType t = t_num / denom;
    const CustomArrangement::DataType u = u_num / denom;

    if (t >= 0.0 && t <= 1.0 && u >= 0.0 && u <= 1.0) 
    {
        const CustomArrangement::DataType ix = A->x + t * dx1;
        const CustomArrangement::DataType iy = A->y + t * dy1;
        return Point(ix, iy, -1);  // -1 since it's not from input, we'll set it later
    }

    return std::nullopt;
}
