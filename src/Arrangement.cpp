#include "./Arrangement.h"

#include <vector>

void CustomArrangement::computeArrangementCustom(Data *data)
{
    //
    // Make an array of points,
    //

    // ASSUMPTION - these are already sorted lexicographically
    std::vector<std::shared_ptr<Point<double>>> points(data->vertexCoordinatesF.size());

    for (int i = 0 ; i < data->vertexCoordinatesF.size() ; i++)
    {
        const double u = data->vertexCoordinatesF[i];
        const double v = data->vertexCoordinatesG[i];

        points[i] = std::make_shared<Point<double>>(Point<double>(u, v, i));
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

    std::vector<std::shared_ptr<Segment<double>>> segments;
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

        segments.push_back(std::make_shared<Segment<double>>(points[edgeVector[0]], points[edgeVector[1]], segments.size()));
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
            std::optional<Point<double>> intersectionPoint = computeIntersection<double>(segments[i], segments[j]);


            // Add the new point of intersection to the list of all arrangement points
            if (intersectionPoint) 
            {
                std::cout << "Intersection at (" << intersectionPoint->x << ", " << intersectionPoint->y << ")\n";

                // Set the index of the point
                intersectionPoint->index = data->customArrangementPoints.size();
                data->customArrangementPoints.push_back(std::make_shared<Point<double>>(*intersectionPoint));

            } else 
            {
                std::cout << "No intersection.\n";
            }

        }

    }
}



template <typename coordType>
std::optional<CustomArrangement::Point<coordType>> CustomArrangement::computeIntersection(const std::shared_ptr<const CustomArrangement::Segment<coordType>> &s1, const std::shared_ptr<const CustomArrangement::Segment<coordType>> &s2)
{
    const auto A = s1->a;
    const auto B = s1->b;
    const auto C = s2->a;
    const auto D = s2->b;

    const coordType dx1 = B->x - A->x;
    const coordType dy1 = B->y - A->y;
    const coordType dx2 = D->x - C->x;
    const coordType dy2 = D->y - C->y;

    const coordType denom = dx1 * dy2 - dy1 * dx2;

    // Lines are parallel or collinear
    if (denom == 0) 
    {
        return std::nullopt;
    }

    const coordType dx = C->x - A->x;
    const coordType dy = C->y - A->y;

    const coordType t_num = dx * dy2 - dy * dx2;
    const coordType u_num = dx * dy1 - dy * dx1;

    const coordType t = t_num / denom;
    const coordType u = u_num / denom;

    if (t >= 0.0 && t <= 1.0 && u >= 0.0 && u <= 1.0) 
    {
        const coordType ix = A->x + t * dx1;
        const coordType iy = A->y + t * dy1;
        return Point<coordType>(ix, iy, -1);  // -1 since it's not from input, we'll set it later
    }

    return std::nullopt;
}

