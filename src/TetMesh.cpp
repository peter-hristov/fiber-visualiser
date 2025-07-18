#include "./TetMesh.h"
#include "./CGALTypedefs.h"

#include <random>

void TetMesh::perturbRangeValues(const float &epsilon)
{
    static std::mt19937 gen(std::random_device{}());
    std::uniform_real_distribution<float> dist(-epsilon, epsilon);

    for (size_t i = 0; i < this->vertexCoordinatesF.size(); i++) 
    {
        this->vertexCoordinatesF[i] += dist(gen);
        this->vertexCoordinatesG[i] += dist(gen);
    }
}

void TetMesh::computeDomainBoundingBox()
{
    assert(!vertexDomainCoordinates.empty() && vertexDomainCoordinates[0].size() == 3);

    const auto &first = vertexDomainCoordinates[0];
    minX = maxX = first[0];
    minY = maxY = first[1];
    minZ = maxZ = first[2];

    for (const auto &v : vertexDomainCoordinates)
    {
        minX = std::min(minX, v[0]);
        maxX = std::max(maxX, v[0]);

        minY = std::min(minY, v[1]);
        maxY = std::max(maxY, v[1]);

        minZ = std::min(minZ, v[2]);
        maxZ = std::max(maxZ, v[2]);
    }
}

void TetMesh::computeRangeBoundingBox()
{
    assert(!vertexCoordinatesF.empty() && !vertexCoordinatesG.empty());

    const auto [minFIt, maxFIt] = std::minmax_element(vertexCoordinatesF.begin(), vertexCoordinatesF.end());
    const auto [minGIt, maxGIt] = std::minmax_element(vertexCoordinatesG.begin(), vertexCoordinatesG.end());

    minF = *minFIt;
    maxF = *maxFIt;
    minG = *minGIt;
    maxG = *maxGIt;
}

void TetMesh::computeBoundingBoxes()
{
    computeDomainBoundingBox();
    computeRangeBoundingBox();
}

void TetMesh::sortVertices()
{
    // Save the unsorted index of every points so that we can rename later
    std::vector<std::pair<CartesianPoint, int>> pointWithIndices(this->vertexCoordinatesG.size());

    for (int i = 0 ; i < pointWithIndices.size() ; i++)
    {
        pointWithIndices[i] = {CartesianPoint(this->vertexCoordinatesF[i], this->vertexCoordinatesG[i]), i};
    }

    // Sort all points
    std::sort(pointWithIndices.begin(), pointWithIndices.end(),
            [](const std::pair<CartesianPoint, int> &a, const std::pair<CartesianPoint, int>& b) {
            return CGAL::compare_xy(a.first, b.first) == CGAL::SMALLER;
            });


    // Set up the inverse index search
    std::vector<int> meshIDtoSortIndex(pointWithIndices.size());
    for (int i = 0 ; i < pointWithIndices.size() ; i++)
    {
        meshIDtoSortIndex[pointWithIndices[i].second] = i;
    }


    // Now we can swap things around

    // Set up copies of the originals for the swap, otherwise editin in place causes errors
    std::vector<std::vector<float>> vertexDomainCoordinatesOriginal = this->vertexDomainCoordinates;
    std::vector<float> vertexCoordinatesFOriginal = this->vertexCoordinatesF;
    std::vector<float> vertexCoordinatesGOriginal = this->vertexCoordinatesG;

    // Swap tet indices
    for (int i = 0 ; i < this->tetrahedra.size() ; i++)
    {
        for (int j = 0 ; j < this->tetrahedra[i].size() ; j++)
        {
            this->tetrahedra[i][j] = meshIDtoSortIndex[this->tetrahedra[i][j]];
        }
    }

    // Swap domain positions
    for (int i = 0 ; i < this->vertexDomainCoordinates.size() ; i++)
    {
        this->vertexDomainCoordinates[i] = vertexDomainCoordinatesOriginal[pointWithIndices[i].second];
    }
    
    // Swap range positions
    for (int i = 0 ; i < this->vertexCoordinatesF.size() ; i++)
    {
        this->vertexCoordinatesF[i] = vertexCoordinatesFOriginal[pointWithIndices[i].second];
        this->vertexCoordinatesG[i] = vertexCoordinatesGOriginal[pointWithIndices[i].second];
    }
}

void TetMesh::printMesh()
{
    // Print vertex domain coordinates
    std::cout << "Vertex domain coordinates: " << std::endl;
    for (int i = 0 ; i < this->vertexDomainCoordinates.size() ; i++)
    {
        printf("%d - (%f, %f, %f)\n", i, this->vertexDomainCoordinates[i][0], this->vertexDomainCoordinates[i][1], this->vertexDomainCoordinates[i][2]);
    }
    // Print vertex range coordinates
    std::cout << "Vertex range coordinates: " << std::endl;
    for (int i = 0 ; i < this->vertexDomainCoordinates.size() ; i++)
    {
        printf("%d - (%f, %f)\n", i, this->vertexCoordinatesF[i], this->vertexCoordinatesG[i]);
    }

    // Print tetrahedra
    std::cout << "Tetrahedra: " << std::endl;
    for (int i = 0 ; i < this->tetrahedra.size() ; i++)
    {
        printf("%d - (%ld, %ld, %ld, %ld)\n", i, this->tetrahedra[i][0], this->tetrahedra[i][1], this->tetrahedra[i][2], this->tetrahedra[i][3]);
    }
}


void TetMesh::computeSingularEdgeTypes()
{
    // Initialization
    for (const auto &edge : this->edges)
    {
        this->edgeSingularTypes[edge] = {};
    }


    // For every edge, save the edges in its link
    std::map<std::array<int, 2>, std::vector<std::array<int, 2>>> linkEdges;

    // otherwise the lower and upper link flip around
    for (int i = 0 ; i <  this->tetrahedra.size() ; i++)
    {
        const auto &tet = this->tetrahedra[i];

        // All pairs give you all six edges
        for (int a = 0 ; a < 4 ; a++)
        {
            for (int b = a + 1 ; b < 4 ; b++)
            {
                // Get the indices of the vertices for the edge
                int aIndex = tet[a];
                int bIndex = tet[b];

                // Make sure the vertices of the edge are in sorted order to have consistent orientation
                if (aIndex > bIndex)
                {
                    std::swap(aIndex, bIndex);
                }

                const std::array<int, 2> edge = {aIndex, bIndex};

                std::vector<int> linkVerticesInTet;

                // Search though the other two unused vertices
                for (int v = 0 ; v < 4 ; v++)
                {
                    const int vIndex = tet[v];

                    // If the currect vertex is not one of two used to define the edge
                    if (vIndex != aIndex && vIndex != bIndex) 
                    {
                        linkVerticesInTet.push_back(vIndex);
                    }
                }

                linkEdges[edge].push_back({linkVerticesInTet[0], linkVerticesInTet[1]});
            }
        }
    }

    // Initialize the disjoint set for the upper and lower link 
    std::map<std::array<int, 2>, DisjointSet<int>> upperLinkComponentsDS;
    for (const auto &[edge, vertices] : this->upperLink)
    {
        upperLinkComponentsDS[edge].initialize(vertices);
    }

    std::map<std::array<int, 2>, DisjointSet<int>> lowerLinkComponentsDS;
    for (const auto &[edge, vertices] : this->lowerLink)
    {
        lowerLinkComponentsDS[edge].initialize(vertices);
    }

    // Add all edges to compute the connected components

    for (const auto &[edge, vertices] : this->upperLink)
    {
        for (const std::array<int, 2> &linkEdge : linkEdges[edge])
        {
            if (vertices.contains(linkEdge[0]) && vertices.contains(linkEdge[1]))
            {
                upperLinkComponentsDS[edge].union_setsTriangle(linkEdge[0], linkEdge[1]);
            }
        }
    }

    for (const auto &[edge, vertices] : this->lowerLink)
    {
        for (const std::array<int, 2> &linkEdge : linkEdges[edge])
        {
            if (vertices.contains(linkEdge[0]) && vertices.contains(linkEdge[1]))
            {
                lowerLinkComponentsDS[edge].union_setsTriangle(linkEdge[0], linkEdge[1]);
            }
        }
    }


    for (auto &[edge, type] : this->edgeSingularTypes)
    {
        printf("Currently at the edge [%d, %d].\n", edge[0], edge[1]);

        int upperLinkComponents = upperLinkComponentsDS[edge].countConnectedComponents();
        int lowerLinkComponents = lowerLinkComponentsDS[edge].countConnectedComponents();

        printf("The upper link has %d components and these vertices: ", upperLinkComponents);
        for (const auto &v : this->upperLink[edge])
        {
            printf("%d ", v);
        }
        printf("\n");

        printf("The lower link has %d components and these vertices: ", lowerLinkComponents);
        for (const auto &v : this->lowerLink[edge])
        {
            printf("%d ", v);
        }
        printf("\n");

        // Definite edge
        if (0 == upperLinkComponents || 0 == lowerLinkComponents)
        {
            type = 0;
        }
        // Regular edge
        else if (1 == upperLinkComponents || 1 == lowerLinkComponents)
        {
            type = 1;
        }
        // Indefinite edge
        else
        {
            type = 2;
        }
    }

}
void TetMesh::computeUpperLowerLinkVertices()
{
    // Initialize all maps with the empty set, so that all edges are covered
    for (auto &edge : this->edges)
    {
        this->upperLink[edge] = {};
        this->lowerLink[edge] = {};
        this->upperStarTriangles[edge] = {};
        this->lowerStarTriangles[edge] = {};
    }

    for (int i = 0 ; i <  this->tetrahedra.size() ; i++)
    {
        const auto &tet = this->tetrahedra[i];

        // All pairs give you all six edges
        for (int a = 0 ; a < 4 ; a++)
        {
            for (int b = a + 1 ; b < 4 ; b++)
            {
                // Get the indices of the vertices for the edge
                int aIndex = tet[a];
                int bIndex = tet[b];

                // Make sure the vertices of the edge are in sorted order to have consistent orientation
                if (aIndex > bIndex)
                {
                    std::swap(aIndex, bIndex);
                }

                const std::array<int, 2> edge = {aIndex, bIndex};

                // Search though the other two unused vertices
                for (int v = 0 ; v < 4 ; v++)
                {
                    const int vIndex = tet[v];

                    // If the currect vertex is not one of two used to define the edge
                    if (vIndex != aIndex && vIndex != bIndex) 
                    {
                        const bool isUpperLink = this->isUpperLinkEdgeVertex(aIndex, bIndex, vIndex);

                        if (true == isUpperLink) {
                            this->upperLink[edge].insert(vIndex);
                            this->upperStarTriangles[edge].insert(triangleToIndex.at({aIndex, bIndex, vIndex}));
                        } else {
                            this->lowerLink[edge].insert(vIndex);
                            this->lowerStarTriangles[edge].insert(triangleToIndex.at({aIndex, bIndex, vIndex}));
                        }
                    }
                }
            }
        }
    }
}


bool TetMesh::isUpperLinkEdgeVertex(int aIndex, int bIndex, int vIndex)
{
    // Make sure the vertices of the edge are in sorted order to have consistent orientation
    if (aIndex > bIndex)
    {
        std::swap(aIndex, bIndex);
    }

    // Define the two points that form the line
    const Point_2 a(this->vertexCoordinatesF[aIndex], this->vertexCoordinatesG[aIndex]);
    const Point_2 b(this->vertexCoordinatesF[bIndex], this->vertexCoordinatesG[bIndex]);

    // Define the test point
    const Point_2 v(this->vertexCoordinatesF[vIndex], this->vertexCoordinatesG[vIndex]);  // Change this to test different locations

    // Determine which half-plane r is in
    const CGAL::Orientation result = CGAL::orientation(a, b, v);

    //printf("Checking line (%d, %d) against vertex %d\n", aIndex, bIndex, vIndex);

    //std::cout << "a coords = " << a << std::endl;
    //std::cout << "b coords = " << b << std::endl;
    //std::cout << "v coords = " << v << std::endl;


    // Upper link = left
    if (result == CGAL::LEFT_TURN) {
        return true;
        //std::cout << "Point r is in the LEFT half-plane.\n";
    // Lower link = right
    } else if (result == CGAL::RIGHT_TURN) {
        return false;
        //std::cout << "Point r is in the RIGHT half-plane.\n";
    // This should not happen for generic maps
    } else {
        //std::cout << "Point r is on the line.\n";
        throw std::runtime_error("Input data is degeneate, a triangle is mapped to a line.");       
        //assert(false);
    }

    // Paranoid assert
    assert(false);
}



void TetMesh::computeTriangleAdjacency()
{
    std::set<std::array<int, 2>> allEdges;

    // Initialize all edges types as regular initially
    for (int i = 0 ; i <  this->tetrahedra.size() ; i++)
    {
        const auto &tet = this->tetrahedra[i];

        // All pairs give you all six edges
        for (int a = 0 ; a < 4 ; a++)
        {
            for (int b = a + 1 ; b < 4 ; b++)
            {
                // Get the indices of the vertices for the edge
                int aIndex = tet[a];
                int bIndex = tet[b];

                // Make sure the vertices of the edge are in sorted order to have consistent orientation
                if (aIndex > bIndex)
                {
                    std::swap(aIndex, bIndex);
                }

                // Initialize edge
                allEdges.insert({aIndex, bIndex});
            }
        }
    }

    for (const auto &edge: allEdges)
    {
        this->edges.push_back(edge);
    }




    std::set<std::set<int>> allTriangles;

    for (const std::array<size_t, 4> tet : this->tetrahedra)
    {
        // The triangles of the tet
        std::set<std::set<int>> triangles;

        // All pairs give you all six edges
        for (int a = 0 ; a < 4 ; a++)
        {
            for (int b = a + 1 ; b < 4 ; b++)
            {
                for (int c = b + 1 ; c < 4 ; c++)
                {
                    int aIndex = tet[a];
                    int bIndex = tet[b];
                    int cIndex = tet[c];
                    allTriangles.insert({aIndex, bIndex, cIndex});
                }
            }
        }
    }

    //
    // Set up the indices for all triangles
    //

    for (const std::set<int> triangle : allTriangles)
    {
        this->indexToTriangle.push_back(triangle);
        this->triangleToIndex[triangle] = this->indexToTriangle.size() - 1;
    }

    this->adjacentTrianglesIndex.resize(allTriangles.size());



    // Compute the adjacency of triangles in the mesh, two triangles are adjacent when they are the faces of the same tet
    for (const std::array<size_t, 4> tet : this->tetrahedra)
    {
        // The triangles of the tet
        std::set<std::set<int>> triangles;

        // All pairs give you all six edges
        for (int a = 0 ; a < 4 ; a++)
        {
            for (int b = a + 1 ; b < 4 ; b++)
            {
                for (int c = b + 1 ; c < 4 ; c++)
                {
                    int aIndex = tet[a];
                    int bIndex = tet[b];
                    int cIndex = tet[c];
                    triangles.insert({aIndex, bIndex, cIndex});
                }
            }
        }

        // Connect all the triangles together
        for(const std::set<int> t1 : triangles)
        {
            for(const std::set<int> t2 : triangles)
            {
                // Create a pair of sets
                std::pair<std::set<int>, std::set<int>> pairOfTriangles = {t1, t2};

                // Insert the pair into the set
                this->connectedTriangles.insert(pairOfTriangles);

                int t1Index = this->triangleToIndex[t1];
                int t2Index = this->triangleToIndex[t2];

                this->adjacentTrianglesIndex[t1Index].push_back(t2Index);
                this->adjacentTrianglesIndex[t2Index].push_back(t1Index);
            }
        }
    }
}
