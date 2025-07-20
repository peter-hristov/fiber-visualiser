#include "./Timer.h"
#include "./Arrangement.h"


void Arrangement::computeArrangement(const TetMesh &tetMesh) 
{
    Timer::start();

    // Add in the vertices of the mesh 
    this->arrangementPoints.resize(tetMesh.vertexCoordinatesF.size());
    for (int i = 0 ; i < tetMesh.vertexCoordinatesF.size() ; i++)
    {
        const float u = tetMesh.vertexCoordinatesF[i];
        const float v = tetMesh.vertexCoordinatesG[i];
        const Point_2 point(u, v);

        this->arrangementPoints[i] = point;
        this->arrangementPointsIdices[point] = i;
    };
    Timer::stop("Converted vertices to points           :");

    Timer::start();
    // Make sure you don't add duplicate edge to the arrangement
    std::map<std::set<int>, bool> uniqueEdges;
    for (int i = 0 ; i < tetMesh.tetrahedra.size() ; i++)
    {
        // Bottom face
        uniqueEdges[std::set<int>({tetMesh.tetrahedra[i][0], tetMesh.tetrahedra[i][1]})] = true;
        uniqueEdges[std::set<int>({tetMesh.tetrahedra[i][1], tetMesh.tetrahedra[i][2]})] = true;
        uniqueEdges[std::set<int>({tetMesh.tetrahedra[i][2], tetMesh.tetrahedra[i][0]})] = true;

        // Connect bottom face to top vertex
        uniqueEdges[std::set<int>({tetMesh.tetrahedra[i][0], tetMesh.tetrahedra[i][3]})] = true;
        uniqueEdges[std::set<int>({tetMesh.tetrahedra[i][1], tetMesh.tetrahedra[i][3]})] = true;
        uniqueEdges[std::set<int>({tetMesh.tetrahedra[i][2], tetMesh.tetrahedra[i][3]})] = true;
    }
    Timer::stop("Computed unique edges                  :");

    Timer::start();
    // Add the unique edges as setments to the arrangement
    std::vector<Segment_2> segments;
    for (const auto& edge : uniqueEdges) 
    {
        // Put in a vector for easy access
        std::vector<int> edgeVector(edge.first.begin(), edge.first.end());
        assert(edgeVector.size() == 2);

        segments.push_back(Segment_2(this->arrangementPoints[edgeVector[0]], this->arrangementPoints[edgeVector[1]]));
        //std::cout << "Adding edge " << edgeVector[0] << " - " << edgeVector[1] << std::endl;
    }
    Timer::stop("Converted edges to segments            :");



    //Timer::start();
    //int intersections = 0;
    //for (int i = 0 ; i < segments.size(); i++)
    //{
        //for (int j = i+1 ; j < segments.size(); j++)
        //{
            //// Check for intersection
            //if (CGAL::do_intersect(segments[i], segments[j])) 
            //{
                //intersections++;
            //}
        //}
    //}

    //std::cout << "There are " << intersections << " out of " << segments.size() << " segments." << std::endl;
    //Timer::stop("Checking intersections                 :");


    Timer::start();
    // Insert all the segments into the arrangement using the range insert method
    CGAL::insert(this->arr, segments.begin(), segments.end());
    Timer::stop("Computed arrangement                   :");

    //Timer::start();
    //Arrangement_2 arr;
    //for (const auto& segment : segments) {
        //CGAL::insert(arr, Curve_2(segment));
    //}
    //Timer::stop("Computed arrangement sequantially      :");


    std::cout << "The arrangement size:"
        << "   |V| = " << this->arr.number_of_vertices()
        << ",  |E| = " << this->arr.number_of_edges()
        << ",  |F| = " << this->arr.number_of_faces() << std::endl << std::endl;

    //std::cout << std::endl << std::endl << "The sequantial arrangement size:\n"
        //<< "   |V| = " << arr.number_of_vertices()
        //<< ",  |E| = " << arr.number_of_edges()
        //<< ",  |F| = " << arr.number_of_faces() << std::endl;

    // Set up the indices and their reverse lookup for all faces
    int counter = 0;
    this->arrangementIndexToFace.resize(this->arr.number_of_faces());
    for (auto f = this->arr.faces_begin(); f != this->arr.faces_end(); ++f) 
    {
        Arrangement_2::Face_const_handle a = f;
        this->arrangementFacesIdices[a] = counter;
        this->arrangementIndexToFace[counter] = f;
        counter++;
    }
}

void Arrangement::computePointLocationDataStructure()
{
    this->pl = std::make_unique<Point_location>(this->arr);
}
