#include "./Timer.h"
#include "./Arrangement.h"


void Arrangement::computeArrangement(const TetMesh &tetMesh) 
{
    // Add in the vertices of the mesh 
    this->arrangementPoints.resize(tetMesh.vertexCoordinatesF.size());
    for (int i = 0 ; i < tetMesh.vertexCoordinatesF.size() ; i++)
    {
        const float u = tetMesh.vertexCoordinatesF[i];
        const float v = tetMesh.vertexCoordinatesG[i];
        const Point_2 point(u, v);

        this->arrangementPoints[i] = point;
        this->arrangementPointIndices[point] = i;
    };

    // Add the unique edges as setments to the arrangement
    std::vector<Segment_2> segments;
    for (const auto& edge : tetMesh.edges) 
    {
        segments.push_back(Segment_2(this->arrangementPoints[edge[0]], this->arrangementPoints[edge[1]]));
    }

    Timer::start();
    CGAL::insert(this->arr, segments.begin(), segments.end());
    Timer::stop("Computed arrangement                   :");

    // Sequential arrangement computationa
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
