#include "./Timer.h"
#include "./Arrangement.h"

Face_const_handle Arrangement::getActiveFace(const std::array<float, 2> fiberPoint)
{
    // The query point (u, v)
    Point_2 query_point(fiberPoint[0], fiberPoint[1]);

    // Locate the point in the arrangement
    CGAL::Object result = pl->locate(query_point);

    // Try to assign to a face, edge or a vertex
    Arrangement_2::Face_const_handle face;
    Arrangement_2::Halfedge_const_handle edge;
    Arrangement_2::Vertex_const_handle vertex;

    if (CGAL::assign(face, result)) 
    {
    } 
    // If we are on an edge, just grad an adjacent face
    else if (CGAL::assign(edge, result)) 
    {
        face = edge->face();
    } 
    // If we are on a vertex grab an indicent edge and get its face
    else if (CGAL::assign(vertex, result)) 
    {
        edge = vertex->incident_halfedges();
        face = edge->face();
    } else 
    {
        assert(false);
    }

    return face;
}

void Arrangement::computeArrangement(const TetMesh &tetMesh, const SegmentMode &segmentMode) 
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
    //for (const auto& edge : tetMesh.edges) 
    for (const auto &[edge, type] : tetMesh.edgeSingularTypes) 
    {
        if (segmentMode == SegmentMode::UseSingularSegments)
        {
            if (type != 1)
            {
                segments.push_back(Segment_2(this->arrangementPoints[edge[0]], this->arrangementPoints[edge[1]]));
            }

        }
        else
        {
            segments.push_back(Segment_2(this->arrangementPoints[edge[0]], this->arrangementPoints[edge[1]]));
        }
    }

    CGAL::insert(this->arr, segments.begin(), segments.end());


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







    //for (auto currentFaceIterator = arr.faces_begin(); currentFaceIterator != arr.faces_end(); ++currentFaceIterator) 
    //{
        //Arrangement_2::Face_const_handle face = currentFaceIterator;
        //if (face->is_unbounded()) { continue; }

        //Arrangement_2::Ccb_halfedge_const_circulator start = face->outer_ccb();
        //Arrangement_2::Ccb_halfedge_const_circulator curr = start;

        //do {
            //// Make sure there is only one originating curve (sanity check)
            //const Segment_2 &segment = *arr.originating_curves_begin(curr);

            //const int aIndex = arrangementPointIndices.at(segment.source());
            //const int bIndex = arrangementPointIndices.at(segment.target());

            //if (tetMesh.edgeSingularTypes.at({aIndex, bIndex}) != 1)
            //{
                //singularFaces.insert(face);
            //}

            //++curr;
        //} while (curr != start);
    //}

    //printf("There are %ld singular faces out of %ld faces. That is %.2f%%.\n", singularFaces.size(), arr.number_of_faces(), 100.0 * (float)singularFaces.size() / (float)arr.number_of_faces());





}

void Arrangement::computePointLocationDataStructure()
{
    this->pl = std::make_unique<Point_location>(this->arr);
}
