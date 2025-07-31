// Printing out exact numbers 
Point_2 ip = std::get<Point_2>(*result);
std::cout << std::endl << CGAL::exact(ip);

// Iterate over all halfedges
//
for (auto he = singularArrangement.arr.halfedges_begin(); he != singularArrangement.arr.halfedges_end(); ++he)
{
    // If you want to not iterate over the same edge twice
    if (he < he->twin())
    {
    }
}

// Iterate over the halfedges of every vertex
for (auto v = regularArrangement.arr.vertices_begin(); v != regularArrangement.arr.vertices_end(); ++v)
{
    if (v->is_isolated()) 
    {
        std::cout << "The vertex (" << v->point() << ") is isolated\n";
        return;
    }

    Arrangement_2::Halfedge_around_vertex_const_circulator first, curr;
    first = curr = v->incident_halfedges();

    std::cout << "The neighbors of the vertex (" << v->point() << ") are:";
    do 
    {
        std::cout << " (" << curr->source()->point() << ")";
    } while (++curr != first);

    std::cout << std::endl;
}


// Segment_2 to edge

const Segment_2 &segment = mySegment.seg;



// Iterate over the halfedges of every faces
// The order is always CCW
for (auto fit = singularArrangement.arr.faces_begin(); fit != singularArrangement.arr.faces_end(); ++fit)
{
    if (fit->is_unbounded())
    {
        std::cout << " (unbounded)" << std::endl;
        continue;
    }

    std::cout << " (bounded):" << std::endl;

    auto circ = fit->outer_ccb();
    auto start = circ;
    do
    {
        const Point_2& p = circ->target()->point();
        std::cout << "  Vertex: (" << p << ")" << std::endl;
        ++circ;
    } while (circ != start);
}




// Get originating segment
const Segment_2 &segment = *arrangement.arr.originating_curves_begin(currentHalfEdge);
//std::cout << "Half-edge   from: " << currentHalfEdge->source()->point() << " to " << currentHalfEdge->target()->point() << std::endl;
//std::cout << "Source-edge from: " << segment.source() << " to " << segment.target() << std::endl;

const int aIndex = arrangement.arrangementPointIndices.at(segment.source());
const int bIndex = arrangement.arrangementPointIndices.at(segment.target());

// Sanity check
assert(aIndex < bIndex);

const std::array<int, 2> edge = {aIndex, bIndex};
