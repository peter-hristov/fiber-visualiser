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
for (auto vit = singularArrangement.arr.vertices_begin(); vit != singularArrangement.arr.vertices_end(); ++vit)
{
    auto circ = vit->incident_halfedges();

    if (circ != nullptr) {
        auto start = circ;

        std::cout << "Vertex at " << vit->point() << " incident halfedges:\n";
        do {
            Halfedge_const_handle heHandle = circ;

            std::cout << "  Halfedge from " << circ->source()->point()
                << " to " << heHandle->target()->point() << "\n";
            ++circ;
        } while (circ != start);
    } else {
        std::cout << "Vertex at " << vit->point() << " has no incident halfedges.\n";
    }
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
