#include "./ReebSpace.h"
#include "./DisjointSet.h"
#include "./CGALTypedefs.h"
#include "./Timer.h"
#include "./utility/indicators.hpp"

#include <cassert>
#include <cstddef>
#include <map>
#include <queue>
#include <set>
#include <iterator>
#include <unordered_map>

std::pair<std::set<int>, std::set<int>> ReebSpace::getMinusPlusTrianglesIndex(const TetMesh &tetMesh, const Arrangement &arrangement, const Arrangement_2::Halfedge_const_handle currentHalfEdge)
{
    // Step 2. Find the edge in the mesh corresponding to the segment corresponding to the half edge
    const Segment_2 &segment = *arrangement.arr.originating_curves_begin(currentHalfEdge);
    //std::cout << "Half-edge   from: " << currentHalfEdge->source()->point() << " to " << currentHalfEdge->target()->point() << std::endl;
    //std::cout << "Source-edge from: " << segment.source() << " to " << segment.target() << std::endl;

    // These will always be sorted, it's how we created the segments
    const int aIndex = arrangement.arrangementPointsIdices.at(segment.source());
    const int bIndex = arrangement.arrangementPointsIdices.at(segment.target());

    // Sanity check
    assert(aIndex < bIndex);

    const std::array<int, 2> edge = {aIndex, bIndex};

    //printf("The original indices are %d and %d", data->arrangementPointsIdices[segment.source()], data->arrangementPointsIdices[segment.target()]);
    //printf("\n");


    // Check to see if the segment and half edge have the same orientation
    const bool isSegmentLeftToRight = segment.source() < segment.target(); 
    const bool isCurrentHalfEdgeLeftToRight = (currentHalfEdge->direction() == CGAL::ARR_LEFT_TO_RIGHT);

    // The half edge has the same direction as the original edge
    if (isSegmentLeftToRight == isCurrentHalfEdgeLeftToRight)
    {
        return {tetMesh.upperStarTriangles.at(edge), tetMesh.lowerStarTriangles.at(edge)};
    }
    else
    {
        return {tetMesh.lowerStarTriangles.at(edge), tetMesh.upperStarTriangles.at(edge)};
    }
}



void ReebSpace::testTraverseArrangement(const Arrangement &arrangement)
{
    // Find the unbounded face (hold the boundary of the arrangement)
    Face_const_handle outerFace;

    // Iterate over all faces and find the unbounded one
    for (Face_const_iterator fit = arrangement.arr.faces_begin(); fit != arrangement.arr.faces_end(); ++fit) 
    {
        if (fit->is_unbounded()) 
        {
            // If this has already been set, we have two outerFace, this should never happen.
            assert(outerFace == Face_const_handle());
            outerFace = fit;
            break;
        }
    }

    // Sanity check, make sure the outer face is simple
    assert(outerFace->number_of_inner_ccbs() == 1);
    assert(outerFace->number_of_outer_ccbs() == 0);
    assert(outerFace->number_of_holes() == 1);
    assert(outerFace->number_of_isolated_vertices() == 0);

    // Get the the first half edge of the outerFace (could be any edge, this is a matter of convention)
    Halfedge_const_handle outerHalfEdge = *outerFace->holes_begin();

    // Make sure there is only one originating curve, something has gone wrong otherwise (edge overlap)
    assert(std::distance(arrangement.arr.originating_curves_begin(outerHalfEdge), arrangement.arr.originating_curves_end(outerHalfEdge)) == 1);

    // Starting from an outer half edge
    std::queue<Arrangement_2::Halfedge_const_handle> traversalQueue;
    traversalQueue.push(outerHalfEdge);

    // Which faces are visited
    std::set<Arrangement_2::Face_const_handle> visited;
    visited.insert(outerHalfEdge->face());

    while (false == traversalQueue.empty())
    {
        // Pop an half edge out
        Arrangement_2::Halfedge_const_handle currentHalfEdge = traversalQueue.front();
        traversalQueue.pop();
        Arrangement_2::Face_const_handle currentFace = currentHalfEdge->face();

        // Get the twin
        Arrangement_2::Halfedge_const_handle twin = currentHalfEdge->twin();
        Arrangement_2::Face_const_handle twinFace = twin->face();

        // Get ids of the current face and the twin face
        int currentFaceID = arrangement.arrangementFacesIdices.at(currentFace);
        int twinFaceID = arrangement.arrangementFacesIdices.at(twinFace);

        // If we have never visited this face, then we have never visited any of the half edges.
        if (visited.find(twinFace) == visited.end())
        {
            //printf("NEW FACE ------------------------------------------ %d -> %d \n", data->arrangementFacesIdices[currentFace], data->arrangementFacesIdices[twinFace]);
            visited.insert(twinFace);

            // Sanity check we should only have onbounded face, the outside face.
            assert(false == twinFace->is_unbounded());


            Arrangement_2::Ccb_halfedge_const_circulator start = twinFace->outer_ccb();
            Arrangement_2::Ccb_halfedge_const_circulator curr = start;

            do {
                traversalQueue.push(curr);

                // Make sure there is only one originating curve (sanity check)
                //const Segment_2 &segment = *data->arr.originating_curves_begin(curr);
                //std::cout << "Half-edge   from: " << curr->source()->point() << " to " << curr->target()->point() << std::endl;
                //std::cout << "Source-edge from: " << segment.source() << " to " << segment.target() << std::endl;
                //printf("The original indices are %d and %d", data->arrangementPointsIdices[segment.source()], data->arrangementPointsIdices[segment.target()]);
                //printf("\n\n");

                ++curr;
            } while (curr != start);
        }
    }
}





void ReebSpace::computeTwinFacePreimageGraph(const TetMesh &tetMesh, const Arrangement &arrangement, const Arrangement_2::Halfedge_const_handle &currentHalfEdge)
{
    // Get the face
    Arrangement_2::Face_const_handle currentFace = currentHalfEdge->face();

    // Get the twin
    Arrangement_2::Halfedge_const_handle twin = currentHalfEdge->twin();
    Arrangement_2::Face_const_handle twinFace = twin->face();

    // Get ids of the current face and the twin face
    int currentFaceID = arrangement.arrangementFacesIdices.at(currentFace);
    int twinFaceID = arrangement.arrangementFacesIdices.at(twinFace);

    auto [minusTriangles, plusTriangles] = ReebSpace::getMinusPlusTrianglesIndex(tetMesh, arrangement, currentHalfEdge);

    //std::cout << "Computed minus/plus triangles..." << std::endl;
    //printf("The current preimage graph has %d triangles...\n", data->preimageGraphs[currentFaceID].data.size());

    // Set the current preimage graph to be the preimage graph of the parent
    std::set<int> preimageGraph;
    for (const auto &[t, id] : preimageGraphs[currentFaceID].data)
    {
        preimageGraph.insert(t);
    }

    //printf("Minus triangles...\n");

    // Add and remove triangles of the upper/lower link of the crossed edge
    for (const auto &triangle: minusTriangles)
    {
        preimageGraph.erase(triangle);
        //std::cout << triangle << std::endl;
    }

    //printf("Plus triangles...\n");
    for (const auto &triangle: plusTriangles)
    {
        preimageGraph.insert(triangle);
        //std::cout << triangle << std::endl;
    }

    //std::cout << "Computed preimage graph soup..." << std::endl;
    //for (const auto &triangle: preimageGraph)
    //{
        //std::cout << triangle << std::endl;
    //}

    //
    // Step 3. Compute disjointSets[twinFaceID]
    //

    this->preimageGraphs[twinFaceID].initialize(preimageGraph);

    for (const auto &[t1, id1] : this->preimageGraphs[twinFaceID].data)
    {
        for (const auto &t2 : tetMesh.tetIncidentTriangles[t1])
        {
            if (this->preimageGraphs[twinFaceID].data.contains(t2))
            {
                this->preimageGraphs[twinFaceID].union_setsTriangle(t1, t2);
            }
        }
    }

    //std::cout << "Unioned sets..." << std::endl;

    // Finaly make sure everyon points to their root
    this->preimageGraphs[twinFaceID].update();

    // Used when drawing the arrangement
    //data->arrangementFiberComponents[twinFaceID] = data->preimageGraphs[twinFaceID].countConnectedComponents();
}

void ReebSpace::computePreimageGraphs(const TetMesh &tetMesh, const Arrangement &arrangement, const bool discardFiberSeedsSets)
{

    using namespace indicators;
    ProgressBar bar{
        option::BarWidth{50},
            option::Start{"["},
            option::Fill{"■"},
            option::Lead{"■"},
            option::Remainder{"-"},
            option::End{" ]"},
            option::PostfixText{"Computing Reeb space."},
            option::ShowPercentage{true},  // Show percentage on the bar
            option::ShowElapsedTime{true},
            //option::ForegroundColor{Color::cyan},
            //option::FontStyles{std::vector<FontStyle>{FontStyle::bold}}
    };









    // Find the unbounded face (hold the boundary of the arrangement)
    Face_const_handle outerFace;

    // Iterate over all faces and find the unbounded one
    for (Face_const_iterator fit = arrangement.arr.faces_begin(); fit != arrangement.arr.faces_end(); ++fit) 
    {
        if (fit->is_unbounded()) 
        {
            // If this has already been set, we have two outerFace, this should never happen.
            assert(outerFace == Face_const_handle());
            outerFace = fit;
            break;
        }
    }

    // Sanity check, make sure the outer face is simple
    assert(outerFace->number_of_inner_ccbs() == 1);
    assert(outerFace->number_of_outer_ccbs() == 0);
    assert(outerFace->number_of_holes() == 1);
    assert(outerFace->number_of_isolated_vertices() == 0);

    // Loop the edges of the boundary
    //Arrangement_2::Ccb_halfedge_const_circulator inner_ccb = *outerFace->holes_begin();
    //Arrangement_2::Ccb_halfedge_const_circulator he = inner_ccb;

    // Iterate over the halfedges forming the inner CCB
    //do {
        //std::cout << he->source()->point() << " -> ";
        //std::cout << he->target()->point() << std::endl;

        // Get the originating curve of the half edge
        // Maybe need a function for this
        //const Segment_2& segment = *data->arr.originating_curves_begin(he);
        //std::cout << data->arrangementPointsIdices[segment.source()] << ", ";
        //std::cout << data->arrangementPointsIdices[segment.target()] << std::endl;
    //} while (++he != inner_ccb);

    //printf("\n\n");

    // Get the the first half edge of the outerFace (could be any edge, this is a matter of convention)
    //Halfedge_const_handle outerHalfEdge = *outerFace->holes_begin();

    // Make sure there is only one originating curve, something has gone wrong otherwise (edge overlap)
    //assert(std::distance(data->arr.originating_curves_begin(outerHalfEdge), data->arr.originating_curves_end(outerHalfEdge)) == 1);

    // Extract the original curve
    //const Segment_2& segment = *data->arr.originating_curves_begin(outerHalfEdge);

    //printf("The half edge came from this edge:\n");
    //std::cout << segment.source() << " -> " << segment.target() << std::endl;

    //printf("The original indices are %d and %d", data->arrangementPointsIdices[segment.source()], data->arrangementPointsIdices[segment.target()]);
    //printf("\n\n");







    // Starting from an outer half edge
    std::queue<Face_const_handle> traversalQueue;
    // This is the order in which faces are added, also servers as a visited flag
    std::vector<int> order(arrangement.arrangementFacesIdices.size(), -1);


    traversalQueue.push(outerFace);

    // This is the order in which a face has been processed, note this is different than level
    // This is used as an index for faces, so that we don't do double work when checking edges for correspondence
    int orderIndex = 0;
    order[arrangement.arrangementFacesIdices.at(outerFace)] = orderIndex;

    // The disjoint set to track the connected components of the preimage graph
    this->preimageGraphs.resize(arrangement.arrangementFacesIdices.size());

    // The number of connected components for each preimage graph (computed from the disjoint set)
    //data->arrangementFiberComponents.resize(data->arrangementFacesIdices.size(), -1);

    // If we want the fiber seeds, initialize them
    if (false == discardFiberSeedsSets)
    {
        this->fiberSeeds.resize(arrangement.arrangementFacesIdices.size());
    }

    int graphsInMemory = 0;
    float averageAraphsInMemory = 0;
    int barTickThreshold = arrangement.arrangementFacesIdices.size() / 50;

    while (false == traversalQueue.empty())
    {
        // Pop an half edge out
        Face_const_handle currentFace = traversalQueue.front();
        traversalQueue.pop();

        // Get ids of the current face and the twin face
        int currentFaceID = arrangement.arrangementFacesIdices.at(currentFace);

        //std::cout << "Current face " << currentFaceID << std::endl;

        // Sanity check, this should always be true
        if (false == currentFace->is_unbounded()) 
        {
            //assert(false == data->preimageGraphs[currentFaceID].isEmpty());
        }

        //
        // For each neighbouring face
        //

        // Unbounded face (the starting one) has a different way of addressing its neighbours
        Halfedge_const_handle start;
        if (currentFace->is_unbounded()) 
        {  
            start = *currentFace->holes_begin();
        }
        else
        {
            start = *currentFace->outer_ccbs_begin();
        }

        Arrangement_2::Ccb_halfedge_const_circulator curr = start;
        do {

            Face_const_handle twinFace = curr->twin()->face();
            const int twinFaceID = arrangement.arrangementFacesIdices.at(twinFace);

            // If the neighbour has not been visited, we enqueue it and also compute its preimage graph
            if (-1 == order[twinFaceID])
            {
                traversalQueue.push(twinFace);
                order[twinFaceID] = ++orderIndex;

                //std::cout << "Computing for neighbour " << twinFaceID << std::endl;

                // Compute the preimage graph of this unvisited face
                ReebSpace::computeTwinFacePreimageGraph(tetMesh, arrangement, curr);
                graphsInMemory++;

                // Get all unique roots and a representative
                const std::vector<std::pair<int, int>> representativesAndRoots = this->preimageGraphs[twinFaceID].getUniqueRepresentativesAndRoots();

                // Initialize the vertices of H with the connected components of the twin graph
                for (const auto &[representative, root] : representativesAndRoots)
                {
                    this->reebSpace.addElements({twinFaceID, root});
                }

                // If we want fiber computation, cache the seeds for the fiber components
                if (false == discardFiberSeedsSets)
                {
                    this->fiberSeeds[twinFaceID] = representativesAndRoots;
                }
            }

            // Sanity check, all graphs should have either been computed before or now
            //assert(false == data->preimageGraphs[twinFaceID].isEmpty());

            // Compute the correspondence with the neighbours, but only if they are at a higher level, or we are at the same level, currentFaceID < twinFaceID is used to avoid double work, we only need it once
            if (order[currentFaceID] < order[twinFaceID])
            {
                //std::cout << "Determining correspondence with neighbour " << twinFaceID << std::endl;
                ReebSpace::determineCorrespondence(tetMesh, arrangement, curr);
            }

            ++curr;
        } while (curr != start);

        averageAraphsInMemory = averageAraphsInMemory + ((float)graphsInMemory - (float)averageAraphsInMemory) / (float)orderIndex;

        //printf("There are %d active preimage graphs with average %f at index %d/%ld.\n", graphsInMemory, averageAraphsInMemory, orderIndex, data->preimageGraphs.size());

        // If the threshold is zero the ticks bugs out, so we don't do it
        if (barTickThreshold > 0 && order[currentFaceID] % barTickThreshold == 0)
        {
            // Update bar state
            if (false == bar.is_completed()) {  bar.tick(); }
            if (false == bar.is_completed()) {  bar.tick(); }
        }


        // Dispose of the preimage graph we will no longer need it
        this->preimageGraphs[currentFaceID].clear();
        graphsInMemory--;
    }

    bar.set_progress(100); // all done
    printf("\n\nThere is an average of %f / %ld active preimage graphs.\n", averageAraphsInMemory, this->preimageGraphs.size());
    printf("The correspondence graphs has %ld vertices and the Reeb space has %ld sheets.\n\n", this->reebSpace.data.size(), this->reebSpace.getUniqueRepresentativesAndRoots().size());
}


void ReebSpace::determineCorrespondence(const TetMesh &tetMesh, const Arrangement &arrangement, const Arrangement_2::Halfedge_const_handle &halfEdge)
{
    // Get the face IDs of the current face and its twin
    Face_const_handle face = halfEdge->face();
    Face_const_handle twinFace = halfEdge->twin()->face();

    const int faceID = arrangement.arrangementFacesIdices.at(face);
    const int twinFaceID = arrangement.arrangementFacesIdices.at(twinFace);

    // Get the originating edge
    const Segment_2 &segment = *arrangement.arr.originating_curves_begin(halfEdge);
    const std::pair<int, int> originatingEdge = {arrangement.arrangementPointsIdices.at(segment.source()), arrangement.arrangementPointsIdices.at(segment.target())};

    // The triangls that are added/removd from face -> twinFace
    const auto& [minusTriangles, plusTriangles] = ReebSpace::getMinusPlusTrianglesIndex(tetMesh, arrangement, halfEdge);

    // See which of the roots (connected components) are active (contain an active triangle)
    std::set<int> activeRootsFace;
    for (int triangle : minusTriangles)
    {
        activeRootsFace.insert(this->preimageGraphs[faceID].findTriangle(triangle));
    }

    std::set<int> activeRootsTwinFace;
    for (int triangle : plusTriangles)
    {
        activeRootsTwinFace.insert(preimageGraphs[twinFaceID].findTriangle(triangle));

    }

    // Definite edge
    //if (activeRootsFace.size() == 0 || activeRootsTwinFace.size() == 0)
    //{

        //data->jacobiType[originatingEdge] = 0;
    //}
    //// Reeb-regular
    //else if (activeRootsFace.size() == 1 && activeRootsTwinFace.size() == 1)
    //{

        //data->jacobiType[originatingEdge] = 1;
    //}
    //// Indefinite edge
    //else
    //{
        //data->jacobiType[originatingEdge] = 2;

    //}

    // This is a regular edge then connect the two active components, for singular edge we skip this, they are considered new
    if (activeRootsFace.size() == 1 && activeRootsTwinFace.size() == 1)
    {
        // Active roots point to each other
        // Add edge [faceId, activeRootsFace.begin()->second()] -> [twinFaceId, activeRootsTwinFace.begin()->second()]
        //connectedFacesAndRoots
        //std::set<int> faceRoot({faceID, *activeRootsFace.begin()->second(});
        //std::pair<std::set<int>, std::set<int>> reebSpaceConnection({t1, t2});


        // @TODO Ulgy, but hard to get the first element out otherwise
        for (const int &rFace :activeRootsFace)
        {
            for (const int &rTwinFace :activeRootsTwinFace)
            {
                //connectedFacesAndRoots.insert(reebSpaceConnection);
                //data->edgesH.push_back({
                        //{faceID, rFace},
                        //{twinFaceID, rTwinFace}
                        //});


                this->reebSpace.union_setsTriangle({faceID, rFace}, {twinFaceID, rTwinFace});



                //int faceVertexHindex = data->vertexHtoIndex[{faceID, rFace}];
                //int twinFaceVertexHindex = data->vertexHtoIndex[{twinFaceID, rTwinFace}];

                    //data->indexToVertexH.push_back(vertexH);
                    //data->vertexHtoIndex[vertexH] = verexHIndex;


                //data->reebSpace.union_setsTriangle(faceVertexHindex, twinFaceVertexHindex);

            }
        }
    }

    //
    // Link together all other connected components
    //
    for (const auto &[t, id] : this->preimageGraphs[faceID].data)
    {
        // The root of the triangle in the face
        const int triangleRootFace = this->preimageGraphs[faceID].findTriangle(t);

        // We have already deal with the active fiber
        if (activeRootsFace.contains(triangleRootFace)) { continue; }

        // The root of the triangle in the twin face
        const int triangleRootTwinFace = this->preimageGraphs[twinFaceID].findTriangle(t);

        //connectedFacesAndRoots.insert(reebSpaceConnection);
        //data->edgesH.push_back({
                //{faceID, triangleRootFace}, 
                //{twinFaceID, triangleRootTwinFace}
                //});
        this->reebSpace.union_setsTriangle({faceID, triangleRootFace}, {twinFaceID, triangleRootTwinFace});

        //int faceVertexHindex = data->vertexHtoIndex[{faceID, triangleRootFace}];
        //int twinFaceVertexHindex = data->vertexHtoIndex[{twinFaceID, triangleRootTwinFace}];
        //data->reebSpace.union_setsTriangle(faceVertexHindex, twinFaceVertexHindex);
    }
}




void ReebSpace::computeCorrespondenceGraph(const TetMesh &tetMesh, const Arrangement &arrangement)
{

    // For evey face
    for (auto face = arrangement.arr.faces_begin(); face != arrangement.arr.faces_end(); ++face) 
    {
        // Skip the outer face
        if (face->is_unbounded()) { continue; }

        // Walk around the boundary of the face
        Arrangement_2::Ccb_halfedge_const_circulator start = face->outer_ccb();
        Arrangement_2::Ccb_halfedge_const_circulator curr = start;
        do {

            ReebSpace::determineCorrespondence(tetMesh, arrangement, curr);
            ++curr;
        } while (curr != start);

    }
}



void ReebSpace::computeReebSpacePostprocess(const TetMesh &tetMesh, const Arrangement &arrangement)
{

    Timer::start();
    // Path compression to make sure every element of H is poiting to its root
    this->reebSpace.update();
    Timer::stop("Updating RS disjoint set               :");


    // For faster lookups cache the sheetd IDs for each face
    std::vector<std::unordered_set<int>> faceSheets(this->fiberSeeds.size());
    for (int i = 0 ; i < this->fiberSeeds.size() ; i++)
    {
        for (const auto &[triangleId, fiberComponentId] : this->fiberSeeds[i])
        {
            const int sheetId = this->reebSpace.findTriangle({i, fiberComponentId});
            faceSheets[i].insert(sheetId);
        }

    }

    Timer::start();
    //
    // In order to compute the polygon of each sheet, first obtain a halfEdge of the arrangement that is on the boundary of the sheet
    //
    std::unordered_map<int, Arrangement_2::Halfedge_const_handle> sheetSeeds;
    for (auto currentFaceIterator = arrangement.arr.faces_begin(); currentFaceIterator != arrangement.arr.faces_end(); ++currentFaceIterator) 
    {
        Arrangement_2::Face_const_handle currentFace = currentFaceIterator;
        if (currentFace->is_unbounded()) { continue; }

        const int currentFaceID = arrangement.arrangementFacesIdices.at(currentFace);

        //printf("\n\nFace %d has these sheets - ", currentFaceID);

        // Which sheets contain the current face
        std::unordered_set<int> currentFaceSheetIds = faceSheets[currentFaceID];
        //for (const auto &[triangleId, fiberComponentId] : data->fiberSeeds[currentFaceID])
        //{
            //const int sheetId = data->reebSpace.findTriangle({currentFaceID, fiberComponentId});
            //currentFaceSheetIds.push_back(sheetId);
            ////printf("%d ", sheetId);
        //}
        //printf("\n");


        Arrangement_2::Ccb_halfedge_const_circulator start = currentFace->outer_ccb();
        Arrangement_2::Ccb_halfedge_const_circulator curr = start;
        do {
            Face_const_handle twinFace = curr->twin()->face();
            if (twinFace->is_unbounded()) { curr++; continue; }

            const int twinFaceID = arrangement.arrangementFacesIdices.at(twinFace);
            std::unordered_set<int> twinFaceSheetIds = faceSheets[twinFaceID];

            // Which sheets are in the currentFace, but NOT in the twin face
            std::vector<int> diff;
            for (const int &sheetId : currentFaceSheetIds) 
            {
                if (false == twinFaceSheetIds.contains(sheetId))
                {
                    diff.push_back(sheetId);
                }
            }

            for (const int &sheetId : diff)
            {
                //std::cout << "Half-edge   from: " << curr->source()->point() << " to " << curr->target()->point() << " is a seed for sheet " << sheetId << std::endl;

                if (false == sheetSeeds.contains(sheetId))
                {
                    sheetSeeds[sheetId] = curr;
                }
            }

            ++curr;
        } while (curr != start);
    }
    Timer::stop("Computing sheet seeds                  :");





    //
    // At this point we have a half edge seed for every sheet, we use it to traverse the boundary of the sheet to obtain its geometry
    // We rely on the half edge data structure provided by the arrangement
    // if ab is an edge between vertices a and b in A and, then we can get the next edge bc as the one of the neighours of b that has
    // the sheet in its face and it doesn't have the sheet in the twin face
    //


    //printf("Done with sheet seeds\n");


    Timer::start();
    for (const auto &[sheetId, halfEdge]: sheetSeeds)
    {
        // Sanity check, make sure the seed we have picked actually is on the boundary of the sheet
        Face_const_handle faceSeedEdge = halfEdge->face();
        const int faceSeedEdgeId = arrangement.arrangementFacesIdices.at(faceSeedEdge);
        if (false == faceSheets[faceSeedEdgeId].contains(sheetId))
        {
            throw std::runtime_error("The face of the seed edge does is not in the sheet.");
        }

        Face_const_handle faceSeedTwinEdge = halfEdge->twin()->face();
        const int faceSeedTwinEdgeId = arrangement.arrangementFacesIdices.at(faceSeedTwinEdge);
        if (true == faceSheets[faceSeedTwinEdgeId].contains(sheetId))
        {
            throw std::runtime_error("The face of the twin of the seed edge does is no in the sheet.");
        }


        //printf("---------------------------------------------------------------------------------------------       At sheet %d \n", sheetId);

        // Starting location, we don't make it visited on purpose, we want to discover it latex to finish the loop
        Vertex_const_handle startVertex = halfEdge->source();

        // Time to find the next vertex
        Vertex_const_handle currentVertex = startVertex;
        std::unordered_set<Vertex_const_handle> visited;
        do 
        {

            // Add the current vertex to the polygon
            this->sheetPolygon[sheetId].push_back({
                    CGAL::to_double(currentVertex->point().x()), 
                    CGAL::to_double(currentVertex->point().y())
                    });


            //printf("-------------------------- The current vertex (%f, %f) \n ", CGAL::to_double(currentVertex->point().x()), CGAL::to_double(currentVertex->point().y()));



            //
            // Find the next vertex on the boundary of the sheet, it will be adjacent to the current vertex
            //
            int sheetBoundaryNeighbours = 0;
            int allNeighbours = 0;
            visited.insert(currentVertex);

            // This loop should always finish
            const auto begin = currentVertex->incident_halfedges();
            auto circ = begin;
            do {

                allNeighbours++;

                // Too many neighourss
                if (allNeighbours > arrangement.arr.number_of_vertices())
                {
                    break;
                }

                const auto twinHalfEdge = circ->twin();
                const auto nextVertex = twinHalfEdge->target();
                //printf("Looking at neighbour (%f, %f) \n ", CGAL::to_double(twinHalfEdge->target()->point().x()), CGAL::to_double(twinHalfEdge->target()->point().y()));

                Face_const_handle faceA = circ->face();
                const int faceAId = arrangement.arrangementFacesIdices.at(faceA);

                Face_const_handle faceB = twinHalfEdge->face();
                const int faceBId = arrangement.arrangementFacesIdices.at(faceB);

                // If one of the face contains the sheet
                const bool faceHalfEdgeContainsSheet = faceSheets[faceAId].contains(sheetId);
                const bool faceTwinHalfEdgeContainsEdge = faceSheets[faceBId].contains(sheetId);

                if (faceHalfEdgeContainsSheet == true && faceTwinHalfEdgeContainsEdge == false)
                {
                    currentVertex = nextVertex;
                    sheetBoundaryNeighbours++;
                }

                ++circ;
            } while (circ != begin);

            // A the polygon of a sheet sheet cannot have more vertices than the number of the vertices in the arrangement, something went wrong.
            if (this->sheetPolygon[sheetId].size() > arrangement.arr.number_of_vertices())
            {
                printf("The boundary of sheet %d is degenerate, it has more vertices %ld than the arrangement.\n", sheetId, this->sheetPolygon[sheetId].size());
                this->incompleteSheets.insert(sheetId);
                break;
            }
            
            if (sheetBoundaryNeighbours !=1)
            {
                printf("The boundary of sheet %d is degenerate, a vertex has %d boundary neighours, should only be one.\n", sheetId, sheetBoundaryNeighbours);
                this->incompleteSheets.insert(sheetId);
                break;
            }

            // I'm not sure if this will ever get triggered, since if we loop back to someone else than the start 
            // that should mean 
            if (currentVertex != startVertex && visited.contains(currentVertex))
            {
                printf("The boundary of sheet %d is degenerate, a vertex loops back to an already visited one.\n", sheetId);
                this->incompleteSheets.insert(sheetId);
                break;

            }

            if (allNeighbours > arrangement.arr.number_of_vertices())
            {
                printf("The boundary of sheet %d is degenerate, the inner loop kept going for too long.\n", sheetId);
                this->incompleteSheets.insert(sheetId);
                break;
            }

        } while (currentVertex != startVertex);
    }
    Timer::stop("Computing sheet boundary polygons      :");

    if (this->incompleteSheets.size() == 0)
    {
        std::cout << "\nThe boundaries of all sheets are okay! No degeneracy.\n" << std::endl;
    }


    // Compute the areas of sheets


    Timer::start();
    for (const auto &[sheetId, polygon] : this->sheetPolygon)
    {
        float area = 0.0;

        // If the sheet is incomplete sum up the faces that make it
        if (this->incompleteSheets.contains(sheetId))
        {
            // Loop through all faces to see which ones are in the sheet
            for (auto f = arrangement.arr.faces_begin(); f != arrangement.arr.faces_end(); ++f) 
            {
                const int currentFaceID = arrangement.arrangementFacesIdices.at(f);

                // For each fiber component in the face, see if one of those is in our sheet
                for (const auto &[triangleId, fiberComponentId] : this->fiberSeeds[currentFaceID])
                {
                    const int componentSheetId = this->reebSpace.findTriangle({currentFaceID, fiberComponentId});

                    // Now we can add the polygon
                    if (componentSheetId == sheetId)
                    {
                        CartesianPolygon_2 poly;

                        typename Arrangement_2::Ccb_halfedge_const_circulator circ = f->outer_ccb();
                        typename Arrangement_2::Ccb_halfedge_const_circulator curr = circ;
                        do {
                            typename Arrangement_2::Halfedge_const_handle e = curr;

                            // Get point from CGAL (and convert to double )
                            const float u = CGAL::to_double(e->source()->point().x());
                            const float v = CGAL::to_double(e->source()->point().y());

                            poly.push_back({u, v});
                        } while (++curr != circ);

                        area += poly.area();
                    }
                }
            }
            printf("Incomplete sheet %d has area %f\n", sheetId, area);
        }
        else
        {
            area = abs(polygon.area());
        }

        this->sheetArea[sheetId] = area;
        
        //printf("The polygon of sheet %d has size %ld and area %f\n", sheetId, polygon.size(), area);
        //printf("The area of sheet %d is %f \n", sheetId, area);
    }
    Timer::stop("Computing sheet areas                  :");
    //Timer::stop("Computing sheet boundary polygons      :");

    Timer::start();
    // Transfer the map entries into a vector of pairs so that we can sort
    std::vector<std::pair<int, double>> sheetAreaSortVector(this->sheetArea.begin(), this->sheetArea.end());

    //printf("The sort vector has size %ld\n", sheetAreaSortVector.size());

    // Sort by area
    std::sort(sheetAreaSortVector.begin(), sheetAreaSortVector.end(), [](const std::pair<int, double>& a, const std::pair<int, double>& b) {
            return a.second > b.second; 
            });

    // Assign a sorted order to each sheet
    for (int i = 0 ; i < sheetAreaSortVector.size() ; i++)
    {
        const auto &[sheetId, area] = sheetAreaSortVector[i];
        this->sheetToColour[sheetId] = i;
    }
    Timer::stop("Sorting sheets and labeling them       :");

    std::cout << "\nHere are the top 20 sheets sorted by range area.\n";
    // Print to debug, at least the first new
    //for (int i = 0 ; i < sheetAreaSortVector.size() ; i++)
    const int maxSheets = std::min((size_t)30, sheetAreaSortVector.size());
    for (int i = 0 ; i < maxSheets ; i++)
    {
        std::cout << i << " -- sheet " << sheetAreaSortVector[i].first << " has area " << sheetAreaSortVector[i].second << std::endl;
    }
}

//void ReebSpace::countIntersectionsTypes(Data *data)
//{
    //int regularRegularIntersections = 0;
    //int singularRegularIntersections = 0;
    //int singularSingularIntersections = 0;
    //int degenerateIntersections = 0;
    //int maxNeihbours = 0;

    //// Now loop over all vertices
    //for (auto currentVertex = data->arrangement.arr.vertices_begin(); currentVertex != data->arrangement.arr.vertices_end(); ++currentVertex)
    //{
        //const Point_2& p = currentVertex->point();
        ////std::cout << "Vertex at: " << p << std::endl;

        //// Skip the vertices, we just want to count the inter intersections.
        //if (data->arrangement.arrangementPointsIdices.contains(p))
        //{
            //continue;
        //}


        //// This loop should always finish
        //const auto begin = currentVertex->incident_halfedges();
        //auto circ = begin;

        //int regularNeighbours = 0;
        //int singularNeighbours = 0;

        //do {

            //const Segment_2 &segment = *data->arrangement.arr.originating_curves_begin(circ);

            //// These will always be sorted, it's how we created the segments
            //const int aIndex = data->arrangement.arrangementPointsIdices[segment.source()];
            //const int bIndex = data->arrangement.arrangementPointsIdices[segment.target()];

            //if (data->jacobiType[{aIndex, bIndex}] == 1)
            //{
                //regularNeighbours++;
            //}
            //else
            //{
                //singularNeighbours++;
            //}

            //++circ;
        //} while (circ != begin);

        //if (regularNeighbours + singularNeighbours > 4)
        //{
            //degenerateIntersections++;
        //}

        //maxNeihbours = std::max(maxNeihbours, regularNeighbours + singularNeighbours);



        //if (regularNeighbours == 0)
        //{
            //singularSingularIntersections++;
        //}
        //else if (singularNeighbours == 0)
        //{
            //regularRegularIntersections++;
        //}
        //else
        //{
            //singularRegularIntersections++;
        //}
    //}


    //const int totalIntersctionPoints = regularRegularIntersections + singularRegularIntersections + singularSingularIntersections;

    //printf("Here is a summary of the intersections types:\n");
    //printf("regular   - regular  intersections: %d the ratio is %.2f%%.\n", regularRegularIntersections, 100.0 * (float)regularRegularIntersections / (float)totalIntersctionPoints);
    //printf("singular  - regular  intersections: %d the ratio is %.2f%%.\n", singularRegularIntersections, 100.0 * (float)singularRegularIntersections / (float)totalIntersctionPoints);
    //printf("singular  - sigular  intersections: %d the ratio is %.2f%%.\n", singularSingularIntersections, 100.0 * (float)singularSingularIntersections / (float)totalIntersctionPoints);
    //printf("The number of degenerate intersections is %d with maximum being %d.\n", degenerateIntersections, maxNeihbours);

//}
