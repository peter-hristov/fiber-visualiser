#include "Data.h"
#include "./utility/Geometry.h"
#include "./utility/MarchingCubes.h"
#include "./utility/NetCDFWRapper.h"
#include <cassert>
#include <cstdio>
#include <netcdf.h>
#include <omp.h>
#include <qpoint.h>
#include <sys/stat.h>

#include "./utility/Histogram.h"
#include "src/utility/MergeTree.h"

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/linestring.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/multi_point.hpp>
#include <boost/geometry/geometries/multi_polygon.hpp>

using namespace std;


void
Data::readNcData(tv9k::InputInformation input)
{
    // This will be the netCDF ID for the file and data variable.
    int ncId;

    // Used to netcdf error handling
    int retval;

    // Open the file.
    if ((retval = nc_open(input.filename.c_str(), NC_NOWRITE, &ncId))) {
        throw netCDFWrapper::errorForCode(retval);
    }

    // Red number of dimensions, variables, attributes and unlimited dims.
    int ndims, nvars, ngatts, unlimdimid;
    if ((retval = nc_inq(ncId, &ndims, &nvars, &ngatts, &unlimdimid))) {
        throw netCDFWrapper::errorForCode(retval);
    }

    // int ndims, nvars, ngatts, unlimdimid;
    // if ((retval = nc_inq(ncId, &ndims, &nvars, &ngatts, &unlimdimid))) {
    // throw netCDFWrapper::errorForCode(retval);
    //}

    // Read Dimensions
    originalXdim = netCDFWrapper::dimension::getLength(ncId, input.xName);
    originalYdim = netCDFWrapper::dimension::getLength(ncId, input.yName);
    originalZdim = netCDFWrapper::dimension::getLength(ncId, input.zName);

    // Only read time if it's there and set the time values
    if (4 == ndims) {
        originalTdim = netCDFWrapper::dimension::getLength(ncId, input.tName);
    } else {
        originalTdim = 1;
    }

    // Crop the dimensions based on input arguments (does not include
    // downsampling, that's done later)
    std::tie(this->xdim, this->ydim, this->zdim, this->tdim) =
      input.cropDimensions(this->originalXdim, this->originalYdim, this->originalZdim, this->originalTdim);

    cout << "Original Dimensions are : " << originalXdim << " , " << originalYdim << " , " << originalZdim << " , "
         << originalTdim << endl;

    cout << "New Dimensions are : " << xdim << " , " << ydim << " , " << zdim << " , " << tdim << endl;

    cout << "Reading input data... " << std::flush;
    // Read the values associated with the time dimensions
    try {
        this->tVals = netCDFWrapper::variable::readDoubleDimensionData(ncId, input.tName, input.tMin, input.tMax);
    } catch (string exception) {
        cerr << exception << ". Default to index values." << endl;
        this->tVals = vector<double>(input.tMax - input.tMin);
        std::iota(std::begin(tVals), std::end(tVals), 0);
    }

    // Read all scalar fields
    for (int i = 0; i < nvars; i++) {

        // Does it have full dimensions?
        if (netCDFWrapper::variable::readDimsNum(ncId, i) == ndims) {
            // Initialize field
            tv9k::utility::ScalarField field;
            field.name = netCDFWrapper::variable::readName(ncId, i);

            // Read in data
            field.values = netCDFWrapper::variable::readFloatVariableData(ncId, i, ndims, input);

            // field.clipDomain(xMin, xMax, yMin, yMax, zMin, zMax, tMin, tMax);
            field.downsampleDomain(input.downsampleX, input.downsampleY, input.downsampleZ);
            field.computeMinMaxValues();

            try {
                field.longName = netCDFWrapper::variable::readAttribute(ncId, i, "longname");
            } catch (string exception) {
                try {
                    field.longName = netCDFWrapper::variable::readAttribute(ncId, i, "long_name");
                } catch (string exception) {
                    field.longName = "N/A";
                }
            }
            try {
                field.units = netCDFWrapper::variable::readAttribute(ncId, i, "units");
            } catch (string exception) {
                field.units = "N/A";
            }

            this->scalarFields[field.name] = field;
        }
    }
    cout << "done" << std::endl;

    // Read the join trees or compute them
    std::string joinTreeFilename = "";
    int fileReadRetureValue = 0;
    auto isoFieldName = input.attributeNames[0];

    // Read tree from data file
    if (true == this->cacheJoinTree) {
        std::string infix = "";
        if (!this->precomputeMergeTrees) {
            infix = std::string(".").append(isoFieldName);
        }

        joinTreeFilename = input.filename.append(infix).append(".tree");

        fileReadRetureValue = this->readTreeData(joinTreeFilename, this->scalarFields);

        if (0 == fileReadRetureValue)
        {
            // Not much to do if it's all good
        }
        else if (1 == fileReadRetureValue)
        {
            cerr << "Cout now open file " << joinTreeFilename << " for reading.";
        }
        else if (2 == fileReadRetureValue)
        {
            cerr << "File " << joinTreeFilename << " does not have the same dimensions as the input data set ";
        }
        else
        {
            assert(false);
        }
    }

    // See which fields don't have join tree computed and compute them
    if (true == this->precomputeMergeTrees) {
        for (auto& field : this->scalarFields) {
            if (field.second.joinTrees.size() == 0)
            {
                field.second.computeJoinTrees();
            }
            else
            {
                // Just in case, make sure we have the correct number of trees
                assert(field.second.joinTrees.size() == tdim);
            }
        }
    }
    // If we're not precomputing sure we at least have the isofield join tree
    else
    {
        // If they haven't been computed so far either they haven't been read in or we haven't precomputed everything
        if (this->scalarFields[isoFieldName].joinTrees.size() == 0){
            this->scalarFields[isoFieldName].computeJoinTrees();
        }
    }

    // If we've failed to read from the file (missing or incorrect dimensions) and we're caching we should overwrite it with a valid file 
    if (0 != fileReadRetureValue && true == this->cacheJoinTree) {
        this->saveTreeData(joinTreeFilename, this->scalarFields);
    }


    // Set the current iso, u and v fields
    this->changeIsoField(isoFieldName, 0);
    this->changeUField(input.attributeNames[1]);
    this->changeVField(input.attributeNames[2]);

    // Make sure we have enough meshes for every timestep
    this->isosurfaceMeshes.resize(this->tdim);
    this->fibersurfaceMeshes.resize(this->tdim);
    this->combinedMeshes.resize(this->tdim);

    for (const auto u : this->scalarFields)
    {

        this->mousePoints.insert({u.first, std::map<std::string, QVector<QPointF>>()});

        for (const auto v : this->scalarFields)
        {
            this->mousePoints[u.first].insert({v.first, QVector<QPointF>()});
            //this->mousePoints[std::pair(u.second,)]

        }
    }

    this->currentSignedDistanceField.resize(this->tdim);

    // Set default currentID and scene centre
    for (int i = 0; i < this->tdim; i++) {
        this->selectedSurfaceType.push_back(SurfaceType::none);
        this->selectedID.push_back({ -1 });
        this->selectedObjectMinMax.push_back({ { 0, 0, 0 }, { 0, 0, 0 } });
    }

    valsF = std::vector<std::vector<std::vector<GLfloat>>>(
            xdim, std::vector<std::vector<GLfloat>>(ydim, std::vector<GLfloat>(zdim, 0)));
}

int
Data::readTreeData(const std::string filename, std::map<std::string, tv9k::utility::ScalarField>& fields)
{
    FILE* file;
    file = fopen(filename.c_str(), "r");
    if (!file) {
        return 1;
        //throw "Could not open binary file for reading.\n";
    }

    int treeTdim, treeXdim, treeYdim, treeZdim;
    fread(&treeTdim, sizeof(treeTdim), 1, file);
    fread(&treeXdim, sizeof(treeXdim), 1, file);
    fread(&treeYdim, sizeof(treeYdim), 1, file);
    fread(&treeZdim, sizeof(treeZdim), 1, file);

    int joinTreeFieldsNumber = 0;
    fread(&joinTreeFieldsNumber, sizeof(joinTreeFieldsNumber), 1, file);

    // Case some useful things
    size_t dataSize = treeXdim * treeYdim * treeZdim;
    float* rawFloatData = new float[dataSize];
    int* rawIntData = new int[dataSize];

    for (int fieldFileIndex = 0; fieldFileIndex < joinTreeFieldsNumber; fieldFileIndex++) {
        //
        // Read the field name
        //
        int fieldNameSize;
        fread(&fieldNameSize, sizeof(fieldNameSize), 1, file);

        char* fieldNameRaw = new char[fieldNameSize + 1];
        fread(fieldNameRaw, sizeof(fieldNameRaw[0]), fieldNameSize, file);

        // Terminate and convertd to std string;
        fieldNameRaw[fieldNameSize] = '\0';
        std::string fieldName(fieldNameRaw);

        // Make sure we have that field
        assert(this->scalarFields.find(fieldName) != this->scalarFields.end());

        // See if the dimensions in the file match our current dimensions
        if (
                treeTdim != fields[fieldName].tDim() ||
                treeXdim != fields[fieldName].xDim() ||
                treeYdim != fields[fieldName].yDim() ||
                treeZdim != fields[fieldName].zDim()
           )
        {
            return 2;
        }

        //assert(treeTdim == fields[fieldName].tDim());
        //assert(treeXdim == fields[fieldName].xDim());
        //assert(treeYdim == fields[fieldName].yDim());
        //assert(treeZdim == fields[fieldName].zDim());

        fields[fieldName].joinTrees.clear();
        for (int t = 0; t < treeTdim; t++) {
            MergeTree currentTree;

            // Set tree dimensions
            currentTree.xdim = treeXdim;
            currentTree.ydim = treeYdim;
            currentTree.zdim = treeZdim;

            //
            // Read data
            //

            // Read raw data
            fread(rawFloatData, sizeof(rawFloatData[0]), dataSize, file);

            // Initialise data
            currentTree.data = std::vector<std::vector<std::vector<float>>>(
                    treeXdim, std::vector<std::vector<float>>(treeYdim, std::vector<float>(treeZdim, 0)));

            // Get data from raw array in tree
            for (int i = 0; i < treeXdim; i++) {
                for (int j = 0; j < treeYdim; j++) {
                    for (int k = 0; k < treeZdim; k++) {

                        // TODO Remove the isovalue mult
                        currentTree.data[i][j][k] = rawFloatData[j * (xdim * zdim) + i * zdim + k];
                        // assert(currentTree.data[i][j][k] == field->joinTrees[t].data[i][j][k]);
                    }
                }
            }

            //
            // Read Merge Tree
            //

            // Read raw data
            fread(rawIntData, sizeof(rawIntData[0]), dataSize, file);

            // Initialize merge tree array
            currentTree.mergeTree = std::vector<int>(dataSize, 0);

            // Transfer values
            for (int i = 0; i < currentTree.mergeTree.size(); i++) {
                currentTree.mergeTree[i] = rawIntData[i];
                // assert(currentTree.mergeTree[i] == field->joinTrees[t].mergeTree[i]);
            }

            //
            // Read Persistence Pairs
            //

            // Read the number of pairs
            int persistencePairsCount;
            fread(&persistencePairsCount, sizeof(persistencePairsCount), 1, file);

            // Read raw data
            std::tuple<int, int, int>* pairs = new std::tuple<int, int, int>[persistencePairsCount];
            fread(pairs, sizeof(pairs[0]), persistencePairsCount, file);

            // Initialize merge tree array
            currentTree.persistencePairs = std::vector<std::tuple<int, int, int>>(persistencePairsCount);

            // Transfer values
            for (int i = 0; i < persistencePairsCount; i++) {
                currentTree.persistencePairs[i] = pairs[i];
                // assert(currentTree.persistencePairs[i] == field->joinTrees[t].persistencePairs[i]);
            }
            delete[] pairs;

            //
            // Read Vis
            //

            // Read raw data
            fread(rawIntData, sizeof(rawIntData[0]), dataSize, file);

            // Initialize vis array
            currentTree.vis = std::vector<std::vector<std::vector<int>>>(
                    treeXdim, std::vector<std::vector<int>>(treeYdim, std::vector<int>(treeZdim, 0)));

            // Transfer data from the raw array to the tree
            for (int i = 0; i < treeXdim; i++) {
                for (int j = 0; j < treeYdim; j++) {
                    for (int k = 0; k < treeZdim; k++) {
                        currentTree.vis[i][j][k] = rawIntData[j * (xdim * zdim) + i * zdim + k];
                        // assert(currentTree.vis[i][j][k] == field->joinTrees[t].vis[i][j][k]);
                    }
                }
            }

            //
            // Read Simplified Vis
            //

            // Read raw data
            fread(rawIntData, sizeof(rawIntData[0]), dataSize, file);

            // Initialize vis array
            currentTree.simplifiedVisited = std::vector<std::vector<std::vector<int>>>(
                    treeXdim, std::vector<std::vector<int>>(treeYdim, std::vector<int>(treeZdim, 0)));

            // Transfer data from the raw array to the tree
            for (int i = 0; i < treeXdim; i++) {
                for (int j = 0; j < treeYdim; j++) {
                    for (int k = 0; k < treeZdim; k++) {
                        currentTree.simplifiedVisited[i][j][k] = rawIntData[j * (xdim * zdim) + i * zdim + k];
                        // assert(currentTree.simplifiedVisited[i][j][k] ==
                        // field->joinTrees[t].simplifiedVisited[i][j][k]);
                    }
                }
            }

            currentTree.isovalueMult = -1;
            fields[fieldName].joinTrees.push_back(currentTree);
        }
    }
    // Clean up
    delete[] rawFloatData;
    delete[] rawIntData;
    fclose(file);

    return 0;
}

    void
Data::saveTreeData(const std::string filename, const std::map<std::string, tv9k::utility::ScalarField> fields)
{
    FILE* file;
    file = fopen(filename.c_str(), "w");
    if (!file) {
        throw "Could not open binary file for writing.\n";
    }

    // Write dimensions
    fwrite(&tdim, sizeof(tdim), 1, file);
    fwrite(&xdim, sizeof(xdim), 1, file);
    fwrite(&ydim, sizeof(ydim), 1, file);
    fwrite(&zdim, sizeof(zdim), 1, file);

    //
    // How fields have computed their join trees?
    //
    int joinTreeFields = 0;
    for (const auto& field : fields) {
        if (field.second.joinTrees.size() != 0) {
            joinTreeFields++;
        }
    }
    fwrite(&joinTreeFields, sizeof(joinTreeFields), 1, file);

    // Cache someuseful things
    size_t dataSize = xdim * ydim * zdim;
    float* rawFloatData = new float[dataSize];
    int* rawIntData = new int[dataSize];

    for (const auto& field : fields) {
        // Skp the ones for which we have not computed the scalar fields
        if (field.second.joinTrees.size() == 0) {
            continue;
        }

        //
        // Write the name of the field
        //
        const int fieldNameSize = field.first.size();
        fwrite(&fieldNameSize, sizeof(fieldNameSize), 1, file);

        const char* fieldName = field.first.c_str();
        fwrite(fieldName, sizeof(fieldName[0]), fieldNameSize, file);

        //
        // Write the join tree
        //
        for (const auto& currentTree : field.second.joinTrees) {
            for (int i = 0; i < xdim; i++) {
                for (int j = 0; j < ydim; j++) {
                    for (int k = 0; k < zdim; k++) {

                        // TODO Remove the isovalue mult
                        rawFloatData[j * (xdim * zdim) + i * zdim + k] = currentTree.data[i][j][k];
                    }
                }
            }
            fwrite(rawFloatData, sizeof(rawFloatData[0]), dataSize, file);

            // Write Merge Tree
            for (int i = 0; i < currentTree.mergeTree.size(); i++) {
                rawIntData[i] = currentTree.mergeTree[i];
            }
            fwrite(rawIntData, sizeof(rawIntData[0]), dataSize, file);

            //
            // Write Persistence Pairs
            //
            int persistencePairsCount = currentTree.persistencePairs.size();
            std::tuple<int, int, int>* pairs = new std::tuple<int, int, int>[persistencePairsCount];

            for (int i = 0; i < persistencePairsCount; i++) {
                pairs[i] = currentTree.persistencePairs[i];
            }
            fwrite(&persistencePairsCount, sizeof(persistencePairsCount), 1, file);
            fwrite(pairs, sizeof(pairs[0]), persistencePairsCount, file);
            delete[] pairs;

            //
            // Write VIS
            //
            for (int i = 0; i < xdim; i++) {
                for (int j = 0; j < ydim; j++) {
                    for (int k = 0; k < zdim; k++) {
                        rawIntData[j * (xdim * zdim) + i * zdim + k] = currentTree.vis[i][j][k];
                    }
                }
            }
            fwrite(rawIntData, sizeof(rawIntData[0]), dataSize, file);

            //
            // Write simplified VIS
            //
            for (int i = 0; i < xdim; i++) {
                for (int j = 0; j < ydim; j++) {
                    for (int k = 0; k < zdim; k++) {
                        rawIntData[j * (xdim * zdim) + i * zdim + k] = currentTree.simplifiedVisited[i][j][k];
                    }
                }
            }
            fwrite(rawIntData, sizeof(rawIntData[0]), dataSize, file);
        }
    }
    // Clean up
    delete[] rawFloatData;
    delete[] rawIntData;
    fclose(file);
}

    void
Data::changeIsoField(const std::string fieldName, const int mergeTreeType)
{
    assert(this->scalarFields.find(fieldName) != this->scalarFields.end());
    this->isoField = &this->scalarFields[fieldName];
}

    void
Data::changeUField(std::string fieldName)
{
    assert(this->scalarFields.find(fieldName) != this->scalarFields.end());
    this->uField = &this->scalarFields[fieldName];
}

    void
Data::changeVField(std::string fieldName)
{
    assert(this->scalarFields.find(fieldName) != this->scalarFields.end());
    this->vField = &this->scalarFields[fieldName];
}

    void
Data::changeTimeStep(int timestep)
{
    assert(timestep >= 0);
    assert(timestep < this->isoField->tDim());

    this->currentTimestep = timestep;
}


    void
Data::computeCombinedMeshes(const std::vector<std::vector<std::vector<GLfloat>>> signedDistanceField, const GLfloat isovalue)
{
    this->combinedMeshes.clear();
    this->combinedMeshes.resize(this->tdim);

#pragma omp parallel
    {
#pragma omp for
        for (int t = 0; t < this->tdim; t++) {

            // Don't compute if there's no isosurface/fibersurface in this timestep
            if (this->fibersurfaceMeshes[t].triangles.size() == 0 || this->isosurfaceMeshes[t].triangles.size() == 0)
            {
                continue;
            }

            int ID = omp_get_thread_num();
            // cout << "HERE IS THE ID OF THE THREAD " << ID << endl;
            // cout << "------------------------------------ Computing FS mesh for field " << this->isoField->name << "
            // at timestep "<< t << " ." << endl;

            std::vector<std::vector<std::vector<GLfloat>>> combinedSignedDistanceField(
                    xdim, std::vector<std::vector<GLfloat>>(ydim, std::vector<GLfloat>(zdim)));

            // cout << "------------------------------------ Computing FS Distance Fields ";
            // Utility::startTimer();

            for (int i = 0; i < xdim; i++) {
                for (int j = 0; j < ydim; j++) {
                    for (int k = 0; k < zdim; k++) {

                        float iDistance = isovalue - this->isoField->values[t][i][j][k];
                        float uvDistance = this->currentSignedDistanceField[t][i][j][k];

                        //cout << "The UV Distance is " << uvDistance << " the isovalue is " << isovalue << " the iDistance is " << iDistance << endl;

                        float combinedDistance = sqrt(uvDistance * uvDistance + iDistance * iDistance);

                        if (this->isoField->values[t][i][j][k] >= isovalue && uvDistance <= 0)
                        {
                            //printf("x = %f, y = %f, xy-Dist = %f, isovalue = %f, i = %f, iDist = %f, combined = %f, signed = %f\n", rescaledX, rescaledY, uvDistance, isovalue, this->isoField->values[t][i][j][k], iDistance, combinedDistance, signedDistanceField[i][j][k]);
                            combinedSignedDistanceField[i][j][k] = -1 * combinedDistance;
                            //combinedSignedDistanceField[i][j][k] = -1;
                        }
                        else
                        {
                            //combinedSignedDistanceField[i][j][k] = 1;
                            combinedSignedDistanceField[i][j][k] = combinedDistance;
                        }



                        // signedDistanceField[i][j][k] *= -1;
                    }
                }
            }

            // printf(" in - "
            //"%3ld.%06ld seconds.\n",
            // Utility::endTimer().first,
            // Utility::endTimer().second);

            // cout << "------------------------------------ Computing FS Volumes ";
            // Utility::startTimer();
            {
                // Compute the connected components
                std::tie(this->combinedMeshes[t].visited, std::ignore, std::ignore) =
                    Utility::computeVolumes(0, combinedSignedDistanceField, this->xdim, this->ydim, this->zdim);
            }
            // printf(" in - "
            //"%3ld.%06ld seconds.\n",
            // Utility::endTimer().first,
            // Utility::endTimer().second);

            // cout << "------------------------------------ Computing FS Meshes ";
            // Utility::startTimer();
            // Compute the mesh with IDs
            tv9k::utility::MarchingCubes::computeTriangles(0,
                    combinedSignedDistanceField,
                    this->uField->values[t],
                    this->vField->values[t],
                    2,
                    this->combinedMeshes[t],
                    1,
                    this);
            // printf("in - %3ld.%06ld seconds with %ld triangles.\n",
            // Utility::endTimer().first,
            // Utility::endTimer().second,
            // this->fibersurfaceMeshes[t].triangles.size());
            // printf("\n");
        }
    }
}

    void
Data::computeFiberMeshes(const GLfloat resolution, const QVector<QPointF> polyPoints)
{
    typedef boost::geometry::model::d2::point_xy<double> point_type;
    typedef boost::geometry::model::polygon<point_type> polygon_type;

    this->fibersurfaceMeshes.clear();
    this->fibersurfaceMeshes.resize(this->tdim);

    // Convert poly points to uv space
    QVector<QPointF> uvPolyPoints;

    polygon_type uvBoostPolyPointsOuter; 
    polygon_type uvBoostPolyPointsInnerOuter; 

    uvBoostPolyPointsInnerOuter.inners().resize (1);

    for (const auto &p : polyPoints)
    {
        const auto uRatio = p.x() / resolution;
        const auto uCoordinate = (1 - uRatio) * this->uField->min + uRatio * this->uField->max;

        const auto vRatio = p.y() / resolution;
        const auto vCoordinate = (1 - vRatio) * this->vField->min + vRatio * this->vField->max;

        uvPolyPoints.push_back({uCoordinate, vCoordinate});

        boost::geometry::append(uvBoostPolyPointsOuter.outer(), point_type{uCoordinate, vCoordinate});

        boost::geometry::append(uvBoostPolyPointsInnerOuter.outer(), point_type{uCoordinate, vCoordinate});
        boost::geometry::append(uvBoostPolyPointsInnerOuter.inners()[0], point_type{uCoordinate, vCoordinate});
    }

    // Ad the first point at the end
    const auto uRatio = polyPoints[0].x() / resolution;
    const auto uCoordinate = (1 - uRatio) * this->uField->min + uRatio * this->uField->max;

    const auto vRatio = polyPoints[0].y() / resolution;
    const auto vCoordinate = (1 - vRatio) * this->vField->min + vRatio * this->vField->max;

    boost::geometry::append(uvBoostPolyPointsOuter.outer(), point_type{uCoordinate, vCoordinate});

    boost::geometry::append(uvBoostPolyPointsInnerOuter.outer(), point_type{uCoordinate, vCoordinate});
    boost::geometry::append(uvBoostPolyPointsInnerOuter.inners()[0], point_type{uCoordinate, vCoordinate});

#pragma omp parallel
    {
#pragma omp for
        for (int t = 0; t < this->tdim; t++) {
            int ID = omp_get_thread_num();
            // cout << "HERE IS THE ID OF THE THREAD " << ID << endl;
            // cout << "------------------------------------ Computing FS mesh for field " << this->isoField->name << "
            // at timestep "<< t << " ." << endl;

            // Compute the scalar signed distance field
            this->currentSignedDistanceField[t] = std::vector<std::vector<std::vector<GLfloat>>>(
                    xdim, std::vector<std::vector<GLfloat>>(ydim, std::vector<GLfloat>(zdim)));

            // cout << "------------------------------------ Computing FS Distance Fields ";
            // Utility::startTimer();

            for (int i = 0; i < xdim; i++) {
                for (int j = 0; j < ydim; j++) {
                    for (int k = 0; k < zdim; k++) {

                        // Project point to the plane
                        float uCoordinate = this->uField->values[t][i][j][k];
                        float vCoordinate = this->vField->values[t][i][j][k];

                        //point_type uvBoostPoint(uCoordinate, vCoordinate);

                        //const auto isInside = boost::geometry::within(uvBoostPoint, uvBoostPolyPointsOuter);
                        //const auto boostDistance = boost::geometry::distance(uvBoostPoint, uvBoostPolyPointsInnerOuter);

                        //if (true == isInside)
                        //{
                            //this->currentSignedDistanceField[t][i][j][k] = -1.0 * boostDistance;
                        //}
                        //else
                        //{
                            //this->currentSignedDistanceField[t][i][j][k] = boostDistance;
                        //}

                        // Get distance
                        std::tie(this->currentSignedDistanceField[t][i][j][k], std::ignore) =
                            tv9k::geometry::getDistancePointPolygon(uvPolyPoints, QPointF(uCoordinate, vCoordinate));

                        //if (this->currentSignedDistanceField[i][j][k] > 0)
                        //{
                            //printf ("Our signed distance is %f, boot distance is %f and it's inside %d.\n", this->currentSignedDistanceField[i][j][k], boostDistance, isInside);
                        //}

                         //this->currentSignedDistanceField[t][i][j][k] *= -1;
                    }
                }
            }

            // printf(" in - "
            //"%3ld.%06ld seconds.\n",
            // Utility::endTimer().first,
            // Utility::endTimer().second);

            // cout << "------------------------------------ Computing FS Volumes ";
            // Utility::startTimer();
            {
                // Compute the connected components
                std::tie(this->fibersurfaceMeshes[t].visited, std::ignore, std::ignore) =
                    Utility::computeVolumes(0, this->currentSignedDistanceField[t], this->xdim, this->ydim, this->zdim);
            }
            // printf(" in - "
            //"%3ld.%06ld seconds.\n",
            // Utility::endTimer().first,
            // Utility::endTimer().second);

            // cout << "------------------------------------ Computing FS Meshes ";
            // Utility::startTimer();
            // Compute the mesh with IDs
            tv9k::utility::MarchingCubes::computeTriangles(0,
                    this->currentSignedDistanceField[t],
                    this->uField->values[t],
                    this->vField->values[t],
                    1,
                    this->fibersurfaceMeshes[t],
                    -1,
                    this);
            // printf("in - %3ld.%06ld seconds with %ld triangles.\n",
            // Utility::endTimer().first,
            // Utility::endTimer().second,
            // this->fibersurfaceMeshes[t].triangles.size());
            // printf("\n");
        }
    }
}

    void
Data::computeMeshes(const GLfloat isovalue)
{
    this->isosurfaceMeshes.clear();
    this->isosurfaceMeshes.resize(this->tdim);

#pragma omp parallel
    {
#pragma omp for
        for (int i = 0; i < this->tdim; i++) {
            int ID = omp_get_thread_num();
            // cout << "HERE IS THE ID OF THE THREAD " << ID << endl;

            // cout << "------------------------------------ Computing Isosurface mesh for field " <<
            // this->isoField->name << " at timestep "<< i << " with isovalue " << isovalue << "." << endl;

            // Compute Volumes
            // cout << "------------------------------------ Computing Isosurface Volumes ";
            // Utility::startTimer();
            {
                // Get the object ID mask for the 3D domain
                this->isosurfaceMeshes[i].visited = this->isoField->joinTrees[i].computeVisitedForIsovalue(isovalue);
            }
            // printf("in - "
            //"%3ld.%06ld seconds.\n",
            // Utility::endTimer().first,
            // Utility::endTimer().second);

            // Run marching cubes
            // cout << "------------------------------------ Computing Isosurface Meshes ";
            // Utility::startTimer();

            tv9k::utility::MarchingCubes::computeTriangles(isovalue,
                    this->isoField->values[i],
                    this->uField->values[i],
                    this->vField->values[i],
                    0,
                    this->isosurfaceMeshes[i],
                    -1,
                    this);

            // printf("in - %3ld.%06ld seconds with %ld triangles.\n",
            // Utility::endTimer().first,
            // Utility::endTimer().second,
            // this->isosurfaceMeshes[i].triangles.size());
            // printf("\n");
        }
    }
}
