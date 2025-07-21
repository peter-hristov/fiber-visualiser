#pragma once

#include <string>

#include "./TetMesh.h"
#include "./Arrangement.h"
#include "./ReebSpace.h"

#include "./FiberPoint.h"

namespace io
{
    TetMesh readData(const std::string&);
    TetMesh readDataTxt(const std::string&);
    TetMesh readDataVtu(const std::string&);

    void saveSheets(const TetMesh &tetMesh, const Arrangement &arrangement, const ReebSpace &reebSpace, const std::string &outputSheetPolygonsFilename);

    void saveFibers(const std::string&, const std::vector<FiberPoint>&);
}
