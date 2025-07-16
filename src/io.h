#pragma once

#include <string>

#include "./TetMesh.h"

namespace io
{
    TetMesh readData(const std::string&);
    TetMesh readDataTxt(const std::string&);
    TetMesh readDataVtu(const std::string&);
}
