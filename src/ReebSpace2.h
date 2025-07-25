#pragma once

#include "./TetMesh.h"
#include "./Arrangement.h"
#include "./CGALTypedefs.h"

namespace ReebSpace2
{
    void compute(const TetMesh &tetMesh, Arrangement &regularArrangement, Arrangement &singularArrangement);
};
