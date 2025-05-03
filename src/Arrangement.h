#pragma once

#include <memory>

#include "./Data.h"
#include "./CustomArrangementTypes.h"

namespace CustomArrangement
{

    void computeArrangementCustom(Data *data);

    std::optional<std::tuple<PointHandler, DataType, DataType>> computeIntersection(ConstSegmentHandler &s1, ConstSegmentHandler &s2);

}
