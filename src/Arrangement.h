#pragma once

#include <memory>

#include "./Data.h"
#include "./CustomArrangementTypes.h"

namespace CustomArrangement
{

    void computeArrangementCustom(Data *data);

    //std::optional<std::tuple<Point, DataType, DataType>> computeIntersection(ConstSegmentHandler &s1, ConstSegmentHandler &s2);
    std::optional<std::tuple<CustomArrangement::Point, CustomArrangement::DataType, CustomArrangement::DataType>> computeIntersection(const Segment &s1, const Segment &s2, const std::vector<Point> &segmentPoints);

}
