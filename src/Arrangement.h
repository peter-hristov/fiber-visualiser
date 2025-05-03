#pragma once

#include <memory>

#include "./Data.h"
#include "./CustomArrangementTypes.h"

namespace CustomArrangement
{

    void computeArrangementCustom(Data *data);

    //std::optional<std::tuple<Point, DataType, DataType>> computeIntersection(ConstSegmentHandler &s1, ConstSegmentHandler &s2);
    std::optional<std::tuple<CustomArrangement::Point, CustomArrangement::DataType, CustomArrangement::DataType>> computeIntersection(Segment &s1, Segment &s2, std::vector<Point> &points);

}
