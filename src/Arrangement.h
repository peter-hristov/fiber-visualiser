#pragma once

#include "./Data.h"
#include "./CustomArrangementTypes.h"

namespace CustomArrangement
{

    void computeArrangementCustom(Data *data);

    template <typename coordType>
    std::optional<Point<coordType>> computeIntersection(const Segment<coordType> &s1, const Segment<coordType> &s2);
}
