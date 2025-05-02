#pragma once

#include "./Data.h"
#include "./CustomArrangementTypes.h"
#include <memory>

namespace CustomArrangement
{

    void computeArrangementCustom(Data *data);

    template <typename coordType>
    std::optional<Point<coordType>> computeIntersection(const std::shared_ptr<const Segment<coordType>> &s1, const std::shared_ptr<const Segment<coordType>> &s2);

}
