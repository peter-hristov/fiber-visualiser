#pragma once

#include <vector>
#include <memory>
#include <gmpxx.h> 

namespace CustomArrangement
{
    template <typename coordType>
        struct CustomPoint {
            coordType x, y;
            int index;

            CustomPoint(coordType _x, coordType _y, int _index)
            {
                this->x = _x;
                this->y = _y;
                this->index = _index;
            }

            bool operator<(const CustomPoint& other) const {
                return (x < other.x) || (x == other.x && y < other.y);
            }

            bool operator==(const CustomPoint& other) const {
                return x == other.x && y == other.y;
            }
        };

    template <typename coordType>
        struct CustomSegment {
            // We assume that *a < *b
            std::shared_ptr<CustomPoint<coordType>> a, b;
            int index;

            CustomSegment(std::shared_ptr<CustomPoint<coordType>> _a, std::shared_ptr<CustomPoint<coordType>> _b, int index)
            {
                this->a = _a;
                this->b = _b;
                this->index = index;
            }

        };

    // Typedef (or using alias) for shared_ptr<Point<coordType>>
    using DataType = mpq_class;

    using Point = CustomPoint<DataType>;
    using Segment = CustomSegment<DataType>;

    using PointHandler = std::shared_ptr<Point>;
    using SegmentHandler = std::shared_ptr<Segment>;

    using ConstPointHandler = const std::shared_ptr<Point>;
    using ConstSegmentHandler = const std::shared_ptr<const Segment>;

}

