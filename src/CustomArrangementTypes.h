#pragma once

#include <vector>
#include <memory>

namespace CustomArrangement
{
    template <typename coordType>
        struct Point {
            coordType x, y;
            int index;

            Point(coordType _x, coordType _y, int _index)
            {
                this->x = _x;
                this->y = _y;
                this->index = _index;
            }

            bool operator<(const Point& other) const {
                return (x < other.x) || (x == other.x && y < other.y);
            }

            bool operator==(const Point& other) const {
                return x == other.x && y == other.y;
            }
        };

    template <typename coordType>
        struct Segment {
            // We assume that *a < *b
            std::shared_ptr<Point<coordType>> a, b;
            int index;

            Segment(std::shared_ptr<Point<coordType>> _a, std::shared_ptr<Point<coordType>> _b, int index)
            {
                this->a = _a;
                this->b = _b;
                this->index = index;
            }

        };
}

