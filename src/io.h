#pragma once

#include <string>

#include "./Data.h"

namespace io
{
    Data* readData(const std::string&);
    Data* readDataTxt(const std::string&);
    Data* readDataVtu(const std::string&);
}
