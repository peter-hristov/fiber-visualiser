#include <cmath>
#include <iostream>
#include <netcdf.h>
#include <string>
#include <vector>

#include "src/GlobalConfig.h"

#define ERR(e)                                                                                                         \
    {                                                                                                                  \
        printf("Error: %s\n", nc_strerror(e));                                                                         \
    }

namespace netCDFWrapper {
std::string
errorForCode(int errorCode)
{
    return (std::string("Error: ") + std::string(nc_strerror(errorCode)));
}

namespace dimension {
int
getId(const int ncId, const std::string dimName)
{
    int retval;
    int dimId;

    if ((retval = nc_inq_dimid(ncId, dimName.c_str(), &dimId))) {
        throw errorForCode(retval);
    }
    return dimId;
}

size_t
getLength(const int ncId, const std::string dimName)
{
    int retval;
    int dimId = getId(ncId, dimName);

    size_t dimLength;
    if ((retval = nc_inq_dimlen(ncId, dimId, &dimLength))) {
        throw errorForCode(retval);
    }

    return dimLength;
}
}

namespace variable {
size_t
readDimsNum(const int ncId, const int varId)
{
    int dims;
    int retval;

    if ((retval = nc_inq_varndims(ncId, varId, &dims))) {
        throw errorForCode(retval);
    }

    return dims;
}

std::string
readName(const int ncId, const int varId)
{
    int retval;

    //
    // Get the variable name
    //
    char* varName = new char[NC_MAX_NAME];
    if ((retval = nc_inq_varname(ncId, varId, varName))) {
        throw errorForCode(retval);
    }

    return std::string(varName);
}

std::string
readAttribute(const int ncId, const int varId, const std::string attrName)
{
    size_t attrLength;
    int retval;

    // Get the length of the variable attribute
    if ((retval = nc_inq_attlen(ncId, varId, attrName.c_str(), &attrLength))) {
        throw errorForCode(retval);
    }

    // Read the contents of the variable attribute
    char* attrNameC = new char[attrLength + 1];
    if ((retval = nc_get_att_text(ncId, varId, attrName.c_str(), attrNameC))) {
        throw errorForCode(retval);
    }

    attrNameC[attrLength] = '\0';

    return std::string(attrNameC);
}

int
getVariableId(const int ncId, const std::string variableName)
{
    int varId;
    int retval;

    if ((retval = nc_inq_varid(ncId, variableName.c_str(), &varId))) {
        throw errorForCode(retval);
    }

    return varId;
}

std::vector<double>
readDoubleDimensionData(const int ncId, const std::string variableName, const size_t min, const size_t max)
{
    int retval;
    int varId = netCDFWrapper::variable::getVariableId(ncId, variableName.c_str());

    size_t start[1] = { min };
    size_t count[1] = { max - min };

    // Raw data read from the file
    double* rawData = new double[count[0]];
    if ((retval = nc_get_vara_double(ncId, varId, start, count, rawData))) {
        throw netCDFWrapper::errorForCode(retval);
    }

    return std::vector<double>(rawData, rawData + count[0]);
}

std::vector<std::vector<std::vector<std::vector<float>>>>
readFloatVariableData(const int ncId, const int varId, const int nDims, tv9k::InputInformation input)
{
    // Used for netcdf error handling
    int retval;

    // Data in cpp format

    const size_t *start, *count;

    if (3 == nDims) {
        start = new (size_t[3]){ input.yMin, input.xMin, input.zMin };
        count = new (size_t[3]){ input.yMax - input.yMin, input.xMax - input.xMin, input.zMax - input.zMin };
    } else if (4 == nDims) {
        start = new (size_t[4]){ input.tMin, input.yMin, input.xMin, input.zMin };
        count = new (size_t[4]){
            input.tMax - input.tMin, input.yMax - input.yMin, input.xMax - input.xMin, input.zMax - input.zMin
        };
    }

    std::vector<std::vector<std::vector<std::vector<float>>>> data(
      count[0],
      std::vector<std::vector<std::vector<float>>>(
        count[2], std::vector<std::vector<float>>(count[1], std::vector<float>(count[3], 0))));

    // Raw data read from the file
    float* rawData = new float[count[0] * count[1] * count[2] * count[3]];
    if ((retval = nc_get_vara_float(ncId, varId, start, count, rawData))) {
        throw netCDFWrapper::errorForCode(retval);
    }

    // Convert raw data to cpp 4D vector format
    for (int i = 0; i < count[2]; i++) {
        for (int j = 0; j < count[1]; j++) {
            for (int k = 0; k < count[3]; k++) {
                for (int t = 0; t < count[0]; t++) {

                    float value =
                      rawData[t * (count[1] * count[2] * count[3]) + j * (count[2] * count[3]) + i * count[3] + k];

                    // @TODO Figure out what to do with NAN values
                    if (std::isnan(value)) {
                        value = 0.0;
                    }

                    data[t][i][j][k] = value;
                }
            }
        }
    }

    delete[] rawData;
    delete[] start;
    delete[] count;

    return data;
}

// std::vector<double> readDimensionValues()
//{

//}

}
};
