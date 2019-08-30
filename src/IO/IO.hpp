#pragma once

// #include <scai/hmemo/HArray.hpp>
// #include <scai/hmemo/WriteAccess.hpp>
// #include <scai/lama.hpp>
#include "../Common/HostPrint.hpp"
#include <scai/lama/Vector.hpp>
//
// #include <cmath>

namespace KITGPI
{
    //! \brief IO namespace
    namespace IO
    {
        using namespace scai;
        /*! \brief Write lama vector to an external file
 *
 *  Write a lama vector to an external file block.
 \param vector lama vector which will be written to filename
 \param filename Name of file in which vector will be written
 \param fileFormat Output file format 1=mtx 2=lmf 3=frv
 */
        template <typename ValueType>
        void writeVector(lama::Vector<ValueType> const &vector, std::string filename, IndexType fileFormat)
        {
            HOST_PRINT(vector.getDistributionPtr()->getCommunicatorPtr(), "writing " << filename << ", fileFormat = " << fileFormat << "\n")

            switch (fileFormat) {
            case 1:
                vector.writeToFile(filename + ".mtx", lama::FileMode::FORMATTED);
                break;

            case 2:
                // write binary file, IndexType as int, ValueType as float, do it via collective I/O
                vector.writeToFile(filename + ".lmf", lama::FileMode::BINARY, common::ScalarType::FLOAT, common::ScalarType::INT);
                break;
            case 3:
                filename += ".frv"; // write binary file with separate header file, done by master process
                vector.writeToFile(filename, lama::FileMode::BINARY, common::ScalarType::FLOAT, common::ScalarType::INT);
                break;

            default:
                COMMON_THROWEXCEPTION("Unexpected fileFormat option!")
                break;
            }
        }

        /*! \brief Read a Vector from file
 * 
 \param vector lama vector which will be written to filename the size of the input vector must be known before reading the data. The vector will be redistributed to the dist of the input vector.
 \param filename Name of file in which vector will be written
 \param fileFormat Output file format 1=mtx 2=lmf 3=frv
 */
        template <typename ValueType>
        void readVector(scai::lama::Vector<ValueType> &vector, std::string filename, IndexType fileFormat)
        {
            HOST_PRINT(vector.getDistributionPtr()->getCommunicatorPtr(), "reading " << filename << ", fileFormat = " << fileFormat << "\n");

            switch (fileFormat) {
            case 1:
                filename += ".mtx";
                break;
            case 2:
                filename += ".lmf";
                break;
            case 3:
                filename += ".frv";
                break;

            default:
                COMMON_THROWEXCEPTION("Unexpected fileFormat option!")
                break;
            }

            vector.readFromFile(filename, vector.getDistributionPtr());
        }
    }
}
