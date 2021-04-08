#pragma once

#include "../Common/HostPrint.hpp"
#include <scai/lama/Vector.hpp>

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
            auto comm = vector.getDistributionPtr()->getCommunicatorPtr();

            switch (fileFormat) {
            case 1:
                filename += ".mtx";
                HOST_PRINT(comm, "", "writing " << filename << " (MatrixMarket)\n")
                vector.writeToFile(filename, lama::FileMode::FORMATTED);
                break;

            case 2:
                // write binary file, IndexType as int, ValueType as float, do it via collective I/O
                filename += ".lmf";
                HOST_PRINT(comm, "", "writing " << filename << " (LAMA proprietary binary, float)\n")
                vector.writeToFile(filename, lama::FileMode::BINARY, common::ScalarType::FLOAT, common::ScalarType::INT);
                break;
            case 3:
                filename += ".frv"; // write binary file with separate header file, done by master process
                HOST_PRINT(comm, "", "writing " << filename << " (binary + separate header)\n")
                vector.writeToFile(filename, lama::FileMode::BINARY);
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

            IndexType numelements_exp=vector.size();
            
            HOST_PRINT(vector.getDistributionPtr()->getCommunicatorPtr(), "", "reading " << filename << "\n");
            vector.readFromFile(filename, vector.getDistributionPtr());
            
            IndexType numelements_read=vector.size();
            
            SCAI_ASSERT(numelements_exp == numelements_read, "Read " << numelements_read << " elements from file: " << filename << ", expected " << numelements_exp << " elements!");
            
        }

        /*! \brief Write lama Matrix to an external file
        *
        *  Write a lama matrix to an external file block.
        \param matrix lama matrix which will be written to filename
        \param filename Name of file in which matrix will be written
        \param fileFormat Output file format 1=mtx 2=lmf 3=frv
        */
        template <typename ValueType>
        void writeMatrix(lama::Matrix<ValueType> const &matrix, std::string filename, IndexType fileFormat)
        {

            switch (fileFormat) {
            case 1:
                filename += ".mtx";
                matrix.writeToFile(filename, lama::FileMode::FORMATTED);
                break;
            case 2:
                // write binary file, IndexType as int, ValueType as float, do it via collective I/O
                filename += ".lmf";
                matrix.writeToFile(filename, lama::FileMode::BINARY, common::ScalarType::FLOAT, common::ScalarType::INT);
                break;
            case 3:
                filename += ".frv"; // write binary file with separate header file, done by master process
                matrix.writeToFile(filename, lama::FileMode::BINARY);
                break;

            default:
                COMMON_THROWEXCEPTION("Unexpected fileFormat option!")
                break;
            }

            HOST_PRINT(matrix.getRowDistributionPtr()->getCommunicatorPtr(), "", "writing " << filename << "\n")

        }

        /*! \brief Read a Matrix from file
        * 
        \param matrix lama matrix which will be written to filename the size of the input matrix must be known before reading the data. The matrix will be redistributed to the dist of the input matrix.
        \param filename Name of file in which matrix will be written
        \param fileFormat Output file format 1=mtx 2=lmf 3=frv
        */
        template <typename ValueType>
        void readMatrix(scai::lama::Matrix<ValueType> &matrix, std::string filename, IndexType fileFormat)
        {

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
            HOST_PRINT(matrix.getRowDistributionPtr()->getCommunicatorPtr(), "", "reading " << filename << "\n");

            IndexType numrows_exp=matrix.getNumRows();
            IndexType numcols_exp=matrix.getNumColumns();
            
            matrix.readFromFile(filename, matrix.getRowDistributionPtr());
            matrix.redistribute(matrix.getRowDistributionPtr(), matrix.getColDistributionPtr());
            
            IndexType numrows_read=matrix.getNumRows();
            IndexType numcols_read=matrix.getNumColumns();           
            
            SCAI_ASSERT(numrows_exp == numrows_read, "Read " << numrows_read << " rows from file: " << filename << " expected " << numrows_exp << " rows!");
            SCAI_ASSERT(numcols_exp == numcols_read, "Read " << numcols_read << " columns from file: " << filename << " expected " << numcols_exp << " columns!");

        }

        /*! \brief Read single row of a Matrix from file
        * 
        \param filename Name of file in which matrix will be written
        \param fileFormat Output file format 1=mtx 2=lmf 3=frv
        */
        template <typename ValueType>
        hmemo::HArray<ValueType> readMatrix(std::string filename, IndexType rowNumber, IndexType fileFormat)
        {

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
            scai::lama::DenseStorage<ValueType> matrix;
            matrix.readFromFile(filename, rowNumber, 1);

            hmemo::HArray<ValueType> localsignal = matrix.getValues();
            return (localsignal);
        }
    }
}
