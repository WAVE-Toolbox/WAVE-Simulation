
#pragma once

#include <scai/lama.hpp>

#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/dmemo/Communicator.hpp>
#include <scai/lama/DenseVector.hpp>
#include <scai/lama/io/PartitionIO.hpp>

#include <iostream>

using namespace scai;
using namespace KITGPI;

namespace KITGPI
{
    namespace PartitionedInOut
    {

        //! Handling parallel input/output
        /*!
         * This class accounts for the parallel input and output.
         * It provides the readin of a single file in block distribution or several files in block distribution and can write to a fileblock 
         * To obtain a fileblock from a single file the  vectorRepartition.exe or matrixRepartition.exe can be executed by 
         *          scai_lama/build/lama/examples/io/vectorRepartition.exe {filename.mtx} 1 {filename.%r.mtx} {NProcessors}
         * To rearrange a fileblock in a single file execute:
         *          scai_lama/build/lama/examples/io/vectorRepartition.exe {filename.%r.mtx} {NProcessors} {filename.mtx} 1
         */
        template <typename ValueType>
        class PartitionedInOut
        {

          public:
            //! Default constructor
            PartitionedInOut(){};

            //! Default destructor
            ~PartitionedInOut(){};

            void readFromDistributedFiles(lama::Vector &vec, std::string const &filename, dmemo::DistributionPtr dist);
            void readFromOneFile(lama::Vector &vec, std::string const &filename, dmemo::DistributionPtr dist);
            void writeToDistributedFiles(lama::Vector const &vec, std::string const &filename);
        };
    }
}

/*! \brief read the distributed file and redistribut the entries
 *
 \param vec Vector which is returned
 \param filename Name of the file Block (without ".%r")
 \param dist Distribution of the vector
 */
template <typename ValueType>
void KITGPI::PartitionedInOut::PartitionedInOut<ValueType>::readFromDistributedFiles(lama::Vector &vec, std::string const &filename, dmemo::DistributionPtr dist)
{
    std::string fileNameBlockIn = filename.substr(0, filename.size() - 4) + ".%r.mtx";
    vec.readFromFile(fileNameBlockIn);
    vec.redistribute(dist);
}

/*! \brief read one file with several processors and redistribut the entries
 *
 \param vec Vector which is returned
 \param filename Name of the file
 \param dist Distribution of the vector
 */
template <typename ValueType>
void KITGPI::PartitionedInOut::PartitionedInOut<ValueType>::readFromOneFile(lama::Vector &vec, std::string const &filename, dmemo::DistributionPtr dist)
{
    vec.readFromFile(filename, "BLOCK");
    vec.redistribute(dist);
}

/*! \brief redistribute vector to block distribution and write to file block
 *
 \param vec Vector which is returned
 \param filename Name of the file
 */
template <typename ValueType>
void KITGPI::PartitionedInOut::PartitionedInOut<ValueType>::writeToDistributedFiles(lama::Vector const &vec, std::string const &filename)
{
    common::unique_ptr<lama::Vector> updatePtr(vec.newVector());
    lama::Vector &tempVector = *updatePtr;
    tempVector = vec;
    std::string distFileName = filename.substr(0, filename.size() - 4) + ".dist.mtx";
    std::string fileNameBlockOut = filename.substr(0, filename.size() - 4) + ".%r.mtx";
    lama::PartitionIO::write(tempVector.getDistribution(), distFileName);

    dmemo::DistributionPtr rowDist(new dmemo::BlockDistribution(tempVector.size(), dmemo::Communicator::getCommunicatorPtr()));
    tempVector.redistribute(rowDist);
    tempVector.writeToFile(fileNameBlockOut);
}
