
#pragma once

#include <scai/lama.hpp>

#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/dmemo/Communicator.hpp>
#include <scai/lama/DenseVector.hpp>
#include <scai/lama/io/PartitionIO.hpp>

#include <iostream>

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

            void readFromDistributedFiles(scai::lama::Vector &vec, std::string const &filename, scai::dmemo::DistributionPtr dist);
            void readFromOneFile(scai::lama::Vector &vec, std::string const &filename, scai::dmemo::DistributionPtr dist);
            void writeToDistributedFiles(scai::lama::Vector const &vec, std::string const &filename);
        };
    }
}
