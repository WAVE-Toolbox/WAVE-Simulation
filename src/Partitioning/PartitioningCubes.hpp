#pragma once

#include "../Acquisition/Coordinates.hpp"
#include "../Configuration/Configuration.hpp"
#include "Partitioning.hpp"

namespace KITGPI
{

    namespace Partitioning
    {

        //! \brief Creating a partition of cubes
        /*!
         * This class can create a partition of the wavefield that consists of cubes.
         *
         * procNX*procNY*procNZ have to be equal to comm->getSize();
         *
         */
        template <typename ValueType>
        class PartitioningCubes : public Partitioning<ValueType>
        {

          public:
            //! Default constructor
            PartitioningCubes() = delete;

            explicit PartitioningCubes(KITGPI::Configuration::Configuration const &config, scai::dmemo::CommunicatorPtr comm);

            //! Default destructor
            ~PartitioningCubes(){};

            scai::dmemo::DistributionPtr getDist() const;

          private:
            scai::dmemo::DistributionPtr calculate(scai::IndexType procNX, scai::IndexType procNY, scai::IndexType procNZ, scai::IndexType NX, scai::IndexType NY, scai::IndexType NZ, scai::dmemo::CommunicatorPtr comm);

            scai::dmemo::DistributionPtr dist_cubes; //!< Distribution
        };
    }
}
