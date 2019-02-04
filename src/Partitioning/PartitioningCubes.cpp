#include "PartitioningCubes.hpp"
using namespace scai;

/*! \brief Getter method for the disttribution pointer *
 *
 */
template <typename ValueType>
dmemo::DistributionPtr KITGPI::Partitioning::PartitioningCubes<ValueType>::getDist() const
{
    SCAI_ASSERT(dist_cubes != nullptr, "Distribution ist not set ");
    return (dist_cubes);
}

/*! \brief Constructor based on the configuration class and the communicator
 *
 \param config Configuration class, which is used to derive all requiered parameters
 \param comm Communicator
 */
template <typename ValueType>
KITGPI::Partitioning::PartitioningCubes<ValueType>::PartitioningCubes(KITGPI::Configuration::Configuration const &config, scai::dmemo::CommunicatorPtr comm)
{
    dist_cubes = calculate(config.get<IndexType>("ProcNX"), config.get<IndexType>("ProcNY"), config.get<IndexType>("ProcNZ"), config.get<IndexType>("NX"), config.get<IndexType>("NY"), config.get<IndexType>("NZ"), comm);
}

//! \brief Calculation of cube partition
/*!
 * This class can create a partition of the wavefield that consists of cubes.
 *
 * procNX*procNY*procNZ have to be equal to comm->getSize();
 *
 \param procNX Number of cores in X-direction
 \param procNY Number of cores in Y-direction
 \param procNZ Number of cores in Z-direction
 \param NX Number of grid points in X-direction
 \param NY Number of grid points in X-direction
 \param NZ Number of grid points in X-direction
 \param comm Communicator for the distribution
 */
template <typename ValueType>
dmemo::DistributionPtr KITGPI::Partitioning::PartitioningCubes<ValueType>::calculate(IndexType procNX, IndexType procNY, IndexType procNZ, IndexType NX, IndexType NY, IndexType NZ, scai::dmemo::CommunicatorPtr comm)
{

    IndexType rank = comm->getRank();
    IndexType numRanks = comm->getSize();

    /* Check some things */
    SCAI_ASSERT(procNX * procNY * procNZ == numRanks, " Number of cores differ between config and actual setting  ")
    SCAI_ASSERT(NX % procNX == 0, " NX % procNX != 0  ");
    SCAI_ASSERT(NY % procNY == 0, " NY % procNY != 0  ");
    SCAI_ASSERT(NZ % procNZ == 0, " NZ % procNZ != 0  ");

    /* Calculate coordinates of the CPU */
    IndexType posNX = rank % procNX;
    IndexType posNZ = rank / (procNX * procNY);
    IndexType posNY = ((rank - (procNX * procNY) * posNZ) / procNX);

    /* Calculate range of grid points */
    IndexType range_x_lower = (NX / procNX) * posNX;
    IndexType range_y_lower = (NY / procNY) * posNY;
    IndexType range_z_lower = (NZ / procNZ) * posNZ;

    IndexType range_x_upper = (NX / procNX) * (posNX + 1);
    IndexType range_y_upper = (NY / procNY) * (posNY + 1);
    IndexType range_z_upper = (NZ / procNZ) * (posNZ + 1);

    IndexType numGlobalGridPoints = NX * NY * NZ;
    IndexType numLocalGridPoints = (numGlobalGridPoints) / numRanks;

    Acquisition::Coordinates coord(NX,NY,NZ);

    /* Determine local indices */
    hmemo::HArray<IndexType> localIndices;
    localIndices.resize(numLocalGridPoints);
    hmemo::WriteAccess<IndexType> write_localIndices(localIndices);
    IndexType i = 0;
    IndexType indice;
    for (IndexType x = 0; x < NX; x++) {
        for (IndexType y = 0; y < NY; y++) {
            for (IndexType z = 0; z < NZ; z++) {
                if (x >= range_x_lower && x < range_x_upper && y >= range_y_lower && y < range_y_upper && z >= range_z_lower && z < range_z_upper) {
                    indice = coord.coordinate2index(x, y, z);
                    write_localIndices[i] = indice;
                    i++;
                }
            }
        }
    }
    SCAI_ASSERT_ERROR(i == numLocalGridPoints, " i != numLocalGridPoints   " << i << numLocalGridPoints);
    write_localIndices.release();

    /* create GeneralDistribution */
    dmemo::DistributionPtr dist_cubus(new dmemo::GeneralDistribution(numGlobalGridPoints, localIndices, true, comm));

    return (dist_cubus);
}

template class KITGPI::Partitioning::PartitioningCubes<double>;
template class KITGPI::Partitioning::PartitioningCubes<float>;
