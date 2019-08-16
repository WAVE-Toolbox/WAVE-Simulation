#include "FreeSurfaceAcoustic.hpp"
using namespace scai;

//! Default destructor
template <typename ValueType>
KITGPI::ForwardSolver::BoundaryCondition::FreeSurfaceAcoustic<ValueType>::~FreeSurfaceAcoustic(){};

/*! \brief Initialitation of the free surface
 *
 \param dist Distribution of wavefields
 \param derivatives Derivative class
 \param modelCoordinates Coordinate class, which eg. maps 3D coordinates to 1D model indices
 \param DT Temporal Sampling
 */
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::FreeSurfaceAcoustic<ValueType>::init(scai::dmemo::DistributionPtr dist, Derivatives::Derivatives<ValueType> &derivatives, Acquisition::Coordinates<ValueType> const &modelCoordinates, ValueType DT)
{
    dmemo::CommunicatorPtr comm = dist->getCommunicatorPtr();

    HOST_PRINT(comm, "", "Initialization of the free surface...\n");

    active = true;

    HOST_PRINT(comm, "", "Finished initializing of the free surface\n\n");
}

template class KITGPI::ForwardSolver::BoundaryCondition::FreeSurfaceAcoustic<float>;
template class KITGPI::ForwardSolver::BoundaryCondition::FreeSurfaceAcoustic<double>;
