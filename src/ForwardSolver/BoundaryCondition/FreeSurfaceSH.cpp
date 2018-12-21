#include "FreeSurfaceSH.hpp"
using namespace scai;

//! Default destructor
template <typename ValueType>
KITGPI::ForwardSolver::BoundaryCondition::FreeSurfaceSH<ValueType>::~FreeSurfaceSH(){};

/*! \brief Initialitation of the free surface
 *
 \param dist Distribution of wavefields
 \param derivatives Derivative class
 \param NX Number of grid points in X-Direction
 \param NY Number of grid points in Y-Direction (Depth)
 \param NZ Number of grid points in Z-Direction
 \param modelCoordinates Coordinate class, which eg. maps 3D coordinates to 1D model indices
 \param DT Temporal Sampling
 \param DH Distance between grid points
 */
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::FreeSurfaceSH<ValueType>::init(scai::dmemo::DistributionPtr /*dist*/, Derivatives::Derivatives<ValueType> &/*derivatives*/, IndexType /*NX*/, IndexType /*NY*/, IndexType /*NZ*/, Acquisition::Coordinates const &/*modelCoordinates*/,ValueType /*DT*/, ValueType /*DH*/)
{
COMMON_THROWEXCEPTION(" Image method is not implemented for Love-Waves ");
}

template class KITGPI::ForwardSolver::BoundaryCondition::FreeSurfaceSH<float>;
template class KITGPI::ForwardSolver::BoundaryCondition::FreeSurfaceSH<double>;
