#include "FreeSurfaceEM.hpp"

/*! \brief Initialitation of the free surface
 *
 \param dist Distribution of wavefields
 \param derivatives Derivative class
 \param modelCoordinates Coordinate class, which eg. maps 3D coordinates to 1D model indices
 \param DT Temporal Sampling
 */
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::FreeSurfaceEM<ValueType>::init(scai::dmemo::DistributionPtr /*dist*/, Derivatives::Derivatives<ValueType> & /*derivatives*/, Acquisition::Coordinates<ValueType> const & /*modelCoordinates*/, ValueType /*DT*/){}

template class KITGPI::ForwardSolver::BoundaryCondition::FreeSurfaceEM<float>;
template class KITGPI::ForwardSolver::BoundaryCondition::FreeSurfaceEM<double>;
