#include "FreeSurfaceSH.hpp"
using namespace scai;

/*! \brief Initialization of the free surface
 *
 \param dist Distribution of wavefields
 \param derivatives Derivative class
 \param modelCoordinates Coordinate class, which eg. maps 3D coordinates to 1D model indices
 \param DT Temporal Sampling
 */
template <typename ValueType>
void KITGPI::ForwardSolver::BoundaryCondition::FreeSurfaceSH<ValueType>::init(scai::dmemo::DistributionPtr /*dist*/, Derivatives::Derivatives<ValueType> & /*derivatives*/, Acquisition::Coordinates<ValueType> const & /*modelCoordinates*/, ValueType /*DT*/){}

template class KITGPI::ForwardSolver::BoundaryCondition::FreeSurfaceSH<float>;
template class KITGPI::ForwardSolver::BoundaryCondition::FreeSurfaceSH<double>;
