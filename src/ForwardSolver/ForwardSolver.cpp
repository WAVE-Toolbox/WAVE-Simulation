#include "ForwardSolver.hpp"

/*! \brief memory estimation wrapper for the boundary conditions
*
*
\param config Configuration
\param dist Distribution of the wave fields
\param modelCoordinates Coordinate class, which eg. maps 3D coordinates to 1D model indices
*/
template <typename ValueType>
ValueType KITGPI::ForwardSolver::ForwardSolver<ValueType>::estimateBoundaryMemory(Configuration::Configuration const &config, scai::dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> const &modelCoordinates, BoundaryCondition::ABS<ValueType> &DampingBoundary, BoundaryCondition::CPML<ValueType> &ConvPML)
{
    if (config.get<IndexType>("DampingBoundary") == 1) {
        return (DampingBoundary.estimateMemory(config.get<scai::IndexType>("BoundaryWidth"), config.get<scai::IndexType>("FreeSurface"), dist, modelCoordinates));
    } else if (config.get<IndexType>("DampingBoundary") == 2) {
        return (ConvPML.estimateMemory(config.get<scai::IndexType>("BoundaryWidth"), config.get<scai::IndexType>("FreeSurface"), dist, modelCoordinates));
    }
    return 0;
}

/*! \brief Initialization of the boundary conditions
 *
 *
 \param config Configuration
 \param modelCoordinates Coordinate class, which eg. maps 3D coordinates to 1D model indices
 \param derivatives Derivatives matrices
 \param dist Distribution of the wave fields
 \param ctx Context
 \param FreeSurface Free surface object
 \param DampingBoundary Absorbing Boundary object
 \param ConvPML PML object
 */
template <typename ValueType>
void KITGPI::ForwardSolver::ForwardSolver<ValueType>::prepareBoundaries(Configuration::Configuration const &config, Acquisition::Coordinates<ValueType> const &modelCoordinates, Derivatives::Derivatives<ValueType> &derivatives, scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, KITGPI::ForwardSolver::BoundaryCondition::FreeSurface<ValueType> &FreeSurface, BoundaryCondition::ABS<ValueType> &DampingBoundary, BoundaryCondition::CPML<ValueType> &ConvPML)
{

    useFreeSurface = config.get<scai::IndexType>("FreeSurface");

    /* Prepare Free Surface */
    if (useFreeSurface == 1) {
        FreeSurface.init(dist, derivatives, modelCoordinates, config.get<ValueType>("DT"));
    }

    /* Prepare Damping Boundary */
    if (config.get<IndexType>("DampingBoundary") == 1) {
        useDampingBoundary = true;
        DampingBoundary.init(dist, ctx, modelCoordinates, config.get<scai::IndexType>("BoundaryWidth"), config.get<ValueType>("DampingCoeff"), useFreeSurface);
    }

    if (config.get<IndexType>("DampingBoundary") == 2) {
        useConvPML = true;
        ConvPML.init(dist, ctx, modelCoordinates, config.get<ValueType>("DT"), config.get<scai::IndexType>("BoundaryWidth"), config.get<ValueType>("NPower"), config.get<ValueType>("CenterFrequencyCPML"), config.get<ValueType>("VMaxCPML"), useFreeSurface);
    }
}

template class KITGPI::ForwardSolver::ForwardSolver<double>;
template class KITGPI::ForwardSolver::ForwardSolver<float>;
