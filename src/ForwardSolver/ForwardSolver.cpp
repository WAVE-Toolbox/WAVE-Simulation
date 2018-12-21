#include "ForwardSolver.hpp"

/*! \brief Initialitation of the boundary conditions
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
void KITGPI::ForwardSolver::ForwardSolver<ValueType>::prepareBoundaries(Configuration::Configuration const &config, Acquisition::Coordinates const &modelCoordinates,  Derivatives::Derivatives<ValueType> &derivatives, scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, KITGPI::ForwardSolver::BoundaryCondition::FreeSurface<ValueType> &FreeSurface,BoundaryCondition::ABS<ValueType> &DampingBoundary,BoundaryCondition::CPML<ValueType> &ConvPML)
{

    useFreeSurface = config.get<scai::IndexType>("FreeSurface");

    /* Prepare Free Surface */
    if (useFreeSurface==1) {
        FreeSurface.init(dist, derivatives, config.get<scai::IndexType>("NX"), config.get<scai::IndexType>("NY"), config.get<scai::IndexType>("NZ"),modelCoordinates, config.get<ValueType>("DT"),config.get<ValueType>("DH"));
    }

    /* Prepare Damping Boundary */
    if (config.get<IndexType>("DampingBoundary") == 1) {
        useDampingBoundary = true;
        DampingBoundary.init(dist, ctx,modelCoordinates, config.get<scai::IndexType>("BoundaryWidth"), config.get<ValueType>("DampingCoeff"), useFreeSurface);
    }

    if (config.get<IndexType>("DampingBoundary") == 2) {
        useConvPML = true;
        ConvPML.init(dist, ctx,modelCoordinates, config.get<ValueType>("DT"), config.get<scai::IndexType>("DH"), config.get<scai::IndexType>("BoundaryWidth"), config.get<ValueType>("NPower"), config.get<ValueType>("KMaxCPML"), config.get<ValueType>("CenterFrequencyCPML"), config.get<ValueType>("VMaxCPML"), useFreeSurface);
    }
}

template class KITGPI::ForwardSolver::ForwardSolver<double>;
template class KITGPI::ForwardSolver::ForwardSolver<float>;


