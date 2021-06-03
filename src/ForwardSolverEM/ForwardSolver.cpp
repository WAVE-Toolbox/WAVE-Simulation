#include "ForwardSolver.hpp"

/*! \brief memory estimation wrapper for the boundarie conditions
             *
             *
             \param config Configuration
             \param dist Distribution of the wave fields
             \param modelCoordinates Coordinate class, which eg. maps 3D coordinates to 1D model indices
             */
template <typename ValueType>
ValueType KITGPI::ForwardSolver::ForwardSolverEM<ValueType>::estimateBoundaryMemory(Configuration::Configuration const &config, scai::dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> const &modelCoordinates, BoundaryCondition::ABSEM<ValueType> &DampingBoundary, BoundaryCondition::CPMLEM<ValueType> &ConvPML)
{
    if (config.get<IndexType>("DampingBoundary") == 1) {
        return (DampingBoundary.estimateMemory(config.get<scai::IndexType>("BoundaryWidth"), config.get<scai::IndexType>("FreeSurface"), dist, modelCoordinates));
    } else if (config.get<IndexType>("DampingBoundary") == 2) {
        return (ConvPML.estimateMemory(config.get<scai::IndexType>("BoundaryWidth"), config.get<scai::IndexType>("FreeSurface"), dist, modelCoordinates));
    }
    return 0;
}

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
void KITGPI::ForwardSolver::ForwardSolverEM<ValueType>::prepareBoundaries(Configuration::Configuration const &config, Acquisition::Coordinates<ValueType> const &modelCoordinates, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> &derivatives, scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, KITGPI::ForwardSolver::BoundaryCondition::FreeSurfaceEM<ValueType> &FreeSurface, BoundaryCondition::ABSEM<ValueType> &DampingBoundary, BoundaryCondition::CPMLEM<ValueType> &ConvPML)
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

/*! \brief calculate averaged AveragedCinv which will be used in Ca, Cb
 *
 \param vecAvDielectricPermittivityEM Averaged dielectricPermittivityEM vector
 \param vecAvConductivityEM Averaged ConductivityEM vector
 \param DT time step
 
 \begin{align*}
  &C_{inv} = \left( 1 + \frac{\sigma_{ei} \Delta t}{2 \varepsilon_{ei}} \right)^{-1}, 
 \end{align*}
 */
template <typename ValueType>
scai::lama::DenseVector<ValueType> KITGPI::ForwardSolver::ForwardSolverEM<ValueType>::getAveragedCinv(scai::lama::Vector<ValueType> const &vecAvDielectricPermittivityEM, scai::lama::Vector<ValueType> const &vecAvConductivityEM, ValueType DT)
{
    scai::lama::DenseVector<ValueType> averagedCinv;
    averagedCinv = 0.5 / vecAvDielectricPermittivityEM;
    averagedCinv *= vecAvConductivityEM;
    averagedCinv *= DT;
    averagedCinv += 1;
    averagedCinv = 1 / averagedCinv;
    
    return averagedCinv;
}

/*! \brief calculate averaged AveragedCa
 *
 \param vecAvDielectricPermittivityEM Averaged dielectricPermittivityEM vector
 \param vecAvConductivityEM Averaged ConductivityEM vector
 \param DT time step
 
 \begin{align*}
  &C_{ai} = \left( 1 - \frac{\sigma_{ei} \Delta t}{2 \varepsilon_{ei}} \right) \left( 1 + \frac{\sigma_{ei} \Delta t}{2 \varepsilon_{ei}} \right)^{-1}
 \end{align*}
 */
template <typename ValueType>
scai::lama::DenseVector<ValueType> KITGPI::ForwardSolver::ForwardSolverEM<ValueType>::getAveragedCa(scai::lama::Vector<ValueType> const &vecAvDielectricPermittivityEM, scai::lama::Vector<ValueType> const &vecAvConductivityEM, ValueType DT)
{
    scai::lama::DenseVector<ValueType> averagedCa;
    averagedCa = 0.5 / vecAvDielectricPermittivityEM;
    averagedCa *= vecAvConductivityEM;
    averagedCa *= DT;
    averagedCa = 1 - averagedCa;
    averagedCa *= getAveragedCinv(vecAvDielectricPermittivityEM, vecAvConductivityEM, DT);
    
    return averagedCa;
}

/*! \brief calculate averaged AveragedCb
 *
 \param vecAvDielectricPermittivityEM Averaged dielectricPermittivityEM vector
 \param vecAvConductivityEM Averaged ConductivityEM vector
 \param DT time step
 
 \begin{align*}
  C_{bi} = \frac{1}{\varepsilon_{ei}} \left( 1 + \frac{\sigma_{ei} \Delta t}{2 \varepsilon_{ei}} \right)^{-1}, \quad
  i = x,y,z\\
 \end{align*}
 */
template <typename ValueType>
scai::lama::DenseVector<ValueType> KITGPI::ForwardSolver::ForwardSolverEM<ValueType>::getAveragedCb(scai::lama::Vector<ValueType> const &vecAvDielectricPermittivityEM, scai::lama::Vector<ValueType> const &vecAvConductivityEM, ValueType DT)
{
    scai::lama::DenseVector<ValueType> averagedCb;
    averagedCb = 1 / vecAvDielectricPermittivityEM;
    averagedCb *= getAveragedCinv(vecAvDielectricPermittivityEM, vecAvConductivityEM, DT);
    
    return averagedCb;
}

/*! \brief calculate Cc
 *
 \param tauDisplacementEM relaxtiontime of electric displacement
 \param DT time step
 
 \begin{align*}
  C_{c l} = \left( 1 - \frac{ \Delta t}{2 \tau_{Dl}} \right) \left( 1 + \frac{ \Delta t}{2 \tau_{Dl}} \right)^{-1}
 \end{align*}
 */
template <typename ValueType>
ValueType KITGPI::ForwardSolver::ForwardSolverEM<ValueType>::getCc(ValueType tauDisplacementEM, ValueType DT)
{
    ValueType Cc;
    Cc = (1 - 0.5 * DT / tauDisplacementEM) / (1 + 0.5 * DT / tauDisplacementEM);
    
    return Cc;
}

/*! \brief calculate averaged AveragedCd
 *
 \param vecAvDielectricPermittivityEMstatic static Dielectric Permittivity
 \param vecAvTauDielectricPermittivityEM relaxtiontime of Dielectric Permittivity
 \param tauDisplacementEM relaxtiontime of electric displacement
 \param DT time step
 
 \begin{align*}
  C_{d i l} = -\varepsilon^s_{i} \frac{\tau_{\varepsilon i l}}{L \tau_{Dl}^2} \left( 1 + \frac{ \Delta t}{2 \tau_{Dl}} \right)^{-1}
 \end{align*}
 */
template <typename ValueType>
scai::lama::DenseVector<ValueType> KITGPI::ForwardSolver::ForwardSolverEM<ValueType>::getAveragedCd(scai::lama::Vector<ValueType> const &vecAvDielectricPermittivityEMstatic, scai::lama::Vector<ValueType> const &vecAvTauDielectricPermittivityEM, scai::IndexType numRelaxationMechanisms, ValueType tauDisplacementEM, ValueType DT)
{
    scai::lama::DenseVector<ValueType> averagedCd;
    ValueType tempValue;
    tempValue = 1 / (1 + 0.5 * DT / tauDisplacementEM);
    tempValue /= (numRelaxationMechanisms * tauDisplacementEM * tauDisplacementEM);
    averagedCd = vecAvTauDielectricPermittivityEM * tempValue;
    averagedCd *= vecAvDielectricPermittivityEMstatic;
    averagedCd *= -DT; // please note here we include DT into the coefficient
    
    return averagedCd;
}

template class KITGPI::ForwardSolver::ForwardSolverEM<double>;
template class KITGPI::ForwardSolver::ForwardSolverEM<float>;
