#include "ForwardSolverSeismic.hpp"

/*! \brief calculate averaged AveragedCinv which will be used in Ca, Cb
 *
 \param vecAvDielectricPermittivity Averaged dielectricPermittivity vector
 \param vecAvElectricConductivity Averaged ElectricConductivity vector
 \param DT time step
 
 \begin{align*}
  &C_{inv} = \left( 1 + \frac{\sigma_{ei} \Delta t}{2 \varepsilon_{ei}} \right)^{-1}, 
 \end{align*}
 */
template <typename ValueType>
scai::lama::DenseVector<ValueType> KITGPI::ForwardSolver::ForwardSolverSeismic<ValueType>::getAveragedCinv(scai::lama::Vector<ValueType> const &vecAvDielectricPermittivity, scai::lama::Vector<ValueType> const &vecAvElectricConductivity, ValueType DT)
{
    scai::lama::DenseVector<ValueType> averagedCinv;
    COMMON_THROWEXCEPTION("There is no averagedCinv in an Seismic modelling")
    
    return averagedCinv;
}

/*! \brief calculate averaged AveragedCa
 *
 \param vecAvDielectricPermittivity Averaged dielectricPermittivity vector
 \param vecAvElectricConductivity Averaged ElectricConductivity vector
 \param DT time step
 
 \begin{align*}
  &C_{ai} = \left( 1 - \frac{\sigma_{ei} \Delta t}{2 \varepsilon_{ei}} \right) \left( 1 + \frac{\sigma_{ei} \Delta t}{2 \varepsilon_{ei}} \right)^{-1}
 \end{align*}
 */
template <typename ValueType>
scai::lama::DenseVector<ValueType> KITGPI::ForwardSolver::ForwardSolverSeismic<ValueType>::getAveragedCa(scai::lama::Vector<ValueType> const &vecAvDielectricPermittivity, scai::lama::Vector<ValueType> const &vecAvElectricConductivity, ValueType DT)
{
    scai::lama::DenseVector<ValueType> averagedCa;
    COMMON_THROWEXCEPTION("There is no averagedCa in an Seismic modelling")
    
    return averagedCa;
}

/*! \brief calculate averaged AveragedCb
 *
 \param vecAvDielectricPermittivity Averaged dielectricPermittivity vector
 \param vecAvElectricConductivity Averaged ElectricConductivity vector
 \param DT time step
 
 \begin{align*}
  C_{bi} = \frac{1}{\varepsilon_{ei}} \left( 1 + \frac{\sigma_{ei} \Delta t}{2 \varepsilon_{ei}} \right)^{-1}, \quad
  i = x,y,z\\
 \end{align*}
 */
template <typename ValueType>
scai::lama::DenseVector<ValueType> KITGPI::ForwardSolver::ForwardSolverSeismic<ValueType>::getAveragedCb(scai::lama::Vector<ValueType> const &vecAvDielectricPermittivity, scai::lama::Vector<ValueType> const &vecAvElectricConductivity, ValueType DT)
{
    scai::lama::DenseVector<ValueType> averagedCb;
    COMMON_THROWEXCEPTION("There is no averagedCb in an Seismic modelling")
    
    return averagedCb;
}

/*! \brief calculate Cc
 *
 \param tauElectricDisplacement relaxtiontime of electric displacement
 \param DT time step
 
 \begin{align*}
  C_{c l} = \left( 1 - \frac{ \Delta t}{2 \tau_{Dl}} \right) \left( 1 + \frac{ \Delta t}{2 \tau_{Dl}} \right)^{-1}
 \end{align*}
 */
template <typename ValueType>
ValueType KITGPI::ForwardSolver::ForwardSolverSeismic<ValueType>::getCc(ValueType tauElectricDisplacement, ValueType DT)
{
    ValueType Cc;
    COMMON_THROWEXCEPTION("There is no Cc in an Seismic modelling")
    
    return Cc;
}

/*! \brief calculate averaged AveragedCd
 *
 \param vecAvDielectricPermittivitystatic static Dielectric Permittivity
 \param vecAvTauDielectricPermittivity relaxtiontime of Dielectric Permittivity
 \param tauElectricDisplacement relaxtiontime of electric displacement
 \param DT time step
 
 \begin{align*}
  C_{d i l} = -\varepsilon^s_{i} \frac{\tau_{\varepsilon i l}}{L \tau_{Dl}^2} \left( 1 + \frac{ \Delta t}{2 \tau_{Dl}} \right)^{-1}
 \end{align*}
 */
template <typename ValueType>
scai::lama::DenseVector<ValueType> KITGPI::ForwardSolver::ForwardSolverSeismic<ValueType>::getAveragedCd(scai::lama::Vector<ValueType> const &vecAvDielectricPermittivitystatic, scai::lama::Vector<ValueType> const &vecAvTauDielectricPermittivity, scai::IndexType numRelaxationMechanisms, ValueType tauElectricDisplacement, ValueType DT)
{
    scai::lama::DenseVector<ValueType> averagedCd;
    COMMON_THROWEXCEPTION("There is no averagedCd in an Seismic modelling")
    
    return averagedCd;
}

template class KITGPI::ForwardSolver::ForwardSolverSeismic<double>;
template class KITGPI::ForwardSolver::ForwardSolverSeismic<float>;
