#include "ForwardSolverEM.hpp"

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
scai::lama::DenseVector<ValueType> KITGPI::ForwardSolver::ForwardSolverEM<ValueType>::getAveragedCinv(scai::lama::Vector<ValueType> const &vecAvDielectricPermittivity, scai::lama::Vector<ValueType> const &vecAvElectricConductivity, ValueType DT)
{
    scai::lama::DenseVector<ValueType> averagedCinv;
    averagedCinv = 0.5 / vecAvDielectricPermittivity;
    averagedCinv *= vecAvElectricConductivity;
    averagedCinv *= DT;
    averagedCinv += 1;
    averagedCinv = 1 / averagedCinv;
    
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
scai::lama::DenseVector<ValueType> KITGPI::ForwardSolver::ForwardSolverEM<ValueType>::getAveragedCa(scai::lama::Vector<ValueType> const &vecAvDielectricPermittivity, scai::lama::Vector<ValueType> const &vecAvElectricConductivity, ValueType DT)
{
    scai::lama::DenseVector<ValueType> averagedCa;
    averagedCa = 0.5 / vecAvDielectricPermittivity;
    averagedCa *= vecAvElectricConductivity;
    averagedCa *= DT;
    averagedCa = 1 - averagedCa;
    averagedCa *= getAveragedCinv(vecAvDielectricPermittivity, vecAvElectricConductivity, DT);
    
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
scai::lama::DenseVector<ValueType> KITGPI::ForwardSolver::ForwardSolverEM<ValueType>::getAveragedCb(scai::lama::Vector<ValueType> const &vecAvDielectricPermittivity, scai::lama::Vector<ValueType> const &vecAvElectricConductivity, ValueType DT)
{
    scai::lama::DenseVector<ValueType> averagedCb;
    averagedCb = 1 / vecAvDielectricPermittivity;
    averagedCb *= getAveragedCinv(vecAvDielectricPermittivity, vecAvElectricConductivity, DT);
    
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
ValueType KITGPI::ForwardSolver::ForwardSolverEM<ValueType>::getCc(ValueType tauElectricDisplacement, ValueType DT)
{
    ValueType Cc;
    Cc = (1 - 0.5 * DT / tauElectricDisplacement) / (1 + 0.5 * DT / tauElectricDisplacement);
    
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
scai::lama::DenseVector<ValueType> KITGPI::ForwardSolver::ForwardSolverEM<ValueType>::getAveragedCd(scai::lama::Vector<ValueType> const &vecAvDielectricPermittivitystatic, scai::lama::Vector<ValueType> const &vecAvTauDielectricPermittivity, scai::IndexType numRelaxationMechanisms, ValueType tauElectricDisplacement, ValueType DT)
{
    scai::lama::DenseVector<ValueType> averagedCd;
    ValueType tempValue;
    tempValue = 1 / (1 + 0.5 * DT / tauElectricDisplacement);
    tempValue /= (numRelaxationMechanisms * tauElectricDisplacement * tauElectricDisplacement);
    averagedCd = vecAvTauDielectricPermittivity * tempValue;
    averagedCd *= vecAvDielectricPermittivitystatic;
    averagedCd *= -DT; // please note here we include DT into the coefficient
    
    return averagedCd;
}

template class KITGPI::ForwardSolver::ForwardSolverEM<double>;
template class KITGPI::ForwardSolver::ForwardSolverEM<float>;
