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
 \param relaxationTime relaxtiontime of electric displacement
 \param DT time step
 
 \begin{align*}
  C_{c l} = \left( 1 - \frac{ \Delta t}{2 \tau_{Dl}} \right) \left( 1 + \frac{ \Delta t}{2 \tau_{Dl}} \right)^{-1}
 \end{align*}
 */
template <typename ValueType>
std::vector<ValueType> KITGPI::ForwardSolver::ForwardSolverEM<ValueType>::getCc(scai::IndexType numRelaxationMechanisms, std::vector<ValueType> relaxationTime, ValueType DT)
{
    std::vector<ValueType> Cc;
    for (int l=0; l<numRelaxationMechanisms; l++) {
        Cc.push_back((1 - 0.5 * DT / relaxationTime[l]) / (1 + 0.5 * DT / relaxationTime[l]));
    }
    
    return Cc;
}

/*! \brief calculate averaged AveragedCd
 *
 \param vecAvDielectricPermittivitystatic static Dielectric Permittivity
 \param vecAvTauDielectricPermittivity relaxtiontime of Dielectric Permittivity
 \param relaxationTime relaxtiontime of electric displacement
 \param DT time step
 
 \begin{align*}
  C_{d i l} = -\varepsilon^s_{i} \frac{\tau_{\varepsilon i l}}{L \tau_{Dl}^2} \left( 1 + \frac{ \Delta t}{2 \tau_{Dl}} \right)^{-1}
 \end{align*}
 */
template <typename ValueType>
std::vector<scai::lama::DenseVector<ValueType>> KITGPI::ForwardSolver::ForwardSolverEM<ValueType>::getAveragedCd(scai::lama::Vector<ValueType> const &vecAvDielectricPermittivitystatic, scai::lama::Vector<ValueType> const &vecAvTauDielectricPermittivity, scai::IndexType numRelaxationMechanisms, std::vector<ValueType> relaxationTime, ValueType DT)
{
    std::vector<scai::lama::DenseVector<ValueType>> averagedCd;
    ValueType tempValue;
    scai::lama::DenseVector<ValueType> temp;
    for (int l=0; l<numRelaxationMechanisms; l++) {
        tempValue = 1 / (1 + 0.5 * DT / relaxationTime[l]);
        tempValue /= (numRelaxationMechanisms * relaxationTime[l] * relaxationTime[l]);
        temp = vecAvTauDielectricPermittivity * tempValue;
        temp *= vecAvDielectricPermittivitystatic;
        temp *= -DT; // please note here we include DT into the coefficient
        averagedCd.push_back(temp);
    }
    
    return averagedCd;
}

/*! \brief Get const reference to electricConductivityEffectiveOptical
 */
template <typename ValueType>
scai::lama::DenseVector<ValueType> const KITGPI::ForwardSolver::ForwardSolverEM<ValueType>::getElectricConductivityEffectiveOptical(scai::lama::Vector<ValueType> const &dielectricPermittivityStatic, scai::lama::Vector<ValueType> const &electricConductivityStatic, scai::lama::Vector<ValueType> const &tauDielectricPermittivityAverage, scai::IndexType numRelaxationMechanisms, std::vector<ValueType> relaxationTime)
{  
    scai::lama::DenseVector<ValueType> electricConductivityEffectiveOptical;
    ValueType sum = 0;
    for (int l=0; l<numRelaxationMechanisms; l++) {
        sum += 1.0 / relaxationTime[l];
    }
    sum /= numRelaxationMechanisms;
    electricConductivityEffectiveOptical = tauDielectricPermittivityAverage * sum;    
    electricConductivityEffectiveOptical *= dielectricPermittivityStatic;
    electricConductivityEffectiveOptical += electricConductivityStatic;
    
    return (electricConductivityEffectiveOptical);
}

/*! \brief Get const reference to dielectricPermittivityEffectiveOptical
 */
template <typename ValueType>
scai::lama::DenseVector<ValueType> const KITGPI::ForwardSolver::ForwardSolverEM<ValueType>::getDielectricPermittivityEffectiveOptical(scai::lama::Vector<ValueType> const &dielectricPermittivityStatic, scai::lama::Vector<ValueType> const &electricConductivityStatic, scai::lama::Vector<ValueType> const &tauDielectricPermittivityAverage, scai::lama::Vector<ValueType> const &tauElectricConductivityAverage, ValueType DielectricPermittivityVacuum)
{
    scai::lama::DenseVector<ValueType> dielectricPermittivityEffectiveOptical;
    scai::lama::DenseVector<ValueType> temp;
    
    dielectricPermittivityEffectiveOptical = 1 - tauDielectricPermittivityAverage;      
    dielectricPermittivityEffectiveOptical *= dielectricPermittivityStatic;
    temp = electricConductivityStatic * tauElectricConductivityAverage;
    dielectricPermittivityEffectiveOptical += temp;
    
    // optical dielectric permittivity should be larger than the dielectric permittivity of vacuum
    Common::searchAndReplace<ValueType>(dielectricPermittivityEffectiveOptical, DielectricPermittivityVacuum, DielectricPermittivityVacuum, 1); 
    
    return (dielectricPermittivityEffectiveOptical);
}

template class KITGPI::ForwardSolver::ForwardSolverEM<double>;
template class KITGPI::ForwardSolver::ForwardSolverEM<float>;
