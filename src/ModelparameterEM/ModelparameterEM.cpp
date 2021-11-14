#include "ModelparameterEM.hpp"
#include "../IO/IO.hpp"

using namespace scai;
using namespace KITGPI;

/*! \brief Prepare modelparameter for inversion
 */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterEM<ValueType>::prepareForInversion(Configuration::Configuration const &config, scai::dmemo::CommunicatorPtr comm)
{
    HOST_PRINT(comm, "", "Preparation of the model parameters inversion\n");
    this->setParameterisation(config.getAndCatch("parameterisation", 0));
    this->setInversionType(config.getAndCatch("inversionType", 1));
    this->setGradientType(config.getAndCatch("gradientType", 0));
    this->setDecomposeType(config.getAndCatch("decomposeType", 0));
    calcElectricConductivityReference(config.get<ValueType>("CenterFrequencyCPML"));
    setArchieFactors(config.get<ValueType>("aArchie"), config.get<ValueType>("mArchie"), config.get<ValueType>("nArchie"));
    HOST_PRINT(comm, "", "Model ready!\n\n");
}

/*! \brief Set relaxation frequency
 */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterEM<ValueType>::setRelaxationFrequency(std::vector<ValueType> const setRelaxationFrequency)
{
    dirtyFlagAveraging = true;    // If modulus will be changed, averaging needs to be redone
    dirtyFlagVelocivityEM = true; // the velocity vector is now dirty
    dirtyFlagDielectricPermittivityEffectiveOptical = true; // the visco Modulus vector is now dirty
    dirtyFlagElectricConductivityEffectiveOptical = true; // the visco Modulus vector is now dirty
    relaxationFrequency = setRelaxationFrequency;
}

/*! \brief Set number of Relaxation Mechanisms
 */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterEM<ValueType>::setNumRelaxationMechanisms(IndexType const setNumRelaxationMechanisms)
{
    dirtyFlagAveraging = true;    // If modulus will be changed, averaging needs to be redone
    dirtyFlagVelocivityEM = true; // the velocity vector is now dirty
    dirtyFlagDielectricPermittivityEffectiveOptical = true; // the visco Modulus vector is now dirty
    dirtyFlagElectricConductivityEffectiveOptical = true; // the visco Modulus vector is now dirty
    numRelaxationMechanisms = setNumRelaxationMechanisms;
}

/*! \brief Set method for Archie factors a, m, n */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterEM<ValueType>::setArchieFactors(ValueType const &aArchie_const, ValueType const &mArchie_const, ValueType const &nArchie_const)
{        
    aArchie = aArchie_const;
    mArchie = mArchie_const;
    nArchie = nArchie_const;
}

/*! \brief Getter method for Archie factors a */
template <typename ValueType>
ValueType const KITGPI::Modelparameter::ModelparameterEM<ValueType>::getArchie_a() const
{        
    return (aArchie);
}

/*! \brief Getter method for Archie factors m */
template <typename ValueType>
ValueType const KITGPI::Modelparameter::ModelparameterEM<ValueType>::getArchie_m() const
{        
    return (mArchie);
}

/*! \brief Getter method for Archie factors n */
template <typename ValueType>
ValueType const KITGPI::Modelparameter::ModelparameterEM<ValueType>::getArchie_n() const
{        
    return (nArchie);
}

/*! \brief Getter method for DielectricPermittivityVacuum */
template <typename ValueType>
ValueType const KITGPI::Modelparameter::ModelparameterEM<ValueType>::getDielectricPermittivityVacuum() const
{        
    return (DielectricPermittivityVacuum);
}

/*! \brief Getter method for ElectricConductivityReference */
template <typename ValueType>
ValueType const KITGPI::Modelparameter::ModelparameterEM<ValueType>::getElectricConductivityReference() const
{        
    return (ElectricConductivityReference);
}

/*! \brief Getter method for TauDielectricPermittivityReference */
template <typename ValueType>
ValueType const KITGPI::Modelparameter::ModelparameterEM<ValueType>::getTauDielectricPermittivityReference() const
{        
    return (TauDielectricPermittivityReference);
}

/*! \brief Getter method for TauElectricConductivityReference */
template <typename ValueType>
ValueType const KITGPI::Modelparameter::ModelparameterEM<ValueType>::getTauElectricConductivityReference() const
{        
    ValueType TauElectricConductivityReference = 10 * TauDielectricPermittivityReference * this->getTauElectricDisplacement(); 
    return (TauElectricConductivityReference);
}

/*! \brief Getter method for RelativeDielectricPermittivityWater */
template <typename ValueType>
ValueType const KITGPI::Modelparameter::ModelparameterEM<ValueType>::getRelativeDielectricPermittivityWater() const
{        
    return (RelativeDielectricPermittivityWater);
}

/*! \brief Getter method for RelativeDielectricPermittivityVacuum */
template <typename ValueType>
ValueType const KITGPI::Modelparameter::ModelparameterEM<ValueType>::getRelativeDielectricPermittivityVacuum() const
{        
    return (RelativeDielectricPermittivityVacuum);
}

/*! \brief Write RockMatrix model to an external file
 *base filename of the model
 \param fileFormat Output file format 1=mtx 2=lmf
 */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterEM<ValueType>::writeRockMatrixParameter(std::string filename, scai::IndexType fileFormat)
{
    IO::writeVector(electricConductivityWater, filename + ".sigmaEMw", fileFormat);
    IO::writeVector(relativeDieletricPeimittivityRockMatrix, filename + ".epsilonEMrma", fileFormat);
};

/*! \brief calculate electricConductivityWater and relativeDielectricPermittivityRockMatrix from porosity and saturation */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterEM<ValueType>::calcRockMatrixParameter(Configuration::Configuration const &config)
{            
    scai::lama::DenseVector<ValueType> temp1;
    scai::lama::DenseVector<ValueType> temp2;
    scai::lama::DenseVector<ValueType> electricConductivityWatertemp;
    scai::lama::DenseVector<ValueType> relativeDielectricPermittivityRockMatrixtemp;  
    
    electricConductivityWatertemp = this->getElectricConductivity();  
    relativeDielectricPermittivityRockMatrixtemp = this->getDielectricPermittivity();      
    relativeDielectricPermittivityRockMatrixtemp /= DielectricPermittivityVacuum;    
    relativeDielectricPermittivityRockMatrixtemp = scai::lama::sqrt(relativeDielectricPermittivityRockMatrixtemp); 
    
    // Based on Archie equation
    temp1 = scai::lama::pow(porosity, -mArchie);
    temp2 = scai::lama::pow(saturation, -nArchie);
    electricConductivityWatertemp *= aArchie;
    electricConductivityWatertemp *= temp1;
    electricConductivityWatertemp *= temp2;
    
    Common::replaceInvalid<ValueType>(electricConductivityWatertemp, 0.0); // in case of air
    Common::searchAndReplace<ValueType>(electricConductivityWatertemp, 0, 0, 1);
    
    // Based on complex refractive index model (CRIM)    
    temp1 = saturation;
    temp1 *= sqrt(RelativeDielectricPermittivityWater);
    temp2 = 1 - saturation;
    temp2 *= sqrt(RelativeDielectricPermittivityVacuum);
    temp1 += temp2;
    temp1 *= porosity;  
    relativeDielectricPermittivityRockMatrixtemp -= temp1;
    temp1 = 1 - porosity;
    temp1 = 1 / temp1;
    relativeDielectricPermittivityRockMatrixtemp *= temp1;
    relativeDielectricPermittivityRockMatrixtemp = scai::lama::pow(relativeDielectricPermittivityRockMatrixtemp, 2);
    
    Common::replaceInvalid<ValueType>(relativeDielectricPermittivityRockMatrixtemp, RelativeDielectricPermittivityVacuum); // in case of air
    Common::searchAndReplace<ValueType>(relativeDielectricPermittivityRockMatrixtemp, RelativeDielectricPermittivityVacuum, RelativeDielectricPermittivityVacuum, 1);
    Common::searchAndReplace<ValueType>(relativeDielectricPermittivityRockMatrixtemp, RelativeDielectricPermittivityWater, RelativeDielectricPermittivityWater, 2);
    
    this->setElectricConductivityWater(electricConductivityWatertemp);
    this->setRelativeDieletricPeimittivityRockMatrix(relativeDielectricPermittivityRockMatrixtemp);
}

/*! \brief calculate ElectricConductivityReference from the source center frequency and DielectricPermittivityVacuum */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterEM<ValueType>::calcElectricConductivityReference(ValueType const CenterFrequencyCPML)
{
    // please note here ElectricConductivityReference is calculated by f_min=0.5*f_c, which is corresponding to the scaling factor = 0.5 in Lavoue et al (2014).
    ElectricConductivityReference = 0.5 * 2 * M_PI * CenterFrequencyCPML * DielectricPermittivityVacuum; 
}

/*! \brief transform porosity and saturation to EM parameters */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterEM<ValueType>::calcWaveModulusFromPetrophysics()
{
    scai::lama::DenseVector<ValueType> temp1;
    scai::lama::DenseVector<ValueType> temp2;
    scai::lama::DenseVector<ValueType> electricConductivitytemp;
    scai::lama::DenseVector<ValueType> dielectricPermittivitytemp;  
    
    electricConductivitytemp = this->getElectricConductivityWater();  
    dielectricPermittivitytemp = this->getRelativeDieletricPeimittivityRockMatrix();   
    
    // Based on Archie equation
    temp1 = scai::lama::pow(porosity, mArchie);
    temp2 = scai::lama::pow(saturation, nArchie);
    electricConductivitytemp /= aArchie;
    electricConductivitytemp *= temp1;
    electricConductivitytemp *= temp2;
    
    // Based on complex refractive index model (CRIM)    
    temp1 = saturation;
    temp1 *= sqrt(RelativeDielectricPermittivityWater);
    temp2 = 1 - saturation;
    temp2 *= sqrt(RelativeDielectricPermittivityVacuum);
    temp1 += temp2;
    temp1 *= porosity;  
    dielectricPermittivitytemp = scai::lama::sqrt(dielectricPermittivitytemp);
    temp2 = 1 - porosity;
    dielectricPermittivitytemp *= temp2;
    dielectricPermittivitytemp += temp1;    
    dielectricPermittivitytemp = scai::lama::pow(dielectricPermittivitytemp, 2);
    
    Common::searchAndReplace<ValueType>(electricConductivitytemp, 0, 0, 1);
    Common::searchAndReplace<ValueType>(dielectricPermittivitytemp, RelativeDielectricPermittivityVacuum, RelativeDielectricPermittivityVacuum, 1);
    Common::searchAndReplace<ValueType>(dielectricPermittivitytemp, RelativeDielectricPermittivityWater, RelativeDielectricPermittivityWater, 2);
    
    dielectricPermittivitytemp *= DielectricPermittivityVacuum;    
    
    this->setElectricConductivity(electricConductivitytemp);
    this->setDielectricPermittivity(dielectricPermittivitytemp);
}

/*! \brief transform EM parameters to porosity and saturation  */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterEM<ValueType>::calcPetrophysicsFromWaveModulus()
{
    scai::lama::DenseVector<ValueType> temp1;
    scai::lama::DenseVector<ValueType> temp2;
    scai::lama::DenseVector<ValueType> porositytemp;
    scai::lama::DenseVector<ValueType> dielectricPermittivitytemp;  
    scai::lama::DenseVector<ValueType> relativeDieletricPeimittivityRockMatrix;  
    
    relativeDieletricPeimittivityRockMatrix = this->getRelativeDieletricPeimittivityRockMatrix();  
    dielectricPermittivitytemp = this->getDielectricPermittivity();   
    dielectricPermittivitytemp /= DielectricPermittivityVacuum;
    porositytemp = this->getPorosity();   
        
    // Based on complex refractive index model (CRIM)        
    temp1 = 1 - porositytemp;     
    temp2 = scai::lama::sqrt(relativeDieletricPeimittivityRockMatrix);
    temp2 *= temp1;
    temp1 = scai::lama::sqrt(dielectricPermittivitytemp);
    temp1 -= temp2;
    temp1 /= porositytemp;
    temp1 -= sqrt(RelativeDielectricPermittivityVacuum);
    temp1 /= (sqrt(RelativeDielectricPermittivityWater) - sqrt(RelativeDielectricPermittivityVacuum));
    Common::replaceInvalid<ValueType>(temp1, 0); // in case of air
    
    this->setSaturation(temp1);      
}

/*! \brief Calculate a modulus from velocity
 *
 *  Calculates Modulus = 1 / (pow(Velocity,2) *  MagneticPermeability)
 *
 \param vecVelocity Velocity-Vector which will be used in the calculation
 \param vecMagneticPermeability MagneticPermeability-Vector which will be used in the calculation
 \param vecDielectricPermittivity Modulus-Vector which is calculated
 */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterEM<ValueType>::calcModulusFromVelocity(scai::lama::Vector<ValueType> &vecVelocity, scai::lama::Vector<ValueType> &vecMagneticPermeability, scai::lama::Vector<ValueType> &vecDielectricPermittivity)
{
    vecDielectricPermittivity = 1 / vecMagneticPermeability;
    vecDielectricPermittivity /= vecVelocity;
    vecDielectricPermittivity /= vecVelocity;
    scai::lama::DenseVector<ValueType> temp;
    temp = vecDielectricPermittivity;
    Common::replaceInvalid<ValueType>(temp, 0.0); // in case of Inf divided by zero.
    vecDielectricPermittivity = temp;
};

/*! \brief Calculate velocities from a modulus
 *
 *  Calculates Velocity = 1 / sqrt( Modulu * MagneticPermeability )
 *
 \param vecDielectricPermittivity Modulus-Vector which will be used in the calculation
 \param vecMagneticPermeability MagneticPermeability-Vector which will be used in the calculation
 \param vecVelocity Velocity-Vector which is calculated
 */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterEM<ValueType>::calcVelocityFromModulus(scai::lama::Vector<ValueType> &vecDielectricPermittivity, scai::lama::Vector<ValueType> &vecMagneticPermeability, scai::lama::Vector<ValueType> &vecVelocity)
{
    vecVelocity = vecDielectricPermittivity * vecMagneticPermeability; 
    vecVelocity = lama::sqrt(vecVelocity);    
    vecVelocity = 1 / vecVelocity;
    scai::lama::DenseVector<ValueType> temp;
    temp = vecVelocity;
    Common::replaceInvalid<ValueType>(temp, 0.0); // in case of Inf divided by zero.
    vecVelocity = temp;
};

/*! \brief Get const reference to EM-wave velocity
 * 
 * If EM-Wave velocity is dirty eg. because the EM-Wave velocity was modified, EM-Wave velocity will be calculated from magneticPermeability and dielectricPermittivity.
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getVelocityEM()
{
    // If the modulus is dirty, then recalculate
    if (dirtyFlagVelocivityEM) {
        HOST_PRINT(dielectricPermittivity.getDistributionPtr()->getCommunicatorPtr(), "", "EM-Wave velocity will be calculated from magneticPermeability and dielectricPermittivity\n");
        calcVelocityFromModulus(dielectricPermittivity, magneticPermeability, velocivityEM);
        dirtyFlagVelocivityEM = false;
    }
    return (velocivityEM);
}

/*! \brief Get const reference to EM-wave velocity
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getVelocityEM() const
{
    SCAI_ASSERT(dirtyFlagVelocivityEM == false, "EM-Wave velocity has to be recalculated! ");
    return (velocivityEM);
}

/*! \brief Get const reference to magneticPermeability model parameter
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getMagneticPermeability() const
{
    return (magneticPermeability);
}

/*! \brief Set magneticPermeability model parameter
 */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterEM<ValueType>::setMagneticPermeability(scai::lama::Vector<ValueType> const &setMagneticPermeability)
{
    dirtyFlagAveraging = true;      // If magneticPermeability will be changed, averaging needs to be redone
    dirtyFlagVelocivityEM = true;
    magneticPermeability = setMagneticPermeability;
}

/*! \brief Get const reference to electricConductivity
 *
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getElectricConductivity() const
{
    return (electricConductivity);
}

/*! \brief Set electricConductivity
 */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterEM<ValueType>::setElectricConductivity(scai::lama::Vector<ValueType> const &setElectricConductivity)
{
    dirtyFlagAveraging = true;    // If modulus will be changed, averaging needs to be redone
    dirtyFlagVelocivityEM = true;   // the velocity vector is now dirty
    dirtyFlagElectricConductivityEffectiveOptical = true; // the visco Modulus vector is now dirty
    dirtyFlagDielectricPermittivityEffectiveOptical = true; // the visco Modulus vector is now dirty
    electricConductivity = setElectricConductivity;
}

/*! \brief Get const reference to dielectricPermittivity
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getDielectricPermittivity() const
{
    return (dielectricPermittivity);
}

/*! \brief Set dielectricPermittivity
 */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterEM<ValueType>::setDielectricPermittivity(scai::lama::Vector<ValueType> const &setDielectricPermittivity)
{
    dirtyFlagAveraging = true;    // If modulus will be changed, averaging needs to be redone
    dirtyFlagVelocivityEM = true; // the velocity vector is now dirty
    dirtyFlagElectricConductivityEffectiveOptical = true; // the visco Modulus vector is now dirty
    dirtyFlagDielectricPermittivityEffectiveOptical = true; // the visco Modulus vector is now dirty
    dielectricPermittivity = setDielectricPermittivity;
}

/*! \brief Get const reference to electricConductivityEffectiveOptical
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getElectricConductivityEffectiveOptical()
{
    // If the modulus is dirty, then recalculate
    if (dirtyFlagElectricConductivityEffectiveOptical) {
        HOST_PRINT(electricConductivity.getDistributionPtr()->getCommunicatorPtr(), "electricConductivityEffectiveOptical will be calculated from\n dielectricPermittivity, electricConductivity, tauDielectricPermittivity\n, tauElectricDisplacement and numRelaxationMechanisms \n");              
        electricConductivityEffectiveOptical = tauDielectricPermittivity / tauElectricDisplacement;    
        electricConductivityEffectiveOptical *= numRelaxationMechanisms; // numRelaxationMechanisms = 1
        electricConductivityEffectiveOptical /= numRelaxationMechanisms;
        electricConductivityEffectiveOptical *= dielectricPermittivity;
        electricConductivityEffectiveOptical += electricConductivity;
        dirtyFlagElectricConductivityEffectiveOptical = false;
    }
    
    return (electricConductivityEffectiveOptical);
}
/*! \brief Get const reference to electricConductivityEffectiveOptical
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getElectricConductivityEffectiveOptical() const
{
    SCAI_ASSERT(dirtyFlagElectricConductivityEffectiveOptical == false, "electricConductivityEffectiveOptical has to be recalculated! ");
    
    return (electricConductivityEffectiveOptical);
}

/*! \brief Set electricConductivityEffectiveOptical model parameter
 */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterEM<ValueType>::setElectricConductivityEffectiveOptical(scai::lama::Vector<ValueType> const &setElectricConductivityEffectiveOptical)
{
    dirtyFlagAveraging = true;    // If modulus will be changed, averaging needs to be redone
    electricConductivityEffectiveOptical = setElectricConductivityEffectiveOptical;
}

/*! \brief Get const reference to dielectricPermittivityEffectiveOptical
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getDielectricPermittivityEffectiveOptical()
{
    // If the modulus is dirty, then recalculate
    if (dirtyFlagDielectricPermittivityEffectiveOptical) {
        HOST_PRINT(dielectricPermittivity.getDistributionPtr()->getCommunicatorPtr(), "dielectricPermittivityEffectiveOptical will be calculated from\n dielectricPermittivity, electricConductivity, tauDielectricPermittivity\n, tauElectricConductivity and numRelaxationMechanisms \n");
        scai::lama::DenseVector<ValueType> temp;
        temp = numRelaxationMechanisms * tauDielectricPermittivity;  // numRelaxationMechanisms = 1
        temp /= numRelaxationMechanisms;
        dielectricPermittivityEffectiveOptical = 1 - temp;      
        dielectricPermittivityEffectiveOptical *= dielectricPermittivity;
        temp = electricConductivity * tauElectricConductivity;
        dielectricPermittivityEffectiveOptical += temp;
        
        // optical dielectric permittivity should be larger than the dielectric permittivity of vacuum
        Common::searchAndReplace<ValueType>(dielectricPermittivityEffectiveOptical, DielectricPermittivityVacuum, DielectricPermittivityVacuum, 1); 
        
        dirtyFlagDielectricPermittivityEffectiveOptical = false;
    }
    
    return (dielectricPermittivityEffectiveOptical);
}

/*! \brief Get const reference to dielectricPermittivityEffectiveOptical
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getDielectricPermittivityEffectiveOptical() const
{
    SCAI_ASSERT(dirtyFlagDielectricPermittivityEffectiveOptical == false, "dielectricPermittivityEffectiveOptical has to be recalculated! ");
    return (dielectricPermittivityEffectiveOptical);
}

/*! \brief Set dielectricPermittivityEffectiveOptical model parameter
 */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterEM<ValueType>::setDielectricPermittivityEffectiveOptical(scai::lama::Vector<ValueType> const &setDielectricPermittivityEffectiveOptical)
{
    dirtyFlagAveraging = true;    // If modulus will be changed, averaging needs to be redone
    dielectricPermittivityEffectiveOptical = setDielectricPermittivityEffectiveOptical;
}

/*! \brief Get const reference to tauElectricConductivity */
template <typename ValueType>
ValueType const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getTauElectricDisplacement()
{
    tauElectricDisplacement = 1.0 / (2.0 * M_PI * this->getRelaxationFrequency()[0]); 
    return (tauElectricDisplacement);
}

/*! \brief Get const reference to tauElectricConductivity */
template <typename ValueType>
ValueType const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getTauElectricDisplacement() const
{
    return (tauElectricDisplacement);
}

/*! \brief Get const reference to tauElectricConductivity */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getTauElectricConductivity() const
{
    return (tauElectricConductivity);
}

/*! \brief Set tauElectricConductivity model parameter */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterEM<ValueType>::setTauElectricConductivity(scai::lama::Vector<ValueType> const &setTauElectricConductivity)
{
    dirtyFlagAveraging = true;    // If modulus will be changed, averaging needs to be redone
    dirtyFlagVelocivityEM = true; // the velocity vector is now dirty
    dirtyFlagDielectricPermittivityEffectiveOptical = true; // the visco Modulus vector is now dirty
    tauElectricConductivity = setTauElectricConductivity;
}

/*! \brief Get const reference to tauDielectricPermittivity */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getTauDielectricPermittivity() const
{
    return (tauDielectricPermittivity);
}

/*! \brief Set tauDielectricPermittivity model parameter */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterEM<ValueType>::setTauDielectricPermittivity(scai::lama::Vector<ValueType> const &setTauDielectricPermittivity)
{
    dirtyFlagAveraging = true;    // If modulus will be changed, averaging needs to be redone
    dirtyFlagVelocivityEM = true; // the velocity vector is now dirty
    dirtyFlagElectricConductivityEffectiveOptical = true; // the visco Modulus vector is now dirty
    dirtyFlagDielectricPermittivityEffectiveOptical = true; // the visco Modulus vector is now dirty
    tauDielectricPermittivity = setTauDielectricPermittivity;
}

/*! \brief calculate reflectivity from permittivity
 */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterEM<ValueType>::calcReflectivity(Acquisition::Coordinates<ValueType> const &modelCoordinates, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives, ValueType DT)
{
    scai::lama::DenseVector<ValueType> impedance;
    scai::lama::DenseVector<ValueType> impedanceAverageY;
    IndexType spatialFDorder = derivatives.getSpatialFDorder();
    IndexType NX = modelCoordinates.getNX();
    IndexType NY = modelCoordinates.getNY();
    IndexType NZ = modelCoordinates.getNZ();
    auto dist = dielectricPermittivity.getDistributionPtr();
    auto ctx = dielectricPermittivity.getContextPtr();
    auto const &Dyf = derivatives.getDyf();
    this->calcAverageMatrixY(modelCoordinates, dist);
    averageMatrixY.setContextPtr(ctx);
    impedance = scai::lama::sqrt(dielectricPermittivity);
    impedance = 1 / impedance;
    this->calcAveragedParameter(impedance, impedanceAverageY, averageMatrixY);
    averageMatrixY.purge();
    impedanceAverageY *= 2;
    reflectivity = Dyf * impedance;
    reflectivity *= modelCoordinates.getDH() / DT;
    reflectivity /= impedanceAverageY;
    for (IndexType i=0; i < spatialFDorder * NX * NZ / 2; i++) {
        reflectivity[i] = 0.0;
        reflectivity[NX*NY*NZ-i-1] = 0.0;
    }
}

/*! \brief Get const reference to electricConductivityWater
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getElectricConductivityWater() const
{
    return (electricConductivityWater);
}

/*! \brief Set electricConductivityWater
 */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterEM<ValueType>::setElectricConductivityWater(scai::lama::Vector<ValueType> const &setElectricConductivityWater)
{
    electricConductivityWater = setElectricConductivityWater;
}

/*! \brief Get const reference to relativeDieletricPeimittivityRockMatrix
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getRelativeDieletricPeimittivityRockMatrix() const
{
    return (relativeDieletricPeimittivityRockMatrix);
}

/*! \brief Set relativeDieletricPeimittivityRockMatrix
 */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterEM<ValueType>::setRelativeDieletricPeimittivityRockMatrix(scai::lama::Vector<ValueType> const &setRelativeDieletricPeimittivityRockMatrix)
{
    relativeDieletricPeimittivityRockMatrix = setRelativeDieletricPeimittivityRockMatrix;
}

/*! \brief Get const reference to averaged magneticPermeability in yz-plane
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getInverseMagneticPermeabilityAverageYZ()
{
    // If Averaging is outdated or has to be calculated for the first time, then recalculate averaging
    if (dirtyFlagAveraging == true) {
        calculateAveraging();
    }
    return (inverseMagneticPermeabilityAverageYZ);
}

/*! \brief Get const reference to averaged magneticPermeability in yz-plane
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getInverseMagneticPermeabilityAverageYZ() const
{
    SCAI_ASSERT(dirtyFlagAveraging == false, "Averaging has to be recalculated! ");
    return (inverseMagneticPermeabilityAverageYZ);
}

/*! \brief Get const reference to averaged magneticPermeability in xz-plane
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getInverseMagneticPermeabilityAverageXZ()
{
    // If Averaging is outdated or has to be calculated for the first time, then recalculate averaging
    if (dirtyFlagAveraging == true) {
        calculateAveraging();
    }
    return (inverseMagneticPermeabilityAverageXZ);
}

/*! \brief Get const reference to averaged magneticPermeability in xz-plane
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getInverseMagneticPermeabilityAverageXZ() const
{
    SCAI_ASSERT(dirtyFlagAveraging == false, "Averaging has to be recalculated! ");
    return (inverseMagneticPermeabilityAverageXZ);
}

/*! \brief Get const reference to averaged magneticPermeability in xy-plane
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getInverseMagneticPermeabilityAverageXY()
{
    // If Averaging is outdated or has to be calculated for the first time, then recalculate averaging
    if (dirtyFlagAveraging == true) {
        calculateAveraging();
    }
    return (inverseMagneticPermeabilityAverageXY);
}

/*! \brief Get const reference to averaged magneticPermeability in xy-plane
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getInverseMagneticPermeabilityAverageXY() const
{
    SCAI_ASSERT(dirtyFlagAveraging == false, "Averaging has to be recalculated! ");
    return (inverseMagneticPermeabilityAverageXY);
}

/*! \brief Get const reference to averaged electricConductivity in x-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getElectricConductivityAverageX()
{
    // If Averaging is outdated or has to be calculated for the first time, then recalculate averaging
    if (dirtyFlagAveraging == true) {
        calculateAveraging();
    }
    return (electricConductivityAverageX);
}

/*! \brief Get const reference to averaged electricConductivity in x-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getElectricConductivityAverageX() const
{
    SCAI_ASSERT(dirtyFlagAveraging == false, "Averaging has to be recalculated! ");
    return (electricConductivityAverageX);
}

/*! \brief Get const reference to averaged electricConductivity in y-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getElectricConductivityAverageY()
{
    // If Averaging is outdated or has to be calculated for the first time, then recalculate averaging
    if (dirtyFlagAveraging == true) {
        calculateAveraging();
    }
    return (electricConductivityAverageY);
}

/*! \brief Get const reference to averaged electricConductivity in y-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getElectricConductivityAverageY() const
{
    SCAI_ASSERT(dirtyFlagAveraging == false, "Averaging has to be recalculated! ");
    return (electricConductivityAverageY);
}

/*! \brief Get const reference to averaged electricConductivity in z-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getElectricConductivityAverageZ()
{
    // If Averaging is outdated or has to be calculated for the first time, then recalculate averaging
    if (dirtyFlagAveraging == true) {
        calculateAveraging();
    }
    return (electricConductivityAverageZ);
}

/*! \brief Get const reference to averaged electricConductivity in z-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getElectricConductivityAverageZ() const
{
    SCAI_ASSERT(dirtyFlagAveraging == false, "Averaging has to be recalculated! ");
    return (electricConductivityAverageZ);
}

/*! \brief Get const reference to averaged dielectricPermittivity in x-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getDielectricPermittivityAverageX()
{
    // If Averaging is outdated or has to be calculated for the first time, then recalculate averaging
    if (dirtyFlagAveraging == true) {
        calculateAveraging();
    }
    return (dielectricPermittivityAverageX);
}

/*! \brief Get const reference to averaged dielectricPermittivity in x-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getDielectricPermittivityAverageX() const
{
    SCAI_ASSERT(dirtyFlagAveraging == false, "Averaging has to be recalculated! ");
    return (dielectricPermittivityAverageX);
}

/*! \brief Get const reference to averaged dielectricPermittivity in y-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getDielectricPermittivityAverageY()
{
    // If Averaging is outdated or has to be calculated for the first time, then recalculate averaging
    if (dirtyFlagAveraging == true) {
        calculateAveraging();
    }
    return (dielectricPermittivityAverageY);
}

/*! \brief Get const reference to averaged dielectricPermittivity in y-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getDielectricPermittivityAverageY() const
{
    SCAI_ASSERT(dirtyFlagAveraging == false, "Averaging has to be recalculated! ");
    return (dielectricPermittivityAverageY);
}

/*! \brief Get const reference to averaged dielectricPermittivity in z-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getDielectricPermittivityAverageZ()
{
    // If Averaging is outdated or has to be calculated for the first time, then recalculate averaging
    if (dirtyFlagAveraging == true) {
        calculateAveraging();
    }
    return (dielectricPermittivityAverageZ);
}

/*! \brief Get const reference to averaged dielectricPermittivity in z-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getDielectricPermittivityAverageZ() const
{
    SCAI_ASSERT(dirtyFlagAveraging == false, "Averaging has to be recalculated! ");
    return (dielectricPermittivityAverageZ);
}

/*! \brief Get const reference to averaged electricConductivityEffectiveOptical in x-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getElectricConductivityEffectiveOpticalAverageX()
{
    // If Averaging is outdated or has to be calculated for the first time, then recalculate averaging
    if (dirtyFlagAveraging == true) {
        calculateAveraging();
    }
    return (electricConductivityEffectiveOpticalAverageX);
}

/*! \brief Get const reference to averaged electricConductivityEffectiveOptical in x-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getElectricConductivityEffectiveOpticalAverageX() const
{
    SCAI_ASSERT(dirtyFlagAveraging == false, "Averaging has to be recalculated! ");
    return (electricConductivityEffectiveOpticalAverageX);
}

/*! \brief Get const reference to averaged electricConductivityEffectiveOptical in y-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getElectricConductivityEffectiveOpticalAverageY()
{
    // If Averaging is outdated or has to be calculated for the first time, then recalculate averaging
    if (dirtyFlagAveraging == true) {
        calculateAveraging();
    }
    return (electricConductivityEffectiveOpticalAverageY);
}

/*! \brief Get const reference to averaged electricConductivityEffectiveOptical in y-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getElectricConductivityEffectiveOpticalAverageY() const
{
    SCAI_ASSERT(dirtyFlagAveraging == false, "Averaging has to be recalculated! ");
    return (electricConductivityEffectiveOpticalAverageY);
}

/*! \brief Get const reference to averaged electricConductivityEffectiveOptical in z-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getElectricConductivityEffectiveOpticalAverageZ()
{
    // If Averaging is outdated or has to be calculated for the first time, then recalculate averaging
    if (dirtyFlagAveraging == true) {
        calculateAveraging();
    }
    return (electricConductivityEffectiveOpticalAverageZ);
}

/*! \brief Get const reference to averaged electricConductivityEffectiveOptical in z-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getElectricConductivityEffectiveOpticalAverageZ() const
{
    SCAI_ASSERT(dirtyFlagAveraging == false, "Averaging has to be recalculated! ");
    return (electricConductivityEffectiveOpticalAverageZ);
}

/*! \brief Get const reference to averaged dielectricPermittivityEffectiveOptical in x-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getDielectricPermittivityEffectiveOpticalAverageX()
{
    // If Averaging is outdated or has to be calculated for the first time, then recalculate averaging
    if (dirtyFlagAveraging == true) {
        calculateAveraging();
    }
    return (dielectricPermittivityEffectiveOpticalAverageX);
}

/*! \brief Get const reference to averaged dielectricPermittivityEffectiveOptical in x-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getDielectricPermittivityEffectiveOpticalAverageX() const
{
    SCAI_ASSERT(dirtyFlagAveraging == false, "Averaging has to be recalculated! ");
    return (dielectricPermittivityEffectiveOpticalAverageX);
}

/*! \brief Get const reference to averaged dielectricPermittivityEffectiveOptical in y-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getDielectricPermittivityEffectiveOpticalAverageY()
{
    // If Averaging is outdated or has to be calculated for the first time, then recalculate averaging
    if (dirtyFlagAveraging == true) {
        calculateAveraging();
    }
    return (dielectricPermittivityEffectiveOpticalAverageY);
}

/*! \brief Get const reference to averaged dielectricPermittivityEffectiveOptical in y-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getDielectricPermittivityEffectiveOpticalAverageY() const
{
    SCAI_ASSERT(dirtyFlagAveraging == false, "Averaging has to be recalculated! ");
    return (dielectricPermittivityEffectiveOpticalAverageY);
}

/*! \brief Get const reference to averaged dielectricPermittivityEffectiveOptical in z-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getDielectricPermittivityEffectiveOpticalAverageZ()
{
    // If Averaging is outdated or has to be calculated for the first time, then recalculate averaging
    if (dirtyFlagAveraging == true) {
        calculateAveraging();
    }
    return (dielectricPermittivityEffectiveOpticalAverageZ);
}

/*! \brief Get const reference to averaged dielectricPermittivityEffectiveOptical in z-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getDielectricPermittivityEffectiveOpticalAverageZ() const
{
    SCAI_ASSERT(dirtyFlagAveraging == false, "Averaging has to be recalculated! ");
    return (dielectricPermittivityEffectiveOpticalAverageZ);
}

/*! \brief Get const reference to averaged tauDielectricPermittivity in x-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getTauDielectricPermittivityAverageX()
{
    // If Averaging is outdated or has to be calculated for the first time, then recalculate averaging
    if (dirtyFlagAveraging == true) {
        calculateAveraging();
    }
    return (tauDielectricPermittivityAverageX);
}

/*! \brief Get const reference to averaged tauDielectricPermittivity in x-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getTauDielectricPermittivityAverageX() const
{
    SCAI_ASSERT(dirtyFlagAveraging == false, "Averaging has to be recalculated! ");
    return (tauDielectricPermittivityAverageX);
}

/*! \brief Get const reference to averaged tauDielectricPermittivity in y-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getTauDielectricPermittivityAverageY()
{
    // If Averaging is outdated or has to be calculated for the first time, then recalculate averaging
    if (dirtyFlagAveraging == true) {
        calculateAveraging();
    }
    return (tauDielectricPermittivityAverageY);
}

/*! \brief Get const reference to averaged tauDielectricPermittivity in y-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getTauDielectricPermittivityAverageY() const
{
    SCAI_ASSERT(dirtyFlagAveraging == false, "Averaging has to be recalculated! ");
    return (tauDielectricPermittivityAverageY);
}

/*! \brief Get const reference to averaged tauDielectricPermittivity in z-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getTauDielectricPermittivityAverageZ()
{
    // If Averaging is outdated or has to be calculated for the first time, then recalculate averaging
    if (dirtyFlagAveraging == true) {
        calculateAveraging();
    }
    return (tauDielectricPermittivityAverageZ);
}

/*! \brief Get const reference to averaged tauDielectricPermittivity in z-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getTauDielectricPermittivityAverageZ() const
{
    SCAI_ASSERT(dirtyFlagAveraging == false, "Averaging has to be recalculated! ");
    return (tauDielectricPermittivityAverageZ);
}



/* override the functions that only used by seismic wave */

/*! \brief calculate BiotCoefficient beta based on Gaussmann equation
 */
template <typename ValueType>
scai::lama::DenseVector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getBiotCoefficient()
{
    COMMON_THROWEXCEPTION("There is no BiotCoefficient in an EM modelling")
    return (BiotCoefficient);    
}

/*! \brief get BiotCoefficient beta based on Gaussmann equation
 */
template <typename ValueType>
scai::lama::DenseVector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getBiotCoefficient() const
{    
    COMMON_THROWEXCEPTION("There is no BiotCoefficient in an EM modelling")
    return (BiotCoefficient);
}

/*! \brief calculate BulkModulus Kf based on Gaussmann equation
 */
template <typename ValueType>
scai::lama::DenseVector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getBulkModulusKf()
{        
    COMMON_THROWEXCEPTION("There is no bulkModulusKf in an EM modelling")
    return (bulkModulusKf); 
}

/*! \brief get BulkModulus Kf based on Gaussmann equation
 */
template <typename ValueType>
scai::lama::DenseVector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getBulkModulusKf() const
{          
    COMMON_THROWEXCEPTION("There is no bulkModulusKf in an EM modelling")
    return (bulkModulusKf);
}

/*! \brief calculate BulkModulus M based on Gaussmann equation
 */
template <typename ValueType>
scai::lama::DenseVector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getBulkModulusM()
{        
    COMMON_THROWEXCEPTION("There is no bulkModulusM in an EM modelling")
    return (bulkModulusM);
}

/*! \brief get BulkModulus M based on Gaussmann equation
 */
template <typename ValueType>
scai::lama::DenseVector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getBulkModulusM() const
{            
    COMMON_THROWEXCEPTION("There is no bulkModulusM in an EM modelling")
    return (bulkModulusM);
}

/*! \brief Get const reference to inverseDensity model parameter
 *
 * If inverseDensity is dirty eg. because the density was modified, inverseDensity will be calculated from density.
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getInverseDensity()
{
    COMMON_THROWEXCEPTION("There is no inverseDensity in an EM modelling")
    return (inverseDensity);
}

/*! \brief Get const reference to inverseDensity model parameter
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getInverseDensity() const
{
    COMMON_THROWEXCEPTION("There is no inverseDensity in an EM modelling")
    return (inverseDensity);
}

/*! \brief Get const reference to density model parameter
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getDensity() const
{
    COMMON_THROWEXCEPTION("There is no density in an EM modelling")
    return (density);
}

/*! \brief Set density model parameter
 */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterEM<ValueType>::setDensity(scai::lama::Vector<ValueType> const &setDensity)
{
    COMMON_THROWEXCEPTION("There is no setDensity in an EM modelling")
}

/*! \brief Get const reference to first Lame model parameter
 *
 * If P-Wave modulus is dirty eg. because the P-Wave velocity was modified, P-Wave modulus will be calculated from density and P-Wave velocity.
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getPWaveModulus()
{
    COMMON_THROWEXCEPTION("There is no pWaveModulus in an EM modelling")
    return (pWaveModulus);
}

/*! \brief Get const reference to first Lame model parameter
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getPWaveModulus() const
{
    COMMON_THROWEXCEPTION("There is no pWaveModulus in an EM modelling")
    return (pWaveModulus);
}

/*! \brief Get const reference to second Lame Parameter sWaveModulus
 *
 * If S-Wave modulus is dirty eg. because the S-Wave Velocity was modified, S-Wave modulus will be calculated from density and S-Wave velocity.
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getSWaveModulus()
{
    COMMON_THROWEXCEPTION("There is no sWaveModulus in an EM modelling")

    return (sWaveModulus);
}

/*! \brief Get const reference to second Lame Parameter sWaveModulus
 *
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getSWaveModulus() const
{
    COMMON_THROWEXCEPTION("There is no sWaveModulus in an EM modelling")
    return (sWaveModulus);
}

/*! \brief Get const reference to P-wave velocity
 *
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getVelocityP() const
{
    COMMON_THROWEXCEPTION("There is no velocityP in an EM modelling")
    return (velocityP);
}

/*! \brief Set P-wave velocity model parameter
 */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterEM<ValueType>::setVelocityP(scai::lama::Vector<ValueType> const &setVelocityP)
{
    COMMON_THROWEXCEPTION("There is no setVelocityP in an EM modelling")
}

/*! \brief Get const reference to S-wave velocity
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getVelocityS() const
{
    COMMON_THROWEXCEPTION("There is no velocityS in an EM modelling")
    return (velocityS);
}

/*! \brief Set S-wave velocity model parameter
 */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterEM<ValueType>::setVelocityS(scai::lama::Vector<ValueType> const &setVelocityS)
{
    COMMON_THROWEXCEPTION("There is no setVelocityS in an EM modelling")
}

/*! \brief Get const reference to S-wave velocity
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getBulkModulusRockMatrix() const
{
    COMMON_THROWEXCEPTION("There is no bulkModulusRockMatrix in an EM modelling")
    return (bulkModulusRockMatrix);
}

/*! \brief Set S-wave velocity model parameter
 */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterEM<ValueType>::setBulkModulusRockMatrix(scai::lama::Vector<ValueType> const &setBulkModulusRockMatrix)
{
    COMMON_THROWEXCEPTION("There is no setBulkModulusRockMatrix in an EM modelling")
}

/*! \brief Get const reference to S-wave velocity
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getShearModulusRockMatrix() const
{
    COMMON_THROWEXCEPTION("There is no shearModulusRockMatrix in an EM modelling")
    return (shearModulusRockMatrix);
}

/*! \brief Set S-wave velocity model parameter
 */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterEM<ValueType>::setShearModulusRockMatrix(scai::lama::Vector<ValueType> const &setShearModulusRockMatrix)
{
    COMMON_THROWEXCEPTION("There is no setShearModulusRockMatrix in an EM modelling")
}

/*! \brief Get const reference to S-wave velocity
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getDensityRockMatrix() const
{
    COMMON_THROWEXCEPTION("There is no densityRockMatrix in an EM modelling")
    return (densityRockMatrix);
}

/*! \brief Set S-wave velocity model parameter
 */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterEM<ValueType>::setDensityRockMatrix(scai::lama::Vector<ValueType> const &setDensityRockMatrix)
{
    COMMON_THROWEXCEPTION("There is no setDensityRockMatrix in an EM modelling")
}

/*! \brief Get const reference to tauP
 *
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getTauP() const
{
    COMMON_THROWEXCEPTION("There is no tauP in an EM modelling")
    return (tauP);
}

/*! \brief Set tauP velocity model parameter
 */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterEM<ValueType>::setTauP(scai::lama::Vector<ValueType> const &setTauP)
{
    COMMON_THROWEXCEPTION("There is no setTauP in an EM modelling")
}

/*! \brief Get const reference to tauS
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getTauS() const
{
    COMMON_THROWEXCEPTION("There is no tauS in an EM modelling")
    return (tauS);
}

/*! \brief Set tauS velocity model parameter
 */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterEM<ValueType>::setTauS(scai::lama::Vector<ValueType> const &setTauS)
{
    COMMON_THROWEXCEPTION("There is no setTauS in an EM modelling")
}

/*! \brief Getter method for DensityWater */
template <typename ValueType>
ValueType const KITGPI::Modelparameter::ModelparameterEM<ValueType>::getDensityWater() const
{
    COMMON_THROWEXCEPTION("There is no DensityWater in an EM modelling")
    return (DensityWater);
}

/*! \brief Getter method for DensityAir */
template <typename ValueType>
ValueType const KITGPI::Modelparameter::ModelparameterEM<ValueType>::getDensityAir() const
{
    COMMON_THROWEXCEPTION("There is no DensityAir in an EM modelling")
    return (DensityAir);
}

/*! \brief Getter method for BulkModulusWater */
template <typename ValueType>
ValueType const KITGPI::Modelparameter::ModelparameterEM<ValueType>::getBulkModulusWater() const
{
    COMMON_THROWEXCEPTION("There is no BulkModulusWater in an EM modelling")
    ValueType BulkModulusWater = DensityWater * VelocityPWater * VelocityPWater; 
    return (BulkModulusWater);
}

/*! \brief Getter method for BulkModulusAir */
template <typename ValueType>
ValueType const KITGPI::Modelparameter::ModelparameterEM<ValueType>::getBulkModulusAir() const
{         
    COMMON_THROWEXCEPTION("There is no BulkModulusAir in an EM modelling")
    ValueType BulkModulusAir = DensityAir * VelocityPAir * VelocityPAir; 
    return (BulkModulusAir);
}

/*! \brief Getter method for CriticalPorosity */
template <typename ValueType>
ValueType const KITGPI::Modelparameter::ModelparameterEM<ValueType>::getCriticalPorosity() const
{
    COMMON_THROWEXCEPTION("There is no CriticalPorosity in an EM modelling")
    return (CriticalPorosity);
}

/*! \brief calculate averaged s-wave modulus
 *
 \param vecSWaveModulus s-wave modulus vector
 \param vecAvSWaveModulus Averaged s-wave modulus vector which is calculated
 \param avSWaveModulusMatrix Averaging matrix which is used to calculate averaged vector
 */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterEM<ValueType>::calcAveragedSWaveModulus(scai::lama::DenseVector<ValueType> &vecSWaveModulus, scai::lama::DenseVector<ValueType> &vecAvSWaveModulus, scai::lama::Matrix<ValueType> &avSWaveModulusMatrix)
{
    COMMON_THROWEXCEPTION("There is no calcAveragedSWaveModulus in an EM modelling")
}

/*! \brief Get const reference to averaged density in x-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getInverseDensityAverageX()
{
    COMMON_THROWEXCEPTION("There is no inverseDensityAverageX in an EM modelling")
    return (inverseDensityAverageX);
}

/*! \brief Get const reference to averaged density in y-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getInverseDensityAverageY()
{
    COMMON_THROWEXCEPTION("There is no inverseDensityAverageY in an EM modelling")
    return (inverseDensityAverageY);
}

/*! \brief Get const reference to averaged density in z-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getInverseDensityAverageZ()
{
    COMMON_THROWEXCEPTION("There is no inverseDensityAverageZ in an EM modelling")
    return (inverseDensityAverageZ);
}

/*! \brief Get const reference to averaged s-wave modulus in xy-plane
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getSWaveModulusAverageXY()
{
    COMMON_THROWEXCEPTION("There is no sWaveModulusAverageXY in an EM modelling")
    return (sWaveModulusAverageXY);
}

/*! \brief Get const reference to averaged s-wave modulus in xz-plane
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getSWaveModulusAverageXZ()
{
    COMMON_THROWEXCEPTION("There is no sWaveModulusAverageXZ in an EM modelling")
    return (sWaveModulusAverageXZ);
}

/*! \brief Get const reference to averaged s-wave modulus in yz-plane
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getSWaveModulusAverageYZ()
{
    COMMON_THROWEXCEPTION("There is no sWaveModulusAverageYZ in an EM modelling")
    return (sWaveModulusAverageYZ);
}

/*! \brief Get const reference to averaged tauS in xy-plane
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getTauSAverageXY()
{
    COMMON_THROWEXCEPTION("There is no tauSAverageXY in an EM modelling")
    return (tauSAverageXY);
}

/*! \brief Get const reference to averaged tauS in xz-plane
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getTauSAverageXZ()
{
    COMMON_THROWEXCEPTION("There is no tauSAverageXZ in an EM modelling")
    return (tauSAverageXZ);
}

/*! \brief Get const reference to averaged tauS in yz-plane
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getTauSAverageYZ()
{
    COMMON_THROWEXCEPTION("There is no tauSAverageYZ in an EM modelling")
    return (tauSAverageYZ);
}

/*! \brief Get const reference to averaged density in x-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getInverseDensityAverageX() const
{
    COMMON_THROWEXCEPTION("There is no inverseDensityAverageX in an EM modelling")
    return (inverseDensityAverageX);
}

/*! \brief Get const reference to averaged density in y-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getInverseDensityAverageY() const
{
    COMMON_THROWEXCEPTION("There is no inverseDensityAverageY in an EM modelling")
    return (inverseDensityAverageY);
}

/*! \brief Get const reference to averaged density in z-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getInverseDensityAverageZ() const
{
    COMMON_THROWEXCEPTION("There is no inverseDensityAverageZ in an EM modelling")
    return (inverseDensityAverageZ);
}

/*! \brief Get const reference to averaged s-wave modulus in xy-plane
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getSWaveModulusAverageXY() const
{
    COMMON_THROWEXCEPTION("There is no sWaveModulusAverageXY in an EM modelling")
    return (sWaveModulusAverageXY);
}

/*! \brief Get const reference to averaged s-wave modulus in xz-plane
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getSWaveModulusAverageXZ() const
{
    COMMON_THROWEXCEPTION("There is no sWaveModulusAverageXZ in an EM modelling")
    return (sWaveModulusAverageXZ);
}

/*! \brief Get const reference to averaged s-wave modulus in yz-plane
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getSWaveModulusAverageYZ() const
{
    COMMON_THROWEXCEPTION("There is no sWaveModulusAverageYZ in an EM modelling")
    return (sWaveModulusAverageYZ);
}

/*! \brief Get const reference to averaged tauS in xy-plane
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getTauSAverageXY() const
{
    COMMON_THROWEXCEPTION("There is no tauSAverageXY in an EM modelling")
    return (tauSAverageXY);
}

/*! \brief Get const reference to averaged tauS in xz-plane
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getTauSAverageXZ() const
{
    COMMON_THROWEXCEPTION("There is no tauSAverageXZ in an EM modelling")
    return (tauSAverageXZ);
}

/*! \brief Get const reference to averaged tauS in yz-plane
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getTauSAverageYZ() const
{
    COMMON_THROWEXCEPTION("There is no tauSAverageYZ in an EM modelling")
    return (tauSAverageYZ);
}

template class KITGPI::Modelparameter::ModelparameterEM<float>;
template class KITGPI::Modelparameter::ModelparameterEM<double>;
