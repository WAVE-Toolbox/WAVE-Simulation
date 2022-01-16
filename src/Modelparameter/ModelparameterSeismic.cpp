#include "ModelparameterSeismic.hpp"
#include "../IO/IO.hpp"

using namespace scai;
using namespace KITGPI;

/*! \brief Prepare modelparameter for inversion
 */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::prepareForInversion(Configuration::Configuration const &config, scai::dmemo::CommunicatorPtr comm)
{
    HOST_PRINT(comm, "", "Preparation of the model parameters inversion\n");
    this->setParameterisation(config.getAndCatch("parameterisation", 0));
    this->setEffectiveParameterisation(config.getAndCatch("effectiveParameterisation", 0));
    this->setInversionType(config.getAndCatch("inversionType", 1));
    this->setGradientType(config.getAndCatch("gradientType", 0));
    this->setDecomposeType(config.getAndCatch("decomposition", 0));
    HOST_PRINT(comm, "", "Model ready!\n\n");
}

/*! \brief Set relaxation frequency
 */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::setRelaxationFrequency(std::vector<ValueType> const setRelaxationFrequency)
{
    dirtyFlagSWaveModulus = true; // the modulus vector is now dirty
    dirtyFlagPWaveModulus = true; // the modulus vector is now dirty
    dirtyFlagAveraging = true;    // If S-Wave velocity will be changed, averaging needs to be redone
    relaxationFrequency = setRelaxationFrequency;
}

/*! \brief Set number of Relaxation Mechanisms
 */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::setNumRelaxationMechanisms(IndexType const setNumRelaxationMechanisms)
{
    dirtyFlagSWaveModulus = true; // the modulus vector is now dirty
    dirtyFlagPWaveModulus = true; // the modulus vector is now dirty
    dirtyFlagAveraging = true;    // If S-Wave velocity will be changed, averaging needs to be redone
    numRelaxationMechanisms = setNumRelaxationMechanisms;
}

/*! \brief calculate BiotCoefficient beta based on Gaussmann equation
 */
template <typename ValueType>
scai::lama::DenseVector<ValueType> const &KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getBiotCoefficient()
{
    BiotCoefficient = porosity / CriticalPorosity;    
    Common::searchAndReplace<ValueType>(BiotCoefficient, 0.0, 0.0, 1);
    Common::searchAndReplace<ValueType>(BiotCoefficient, 1.0, 1.0, 2);

    return (BiotCoefficient);    
}

/*! \brief get BiotCoefficient beta based on Gaussmann equation
 */
template <typename ValueType>
scai::lama::DenseVector<ValueType> const &KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getBiotCoefficient() const
{    
    return (BiotCoefficient);
}

/*! \brief calculate BulkModulus Kf based on Gaussmann equation
 */
template <typename ValueType>
scai::lama::DenseVector<ValueType> const &KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getBulkModulusKf()
{        
    scai::lama::DenseVector<ValueType> temp;
                        
    bulkModulusKf = saturation / this->getBulkModulusWater();
    temp = 1 - saturation;
    temp /= this->getBulkModulusAir();    
    bulkModulusKf += temp; 
    bulkModulusKf = 1 / bulkModulusKf;
    Common::replaceInvalid<ValueType>(bulkModulusKf, 0.0); 
    Common::searchAndReplace<ValueType>(bulkModulusKf, 0.0, 0.0, 1);
    
    return (bulkModulusKf); 
}

/*! \brief get BulkModulus Kf based on Gaussmann equation
 */
template <typename ValueType>
scai::lama::DenseVector<ValueType> const &KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getBulkModulusKf() const
{          
    return (bulkModulusKf);
}

/*! \brief calculate BulkModulus M based on Gaussmann equation
 */
template <typename ValueType>
scai::lama::DenseVector<ValueType> const &KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getBulkModulusM()
{        
    scai::lama::DenseVector<ValueType> Kf;
    scai::lama::DenseVector<ValueType> K_ma;
    scai::lama::DenseVector<ValueType> beta;    
    scai::lama::DenseVector<ValueType> temp;
    
    K_ma = this->getBulkModulusRockMatrix();                     
    beta = this->getBiotCoefficient();
    Kf = this->getBulkModulusKf();
    
    bulkModulusM = porosity / Kf;
    temp = beta - porosity;
    temp /= K_ma;
    bulkModulusM += temp;
    bulkModulusM = 1 / bulkModulusM;
    Common::replaceInvalid<ValueType>(bulkModulusM, 0.0);
    Common::searchAndReplace<ValueType>(bulkModulusM, 0.0, 0.0, 1);
    
    return (bulkModulusM);
}

/*! \brief get BulkModulus M based on Gaussmann equation
 */
template <typename ValueType>
scai::lama::DenseVector<ValueType> const &KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getBulkModulusM() const
{            
    return (bulkModulusM);
}

/*! \brief Calculate a modulus from velocity
 *
 *  Calculates Modulus = pow(Velocity,2) *  Density
 *
 \param vecVelocity Velocity-Vector which will be used in the calculation
 \param vecDensity Density-Vector which will be used in the calculation
 \param vecModulus Modulus-Vector which is calculated
 */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::calcModulusFromVelocity(scai::lama::Vector<ValueType> &vecVelocity, scai::lama::Vector<ValueType> &vecDensity, scai::lama::Vector<ValueType> &vecModulus)
{
    vecModulus = vecDensity;
    vecModulus *= vecVelocity;
    vecModulus *= vecVelocity;
};

/*! \brief Calculate velocities from a modulus
 *
 *  Calculates Velocity = sqrt( Modulu / Density )
 *
 \param vecModulus Modulus-Vector which will be used in the calculation
 \param vecDensity Density-Vector which will be used in the calculation
 \param vecVelocity Velocity-Vector which is calculated
 */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::calcVelocityFromModulus(scai::lama::Vector<ValueType> &vecModulus, scai::lama::Vector<ValueType> &vecDensity, scai::lama::Vector<ValueType> &vecVelocity)
{
    /* Modulus = pow(velocity,2) * Density */
    /* Velocity = sqrt( Modulus / Density )  */
    vecVelocity = vecModulus / vecDensity; /* = Modulus / Density */
    vecVelocity = lama::sqrt(vecVelocity);    /* = sqrt( Modulus / Density ) */
    scai::lama::DenseVector<ValueType> temp;
    temp = vecVelocity;
    Common::replaceInvalid<ValueType>(temp, 0.0); // in case of density = 0 in vacuum.
    vecVelocity = temp;
};

/*! \brief Get const reference to inverseDensity model parameter
 *
 * If inverseDensity is dirty eg. because the density was modified, inverseDensity will be calculated from density.
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getInverseDensity()
{
    if (dirtyFlagInverseDensity) {
        HOST_PRINT(density.getDistributionPtr()->getCommunicatorPtr(), "", "Inverse density will be calculated from density\n");
        inverseDensity = 1 / density;
        dirtyFlagInverseDensity = false;
    }
    return (inverseDensity);
}

/*! \brief Get const reference to inverseDensity model parameter
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getInverseDensity() const
{
    SCAI_ASSERT(dirtyFlagInverseDensity == false, "Inverse density has to be recalculated, prepareForModelling before run forward simulation! ");
    return (inverseDensity);
}

/*! \brief Get const reference to density model parameter
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getDensity() const
{
    return (density);
}

/*! \brief Set density model parameter
 */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::setDensity(scai::lama::Vector<ValueType> const &setDensity)
{
    dirtyFlagPWaveModulus = true;
    dirtyFlagSWaveModulus = true;
    dirtyFlagInverseDensity = true; // If density will be changed, the inverse has to be refreshed if it is accessed
    dirtyFlagAveraging = true;      // If inverseDensity will be changed, averaging needs to be redone
    density = setDensity;
}

/*! \brief Get const reference to first Lame model parameter
 *
 * If P-Wave modulus is dirty eg. because the P-Wave velocity was modified, P-Wave modulus will be calculated from density and P-Wave velocity.
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getPWaveModulus()
{
    // If the modulus is dirty, than recalculate
    if (dirtyFlagPWaveModulus) {
        HOST_PRINT(velocityP.getDistributionPtr()->getCommunicatorPtr(), "", "P-Wave modulus will be calculated from density and P-Wave velocity\n");
        calcModulusFromVelocity(velocityP, density, pWaveModulus);
        dirtyFlagPWaveModulus = false;
    }

    return (pWaveModulus);
}

/*! \brief Get const reference to first Lame model parameter
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getPWaveModulus() const
{
    SCAI_ASSERT(dirtyFlagPWaveModulus == false, "P-Wave Modulus has to be recalculated! ");
    return (pWaveModulus);
}

/*! \brief Get const reference to second Lame Parameter sWaveModulus
 *
 * If S-Wave modulus is dirty eg. because the S-Wave Velocity was modified, S-Wave modulus will be calculated from density and S-Wave velocity.
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getSWaveModulus()
{

    // If the modulus is dirty, than recalculate
    if (dirtyFlagSWaveModulus) {
        HOST_PRINT(velocityS.getDistributionPtr()->getCommunicatorPtr(), "", "S-Wave modulus will be calculated from density and S-Wave velocity\n");
        calcModulusFromVelocity(velocityS, density, sWaveModulus);
        dirtyFlagSWaveModulus = false;
    }

    return (sWaveModulus);
}

/*! \brief Get const reference to second Lame Parameter sWaveModulus
 *
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getSWaveModulus() const
{
    SCAI_ASSERT(dirtyFlagSWaveModulus == false, "Modulus has to be recalculated! ");
    return (sWaveModulus);
}

/*! \brief Get const reference to P-wave velocity
 *
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getVelocityP() const
{
    return (velocityP);
}

/*! \brief Set P-wave velocity model parameter
 */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::setVelocityP(scai::lama::Vector<ValueType> const &setVelocityP)
{
    dirtyFlagPWaveModulus = true; // the modulus vector is now dirty
    velocityP = setVelocityP;
}

/*! \brief Get const reference to S-wave velocity
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getVelocityS() const
{
    return (velocityS);
}

/*! \brief Set S-wave velocity model parameter
 */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::setVelocityS(scai::lama::Vector<ValueType> const &setVelocityS)
{
    dirtyFlagSWaveModulus = true; // the modulus vector is now dirty
    dirtyFlagAveraging = true;    // If S-Wave velocity will be changed, averaging needs to be redone
    velocityS = setVelocityS;
}

/*! \brief Get const reference to S-wave velocity
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getBulkModulusRockMatrix() const
{
    return (bulkModulusRockMatrix);
}

/*! \brief Set S-wave velocity model parameter
 */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::setBulkModulusRockMatrix(scai::lama::Vector<ValueType> const &setBulkModulusRockMatrix)
{
    bulkModulusRockMatrix = setBulkModulusRockMatrix;
}

/*! \brief Get const reference to S-wave velocity
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getShearModulusRockMatrix() const
{
    return (shearModulusRockMatrix);
}

/*! \brief Set S-wave velocity model parameter
 */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::setShearModulusRockMatrix(scai::lama::Vector<ValueType> const &setShearModulusRockMatrix)
{
    shearModulusRockMatrix = setShearModulusRockMatrix;
}

/*! \brief Get const reference to S-wave velocity
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getDensityRockMatrix() const
{
    return (densityRockMatrix);
}

/*! \brief Set S-wave velocity model parameter
 */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::setDensityRockMatrix(scai::lama::Vector<ValueType> const &setDensityRockMatrix)
{
    densityRockMatrix = setDensityRockMatrix;
}

/*! \brief Get const reference to tauP
 *
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getTauP() const
{
    return (tauP);
}

/*! \brief Set tauP velocity model parameter
 */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::setTauP(scai::lama::Vector<ValueType> const &setTauP)
{
    dirtyFlagPWaveModulus = true; // the modulus vector is now dirty
    tauP = setTauP;
}

/*! \brief Get const reference to tauS
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getTauS() const
{
    return (tauS);
}

/*! \brief Set tauS velocity model parameter
 */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::setTauS(scai::lama::Vector<ValueType> const &setTauS)
{
    tauS = setTauS;
    dirtyFlagSWaveModulus = true; // the modulus vector is now dirty
    dirtyFlagAveraging = true;    // If S-Wave velocity will be changed, averaging needs to be redone
}

/*! \brief Getter method for DensityWater */
template <typename ValueType>
ValueType const KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getDensityWater() const
{
    return (DensityWater);
}

/*! \brief Getter method for DensityAir */
template <typename ValueType>
ValueType const KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getDensityAir() const
{
    return (DensityAir);
}

/*! \brief Getter method for BulkModulusWater */
template <typename ValueType>
ValueType const KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getBulkModulusWater() const
{
    ValueType BulkModulusWater = DensityWater * VelocityPWater * VelocityPWater; 
    return (BulkModulusWater);
}

/*! \brief Getter method for BulkModulusAir */
template <typename ValueType>
ValueType const KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getBulkModulusAir() const
{         
    ValueType BulkModulusAir = DensityAir * VelocityPAir * VelocityPAir; 
    return (BulkModulusAir);
}

/*! \brief Getter method for CriticalPorosity */
template <typename ValueType>
ValueType const KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getCriticalPorosity() const
{
    return (CriticalPorosity);
}

/*! \brief calculate averaged s-wave modulus
 *
 \param vecSWaveModulus s-wave modulus vector
 \param vecAvSWaveModulus Averaged s-wave modulus vector which is calculated
 \param avSWaveModulusMatrix Averaging matrix which is used to calculate averaged vector
 */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::calcAveragedSWaveModulus(scai::lama::DenseVector<ValueType> &vecSWaveModulus, scai::lama::DenseVector<ValueType> &vecAvSWaveModulus, scai::lama::Matrix<ValueType> &avSWaveModulusMatrix)
{
    //replace values smaller than 1.0 with 1.0 (to avoid infs)
    Common::searchAndReplace<ValueType>(vecSWaveModulus, 1.0, 1.0, 1);

    vecAvSWaveModulus = 1 / vecSWaveModulus;
    auto temp = lama::eval<lama::DenseVector<ValueType>>(avSWaveModulusMatrix * vecAvSWaveModulus);
    vecAvSWaveModulus = 1 / temp;

    // replace values smaller than 4.0 with 0.0 (improved vacuum formulation)
    Common::searchAndReplace<ValueType>(vecAvSWaveModulus, 4.0, 0.0, 1);
}

/*! \brief Get const reference to averaged density in x-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getInverseDensityAverageX()
{
    // If Averaging is outdated or has to be calculated for the first time, than recalculate averaging
    if (dirtyFlagAveraging == true) {
        calculateAveraging();
    }
    return (inverseDensityAverageX);
}

/*! \brief Get const reference to averaged density in y-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getInverseDensityAverageY()
{
    // If Averaging is outdated or has to be calculated for the first time, than recalculate averaging
    if (dirtyFlagAveraging == true) {
        calculateAveraging();
    }
    return (inverseDensityAverageY);
}

/*! \brief Get const reference to averaged density in z-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getInverseDensityAverageZ()
{
    // If Averaging is outdated or has to be calculated for the first time, than recalculate averaging
    if (dirtyFlagAveraging == true) {
        calculateAveraging();
    }
    return (inverseDensityAverageZ);
}

/*! \brief Get const reference to averaged s-wave modulus in xy-plane
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getSWaveModulusAverageXY()
{
    // If Averaging is outdated or has to be calculated for the first time, than recalculate averaging
    if (dirtyFlagAveraging == true) {
        calculateAveraging();
    }
    return (sWaveModulusAverageXY);
}

/*! \brief Get const reference to averaged s-wave modulus in xz-plane
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getSWaveModulusAverageXZ()
{
    // If Averaging is outdated or has to be calculated for the first time, than recalculate averaging
    if (dirtyFlagAveraging == true) {
        calculateAveraging();
    }
    return (sWaveModulusAverageXZ);
}

/*! \brief Get const reference to averaged s-wave modulus in yz-plane
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getSWaveModulusAverageYZ()
{
    // If Averaging is outdated or has to be calculated for the first time, than recalculate averaging
    if (dirtyFlagAveraging == true) {
        calculateAveraging();
    }
    return (sWaveModulusAverageYZ);
}

/*! \brief Get const reference to averaged tauS in xy-plane
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getTauSAverageXY()
{
    // If Averaging is outdated or has to be calculated for the first time, than recalculate averaging
    if (dirtyFlagAveraging == true) {
        calculateAveraging();
    }
    return (tauSAverageXY);
}

/*! \brief Get const reference to averaged tauS in xz-plane
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getTauSAverageXZ()
{
    // If Averaging is outdated or has to be calculated for the first time, than recalculate averaging
    if (dirtyFlagAveraging == true) {
        calculateAveraging();
    }
    return (tauSAverageXZ);
}

/*! \brief Get const reference to averaged tauS in yz-plane
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getTauSAverageYZ()
{
    // If Averaging is outdated or has to be calculated for the first time, than recalculate averaging
    if (dirtyFlagAveraging == true) {
        calculateAveraging();
    }
    return (tauSAverageYZ);
}

/*! \brief Get const reference to averaged density in x-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getInverseDensityAverageX() const
{
    SCAI_ASSERT(dirtyFlagAveraging == false, "Averaging has to be recalculated! ");
    return (inverseDensityAverageX);
}

/*! \brief Get const reference to averaged density in y-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getInverseDensityAverageY() const
{
    SCAI_ASSERT(dirtyFlagAveraging == false, "Averaging has to be recalculated! ");
    return (inverseDensityAverageY);
}

/*! \brief Get const reference to averaged density in z-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getInverseDensityAverageZ() const
{
    SCAI_ASSERT(dirtyFlagAveraging == false, "Averaging has to be recalculated! ");
    return (inverseDensityAverageZ);
}

/*! \brief Get const reference to averaged s-wave modulus in xy-plane
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getSWaveModulusAverageXY() const
{
    SCAI_ASSERT(dirtyFlagAveraging == false, "Averaging has to be recalculated! ");
    return (sWaveModulusAverageXY);
}

/*! \brief Get const reference to averaged s-wave modulus in xz-plane
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getSWaveModulusAverageXZ() const
{
    SCAI_ASSERT(dirtyFlagAveraging == false, "Averaging has to be recalculated! ");
    return (sWaveModulusAverageXZ);
}

/*! \brief Get const reference to averaged s-wave modulus in yz-plane
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getSWaveModulusAverageYZ() const
{
    SCAI_ASSERT(dirtyFlagAveraging == false, "Averaging has to be recalculated! ");
    return (sWaveModulusAverageYZ);
}

/*! \brief Get const reference to averaged tauS in xy-plane
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getTauSAverageXY() const
{
    SCAI_ASSERT(dirtyFlagAveraging == false, "Averaging has to be recalculated! ");
    return (tauSAverageXY);
}

/*! \brief Get const reference to averaged tauS in xz-plane
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getTauSAverageXZ() const
{
    SCAI_ASSERT(dirtyFlagAveraging == false, "Averaging has to be recalculated! ");
    return (tauSAverageXZ);
}

/*! \brief Get const reference to averaged tauS in yz-plane
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getTauSAverageYZ() const
{
    SCAI_ASSERT(dirtyFlagAveraging == false, "Averaging has to be recalculated! ");
    return (tauSAverageYZ);
}




/* override the functions that only used by EM wave */

/*! \brief Set method for Archie factors a, m, n */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::setArchieFactors(ValueType const &aArchie_const, ValueType const &mArchie_const, ValueType const &nArchie_const)
{        
    COMMON_THROWEXCEPTION("There is no setArchieFactors in an Seismic modelling")
}

/*! \brief Getter method for Archie factors a */
template <typename ValueType>
ValueType const KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getArchie_a() const
{        
    COMMON_THROWEXCEPTION("There is no aArchie in an Seismic modelling")
    return (aArchie);
}

/*! \brief Getter method for Archie factors m */
template <typename ValueType>
ValueType const KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getArchie_m() const
{        
    COMMON_THROWEXCEPTION("There is no mArchie in an Seismic modelling")
    return (mArchie);
}

/*! \brief Getter method for Archie factors n */
template <typename ValueType>
ValueType const KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getArchie_n() const
{        
    COMMON_THROWEXCEPTION("There is no nArchie in an Seismic modelling")
    return (nArchie);
}

/*! \brief Getter method for DielectricPermittivityVacuum */
template <typename ValueType>
ValueType const KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getDielectricPermittivityVacuum() const
{        
    COMMON_THROWEXCEPTION("There is no DielectricPermittivityVacuum in an Seismic modelling")
    return (DielectricPermittivityVacuum);
}

/*! \brief Getter method for ElectricConductivityReference */
template <typename ValueType>
ValueType const KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getElectricConductivityReference() const
{        
    COMMON_THROWEXCEPTION("There is no ElectricConductivityReference in an Seismic modelling")
    return (ElectricConductivityReference);
}

/*! \brief Getter method for TauDielectricPermittivityReference */
template <typename ValueType>
ValueType const KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getTauDielectricPermittivityReference() const
{        
    COMMON_THROWEXCEPTION("There is no TauDielectricPermittivityReference in an Seismic modelling")
    return (TauDielectricPermittivityReference);
}

/*! \brief Getter method for TauElectricConductivityReference */
template <typename ValueType>
ValueType const KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getTauElectricConductivityReference() const
{        
    COMMON_THROWEXCEPTION("There is no getTauElectricConductivityReference in an Seismic modelling")
    return (0);
}

/*! \brief Getter method for RelativeDielectricPermittivityWater */
template <typename ValueType>
ValueType const KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getRelativeDielectricPermittivityWater() const
{        
    COMMON_THROWEXCEPTION("There is no RelativeDielectricPermittivityWater in an Seismic modelling")
    return (RelativeDielectricPermittivityWater);
}

/*! \brief Getter method for RelativeDielectricPermittivityVacuum */
template <typename ValueType>
ValueType const KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getRelativeDielectricPermittivityVacuum() const
{        
    COMMON_THROWEXCEPTION("There is no RelativeDielectricPermittivityVacuum in an Seismic modelling")
    return (RelativeDielectricPermittivityVacuum);
}

/*! \brief calculate ElectricConductivityReference from the source center frequency and DielectricPermittivityVacuum */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::calcElectricConductivityReference(ValueType const CenterFrequencyCPML)
{
    COMMON_THROWEXCEPTION("There is no calcElectricConductivityReference in an Seismic modelling") 
}

/*! \brief Get const reference to EM-wave velocity
 * 
 * If EM-Wave velocity is dirty eg. because the EM-Wave velocity was modified, EM-Wave velocity will be calculated from magneticPermeability and dielectricPermittivity.
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getVelocityEM()
{
    COMMON_THROWEXCEPTION("There is no velocivityEM in an Seismic modelling")
    return (velocivityEM);
}

/*! \brief Get const reference to EM-wave velocity
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getVelocityEM() const
{
    COMMON_THROWEXCEPTION("There is no velocivityEM in an Seismic modelling")
    return (velocivityEM);
}

/*! \brief Get const reference to magneticPermeability model parameter
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getMagneticPermeability() const
{
    COMMON_THROWEXCEPTION("There is no magneticPermeability in an Seismic modelling")
    return (magneticPermeability);
}

/*! \brief Set magneticPermeability model parameter
 */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::setMagneticPermeability(scai::lama::Vector<ValueType> const &setMagneticPermeability)
{
    COMMON_THROWEXCEPTION("There is no setMagneticPermeability in an Seismic modelling")
}

/*! \brief Get const reference to electricConductivity
 *
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getElectricConductivity() const
{
    COMMON_THROWEXCEPTION("There is no electricConductivity in an Seismic modelling")
    return (electricConductivity);
}

/*! \brief Set electricConductivity
 */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::setElectricConductivity(scai::lama::Vector<ValueType> const &setElectricConductivity)
{
    COMMON_THROWEXCEPTION("There is no setElectricConductivity in an Seismic modelling")
}

/*! \brief Get const reference to dielectricPermittivity
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getDielectricPermittivity() const
{
    COMMON_THROWEXCEPTION("There is no dielectricPermittivity in an Seismic modelling")
    return (dielectricPermittivity);
}

/*! \brief Set dielectricPermittivity
 */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::setDielectricPermittivity(scai::lama::Vector<ValueType> const &setDielectricPermittivity)
{
    COMMON_THROWEXCEPTION("There is no setDielectricPermittivity in an Seismic modelling")
}

/*! \brief Get const reference to electricConductivityRealEffective
 */
template <typename ValueType>
scai::lama::DenseVector<ValueType> const KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getElectricConductivityRealEffective() const
{
    scai::lama::DenseVector<ValueType> electricConductivityRealEffective;
    COMMON_THROWEXCEPTION("There is no electricConductivityRealEffective in an Seismic modelling")
    
    return (electricConductivityRealEffective);
}

/*! \brief Get const reference to dielectricPermittivityRealEffective
 */
template <typename ValueType>
scai::lama::DenseVector<ValueType> const KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getDielectricPermittivityRealEffective() const
{
    scai::lama::DenseVector<ValueType> dielectricPermittivityRealEffective;
    COMMON_THROWEXCEPTION("There is no dielectricPermittivityRealEffective in an Seismic modelling")
    
    return (dielectricPermittivityRealEffective);
}

/*! \brief Get const reference to electricConductivityStatic
 */
template <typename ValueType>
scai::lama::DenseVector<ValueType> const KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getElectricConductivityStatic(scai::lama::DenseVector<ValueType> const dielectricPermittivityRealEffective, scai::lama::DenseVector<ValueType> const electricConductivityRealEffective, scai::IndexType calculateType)
{
    scai::lama::DenseVector<ValueType> electricConductivityStatic;
    COMMON_THROWEXCEPTION("There is no electricConductivityStatic in an Seismic modelling")
    
    return (electricConductivityStatic);
}

/*! \brief Get const reference to dielectricPermittivityStatic
 */
template <typename ValueType>
scai::lama::DenseVector<ValueType> const KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getDielectricPermittivityStatic(scai::lama::DenseVector<ValueType> const dielectricPermittivityRealEffective, scai::lama::DenseVector<ValueType> const electricConductivityRealEffective, scai::IndexType calculateType)
{
    scai::lama::DenseVector<ValueType> dielectricPermittivityStatic;
    COMMON_THROWEXCEPTION("There is no dielectricPermittivityStatic in an Seismic modelling")
    
    return (dielectricPermittivityStatic);
}

/*! \brief Get const reference to tauElectricConductivity */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getTauElectricConductivity() const
{
    COMMON_THROWEXCEPTION("There is no tauElectricConductivity in an Seismic modelling")
    return (tauElectricConductivity);
}

/*! \brief Set tauElectricConductivity model parameter */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::setTauElectricConductivity(scai::lama::Vector<ValueType> const &setTauElectricConductivity)
{
    COMMON_THROWEXCEPTION("There is no setTauElectricConductivity in an Seismic modelling")
}

/*! \brief Get const reference to tauDielectricPermittivity */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getTauDielectricPermittivity() const
{
    COMMON_THROWEXCEPTION("There is no tauDielectricPermittivity in an Seismic modelling")
    return (tauDielectricPermittivity);
}

/*! \brief Set tauDielectricPermittivity model parameter */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::setTauDielectricPermittivity(scai::lama::Vector<ValueType> const &setTauDielectricPermittivity)
{
    COMMON_THROWEXCEPTION("There is no setTauDielectricPermittivity in an Seismic modelling")
}

/*! \brief Get const reference to electricConductivityWater
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getElectricConductivityWater() const
{
    COMMON_THROWEXCEPTION("There is no electricConductivityWater in an Seismic modelling")
    return (electricConductivityWater);
}

/*! \brief Set electricConductivityWater
 */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::setElectricConductivityWater(scai::lama::Vector<ValueType> const &setElectricConductivityWater)
{
    COMMON_THROWEXCEPTION("There is no setElectricConductivityWater in an Seismic modelling")
}

/*! \brief Get const reference to relativeDieletricPeimittivityRockMatrix
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getRelativeDieletricPeimittivityRockMatrix() const
{
    COMMON_THROWEXCEPTION("There is no relativeDieletricPeimittivityRockMatrix in an Seismic modelling")
    return (relativeDieletricPeimittivityRockMatrix);
}

/*! \brief Set relativeDieletricPeimittivityRockMatrix
 */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::setRelativeDieletricPeimittivityRockMatrix(scai::lama::Vector<ValueType> const &setRelativeDieletricPeimittivityRockMatrix)
{
    COMMON_THROWEXCEPTION("There is no setRelativeDieletricPeimittivityRockMatrix in an Seismic modelling")
}

/*! \brief Get const reference to averaged magneticPermeability in yz-plane
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getInverseMagneticPermeabilityAverageYZ()
{
    COMMON_THROWEXCEPTION("There is no inverseMagneticPermeabilityAverageYZ in an Seismic modelling")
    return (inverseMagneticPermeabilityAverageYZ);
}

/*! \brief Get const reference to averaged magneticPermeability in yz-plane
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getInverseMagneticPermeabilityAverageYZ() const
{
    COMMON_THROWEXCEPTION("There is no inverseMagneticPermeabilityAverageYZ in an Seismic modelling")
    return (inverseMagneticPermeabilityAverageYZ);
}

/*! \brief Get const reference to averaged magneticPermeability in xz-plane
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getInverseMagneticPermeabilityAverageXZ()
{
    COMMON_THROWEXCEPTION("There is no inverseMagneticPermeabilityAverageXZ in an Seismic modelling")
    return (inverseMagneticPermeabilityAverageXZ);
}

/*! \brief Get const reference to averaged magneticPermeability in xz-plane
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getInverseMagneticPermeabilityAverageXZ() const
{
    COMMON_THROWEXCEPTION("There is no inverseMagneticPermeabilityAverageXZ in an Seismic modelling")
    return (inverseMagneticPermeabilityAverageXZ);
}

/*! \brief Get const reference to averaged magneticPermeability in xy-plane
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getInverseMagneticPermeabilityAverageXY()
{
    COMMON_THROWEXCEPTION("There is no inverseMagneticPermeabilityAverageXY in an Seismic modelling")
    return (inverseMagneticPermeabilityAverageXY);
}

/*! \brief Get const reference to averaged magneticPermeability in xy-plane
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getInverseMagneticPermeabilityAverageXY() const
{
    COMMON_THROWEXCEPTION("There is no inverseMagneticPermeabilityAverageXY in an Seismic modelling")
    return (inverseMagneticPermeabilityAverageXY);
}

/*! \brief Get const reference to averaged electricConductivity in x-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getElectricConductivityAverageX()
{
    COMMON_THROWEXCEPTION("There is no electricConductivityAverageX in an Seismic modelling")
    return (electricConductivityAverageX);
}

/*! \brief Get const reference to averaged electricConductivity in x-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getElectricConductivityAverageX() const
{
    COMMON_THROWEXCEPTION("There is no electricConductivityAverageX in an Seismic modelling")
    return (electricConductivityAverageX);
}

/*! \brief Get const reference to averaged electricConductivity in y-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getElectricConductivityAverageY()
{
    COMMON_THROWEXCEPTION("There is no electricConductivityAverageY in an Seismic modelling")
    return (electricConductivityAverageY);
}

/*! \brief Get const reference to averaged electricConductivity in y-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getElectricConductivityAverageY() const
{
    COMMON_THROWEXCEPTION("There is no electricConductivityAverageY in an Seismic modelling")
    return (electricConductivityAverageY);
}

/*! \brief Get const reference to averaged electricConductivity in z-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getElectricConductivityAverageZ()
{
    COMMON_THROWEXCEPTION("There is no electricConductivityAverageZ in an Seismic modelling")
    return (electricConductivityAverageZ);
}

/*! \brief Get const reference to averaged electricConductivity in z-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getElectricConductivityAverageZ() const
{
    COMMON_THROWEXCEPTION("There is no electricConductivityAverageZ in an Seismic modelling")
    return (electricConductivityAverageZ);
}

/*! \brief Get const reference to averaged dielectricPermittivity in x-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getDielectricPermittivityAverageX()
{
    COMMON_THROWEXCEPTION("There is no dielectricPermittivityAverageX in an Seismic modelling")
    return (dielectricPermittivityAverageX);
}

/*! \brief Get const reference to averaged dielectricPermittivity in x-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getDielectricPermittivityAverageX() const
{
    COMMON_THROWEXCEPTION("There is no dielectricPermittivityAverageX in an Seismic modelling")
    return (dielectricPermittivityAverageX);
}

/*! \brief Get const reference to averaged dielectricPermittivity in y-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getDielectricPermittivityAverageY()
{
    COMMON_THROWEXCEPTION("There is no dielectricPermittivityAverageY in an Seismic modelling")
    return (dielectricPermittivityAverageY);
}

/*! \brief Get const reference to averaged dielectricPermittivity in y-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getDielectricPermittivityAverageY() const
{
    COMMON_THROWEXCEPTION("There is no dielectricPermittivityAverageY in an Seismic modelling")
    return (dielectricPermittivityAverageY);
}

/*! \brief Get const reference to averaged dielectricPermittivity in z-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getDielectricPermittivityAverageZ()
{
    COMMON_THROWEXCEPTION("There is no dielectricPermittivityAverageZ in an Seismic modelling")
    return (dielectricPermittivityAverageZ);
}

/*! \brief Get const reference to averaged dielectricPermittivity in z-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getDielectricPermittivityAverageZ() const
{
    COMMON_THROWEXCEPTION("There is no dielectricPermittivityAverageZ in an Seismic modelling")
    return (dielectricPermittivityAverageZ);
}

/*! \brief Get const reference to averaged tauElectricConductivity in x-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getTauElectricConductivityAverageX()
{
    COMMON_THROWEXCEPTION("There is no tauElectricConductivityAverageX in an Seismic modelling")
    return (tauElectricConductivityAverageX);
}

/*! \brief Get const reference to averaged tauElectricConductivity in x-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getTauElectricConductivityAverageX() const
{
    COMMON_THROWEXCEPTION("There is no tauElectricConductivityAverageX in an Seismic modelling")
    return (tauElectricConductivityAverageX);
}

/*! \brief Get const reference to averaged tauElectricConductivity in y-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getTauElectricConductivityAverageY()
{
    COMMON_THROWEXCEPTION("There is no tauElectricConductivityAverageY in an Seismic modelling")
    return (tauElectricConductivityAverageY);
}

/*! \brief Get const reference to averaged tauElectricConductivity in y-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getTauElectricConductivityAverageY() const
{
    COMMON_THROWEXCEPTION("There is no tauElectricConductivityAverageY in an Seismic modelling")
    return (tauElectricConductivityAverageY);
}

/*! \brief Get const reference to averaged tauElectricConductivity in z-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getTauElectricConductivityAverageZ()
{
    COMMON_THROWEXCEPTION("There is no tauElectricConductivityAverageZ in an Seismic modelling")
    return (tauElectricConductivityAverageZ);
}

/*! \brief Get const reference to averaged tauElectricConductivity in z-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getTauElectricConductivityAverageZ() const
{
    COMMON_THROWEXCEPTION("There is no tauElectricConductivityAverageZ in an Seismic modelling")
    return (tauElectricConductivityAverageZ);
}

/*! \brief Get const reference to averaged tauDielectricPermittivity in x-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getTauDielectricPermittivityAverageX()
{
    COMMON_THROWEXCEPTION("There is no tauDielectricPermittivityAverageX in an Seismic modelling")
    return (tauDielectricPermittivityAverageX);
}

/*! \brief Get const reference to averaged tauDielectricPermittivity in x-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getTauDielectricPermittivityAverageX() const
{
    COMMON_THROWEXCEPTION("There is no tauDielectricPermittivityAverageX in an Seismic modelling")
    return (tauDielectricPermittivityAverageX);
}

/*! \brief Get const reference to averaged tauDielectricPermittivity in y-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getTauDielectricPermittivityAverageY()
{
    COMMON_THROWEXCEPTION("There is no tauDielectricPermittivityAverageY in an Seismic modelling")
    return (tauDielectricPermittivityAverageY);
}

/*! \brief Get const reference to averaged tauDielectricPermittivity in y-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getTauDielectricPermittivityAverageY() const
{
    COMMON_THROWEXCEPTION("There is no tauDielectricPermittivityAverageY in an Seismic modelling")
    return (tauDielectricPermittivityAverageY);
}

/*! \brief Get const reference to averaged tauDielectricPermittivity in z-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getTauDielectricPermittivityAverageZ()
{
    COMMON_THROWEXCEPTION("There is no tauDielectricPermittivityAverageZ in an Seismic modelling")
    return (tauDielectricPermittivityAverageZ);
}

/*! \brief Get const reference to averaged tauDielectricPermittivity in z-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterSeismic<ValueType>::getTauDielectricPermittivityAverageZ() const
{
    COMMON_THROWEXCEPTION("There is no tauDielectricPermittivityAverageZ in an Seismic modelling")
    return (tauDielectricPermittivityAverageZ);
}

template class KITGPI::Modelparameter::ModelparameterSeismic<float>;
template class KITGPI::Modelparameter::ModelparameterSeismic<double>;
