#include "Modelparameter.hpp"
#include "../IO/IO.hpp"

using namespace scai;
using namespace KITGPI;

template <typename ValueType>
ValueType KITGPI::Modelparameter::Modelparameter<ValueType>::getMemoryUsage(scai::dmemo::DistributionPtr dist, scai::IndexType numParameter)
{
    ValueType size = getMemoryModel(dist) / 1024 / 1024 * numParameter;
    return size;
}

//! \brief calculate and return memory usage the of a single ModelParameter
/*!
 */
template <typename ValueType>
ValueType KITGPI::Modelparameter::Modelparameter<ValueType>::getMemoryModel(scai::dmemo::DistributionPtr dist)
{
    /* size of a wavefield is the size of a densevector = numGridpoints*size of Valuetype*/
    return (dist->getGlobalSize() * sizeof(ValueType));
}

/*! \brief Getter method for parameterisation */
template <typename ValueType>
IndexType KITGPI::Modelparameter::Modelparameter<ValueType>::getParameterisation() const
{
    return (parameterisation);
}

/*! \brief Set method for parameterisation */
template <typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::setParameterisation(scai::IndexType const setParameterisation)
{
    parameterisation = setParameterisation;
}

/*! \brief Getter method for inversionType */
template <typename ValueType>
IndexType KITGPI::Modelparameter::Modelparameter<ValueType>::getInversionType() const
{
    return (inversionType);
}

/*! \brief Set method for inversionType */
template <typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::setInversionType(scai::IndexType const setInversionType)
{
    inversionType = setInversionType;
}

/*! \brief Get matrix that multiplies with model matrices to get a pershot
 \param dist Distribution of the pershot
 \param distBig Distribution of the big model
 \param modelCoordinates coordinate class object of the pershot
 \param modelCoordinatesBig coordinate class object of the big model
 \param cutCoordinate coordinate where to cut the pershot
 */
template <typename ValueType>
scai::lama::CSRSparseMatrix<ValueType> KITGPI::Modelparameter::Modelparameter<ValueType>::getShrinkMatrix(scai::dmemo::DistributionPtr dist, scai::dmemo::DistributionPtr distBig, Acquisition::Coordinates<ValueType> const &modelCoordinates, Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, Acquisition::coordinate3D const cutCoordinate)
{
    SparseFormat shrinkMatrix; //!< Shrink Multiplication matrix
    shrinkMatrix.allocate(dist, distBig);
      
    hmemo::HArray<IndexType> ownedIndexes; // all (global) points owned by this process
    dist->getOwnedIndexes(ownedIndexes);
    lama::MatrixAssembly<ValueType> assembly;
    IndexType indexBig = 0;
 
    for (IndexType ownedIndex : hmemo::hostReadAccess(ownedIndexes)) { // loop over all indices
        Acquisition::coordinate3D coordinate = modelCoordinates.index2coordinate(ownedIndex); // get submodel coordinate from singleIndex
        coordinate.x += cutCoordinate.x; // offset depends on shot number
        coordinate.y += cutCoordinate.y;
        coordinate.z += cutCoordinate.z;
        indexBig = modelCoordinatesBig.coordinate2index(coordinate);
        assembly.push(ownedIndex, indexBig, 1.0);
    }
    shrinkMatrix.fillFromAssembly(assembly);
    return shrinkMatrix;
}

/*! \brief Get erase-vector that erases the old values in the big model
 \param dist Distribution of the pershot
 \param distBig Distribution of the big model
 \param modelCoordinates coordinate class object of the pershot
 \param modelCoordinatesBig coordinate class object of the big model
 \param cutCoordinate coordinate where to cut the pershot
 */
template <typename ValueType>
scai::lama::SparseVector<ValueType> KITGPI::Modelparameter::Modelparameter<ValueType>::getEraseVector(scai::dmemo::DistributionPtr dist, scai::dmemo::DistributionPtr distBig, Acquisition::Coordinates<ValueType> const &modelCoordinates, Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, Acquisition::coordinate3D const cutCoordinate, scai::IndexType boundaryWidth)
{
    scai::lama::SparseVector<ValueType> eraseVector(distBig, 1.0); //!< Shrink Multiplication matrix
      
    hmemo::HArray<IndexType> ownedIndexes; // all (global) points owned by this process
    dist->getOwnedIndexes(ownedIndexes);
    lama::VectorAssembly<ValueType> assembly;
    IndexType indexBig = 0;
 
    for (IndexType ownedIndex : hmemo::hostReadAccess(ownedIndexes)) { // loop over all indices
        Acquisition::coordinate3D coordinate = modelCoordinates.index2coordinate(ownedIndex); // get submodel coordinate from singleIndex
        coordinate.x += cutCoordinate.x; // offset depends on shot number
        coordinate.y += cutCoordinate.y;
        coordinate.z += cutCoordinate.z;
        indexBig = modelCoordinatesBig.coordinate2index(coordinate);
        assembly.push(indexBig, 0.0);
    }
    eraseVector.fillFromAssembly(assembly);
    
    // damp the boundary boarders
    for (IndexType y = 0; y < modelCoordinatesBig.getNY(); y++) {
        for (IndexType i = 0; i < boundaryWidth; i++) {
            ValueType tmp = (ValueType)1.0 - (i + 1) / (ValueType)boundaryWidth;
            eraseVector[modelCoordinatesBig.coordinate2index(cutCoordinate.x+i, y, 0)] = tmp;
            eraseVector[modelCoordinatesBig.coordinate2index(cutCoordinate.x+modelCoordinates.getNX()-1-i, y, 0)] = tmp;
        }
    }
    return eraseVector;
}

/*! \brief Prepare modelparameter for inversion
 */
template <typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::prepareForInversion(Configuration::Configuration const &config, scai::dmemo::CommunicatorPtr comm)
{
    HOST_PRINT(comm, "", "Preparation of the model parameters inversion\n");
    setParameterisation(config.get<IndexType>("parameterisation"));
    setInversionType(config.get<IndexType>("inversionType"));
    HOST_PRINT(comm, "", "Model ready!\n\n");
}

/*! \brief Getter method for relaxation frequency */
template <typename ValueType>
ValueType KITGPI::Modelparameter::Modelparameter<ValueType>::getRelaxationFrequency() const
{
    return (relaxationFrequency);
}

/*! \brief Set relaxation frequency
 */
template <typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::setRelaxationFrequency(ValueType const setRelaxationFrequency)
{
    dirtyFlagSWaveModulus = true; // the modulus vector is now dirty
    dirtyFlagPWaveModulus = true; // the modulus vector is now dirty
    dirtyFlagAveraging = true;    // If S-Wave velocity will be changed, averaging needs to be redone
    relaxationFrequency = setRelaxationFrequency;
}

/*! \brief Getter method for number of relaxation mechanisms */
template <typename ValueType>
IndexType KITGPI::Modelparameter::Modelparameter<ValueType>::getNumRelaxationMechanisms() const
{
    return (numRelaxationMechanisms);
}

/*! \brief Set number of Relaxation Mechanisms
 */
template <typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::setNumRelaxationMechanisms(IndexType const setNumRelaxationMechanisms)
{
    dirtyFlagSWaveModulus = true; // the modulus vector is now dirty
    dirtyFlagPWaveModulus = true; // the modulus vector is now dirty
    dirtyFlagAveraging = true;    // If S-Wave velocity will be changed, averaging needs to be redone
    numRelaxationMechanisms = setNumRelaxationMechanisms;
}

/*! \brief Init a single modelparameter by a constant value
 *
 \param vector Singel modelparameter which will be initialized
 \param ctx Context
 \param dist Distribution
 \param value Value which will be used to initialize the single modelparameter to a homogenoeus model
 */
template <typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::initModelparameter(scai::lama::Vector<ValueType> &vector, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, ValueType value)
{
    allocateModelparameter(vector, ctx, dist);

    vector = value;
}

/*! \brief Init a single modelparameter by reading a model from an external file
 *
 *  Reads a single model from an external vector file.
 \param vector Singel modelparameter which will be initialized
 \param ctx Context
 \param dist Distribution
 \param filename Location of external file which will be read in
 \param fileFormat Input file format
 */
template <typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::initModelparameter(scai::lama::Vector<ValueType> &vector, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, std::string filename, IndexType fileFormat)
{
    HOST_PRINT(dist->getCommunicatorPtr(), "", "initModelParameter from file " << filename << "\n")
    allocateModelparameter(vector, ctx, dist);

    IO::readVector(vector, filename, fileFormat);
}

/*! \brief Allocate a single modelparameter
 */
template <typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::allocateModelparameter(scai::lama::Vector<ValueType> &vector, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist)
{
    vector.setContextPtr(ctx);
    vector.allocate(dist);
};

/*! \brief calculate BiotCoefficient beta based on Gaussmann equation
 */
template <typename ValueType>
scai::lama::DenseVector<ValueType> const &KITGPI::Modelparameter::Modelparameter<ValueType>::getBiotCoefficient()
{
    BiotCoefficient = porosity / CriticalPorosity;    
    Common::searchAndReplace<ValueType>(BiotCoefficient, 0.0, 0.0, 1);
    Common::searchAndReplace<ValueType>(BiotCoefficient, 1.0, 1.0, 2);

    return (BiotCoefficient);    
}

/*! \brief get BiotCoefficient beta based on Gaussmann equation
 */
template <typename ValueType>
scai::lama::DenseVector<ValueType> const &KITGPI::Modelparameter::Modelparameter<ValueType>::getBiotCoefficient() const
{    
    return (BiotCoefficient);
}

/*! \brief calculate BulkModulus Kf based on Gaussmann equation
 */
template <typename ValueType>
scai::lama::DenseVector<ValueType> const &KITGPI::Modelparameter::Modelparameter<ValueType>::getBulkModulusKf()
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
scai::lama::DenseVector<ValueType> const &KITGPI::Modelparameter::Modelparameter<ValueType>::getBulkModulusKf() const
{          
    return (bulkModulusKf);
}

/*! \brief calculate BulkModulus M based on Gaussmann equation
 */
template <typename ValueType>
scai::lama::DenseVector<ValueType> const &KITGPI::Modelparameter::Modelparameter<ValueType>::getBulkModulusM()
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
scai::lama::DenseVector<ValueType> const &KITGPI::Modelparameter::Modelparameter<ValueType>::getBulkModulusM() const
{            
    return (bulkModulusM);
}

/*! \brief Calculate a modulus from velocity
 *
 *  Calculates Modulus = pow(Velocity,2) *  Density
 *
 \param vecVelocity Velocity-Vector which will be used in the calculation
 \param vecDensity Density-Vector which will be used in the calculation
 \param vectorModulus Modulus-Vector which is calculated
 */
template <typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::calcModulusFromVelocity(scai::lama::Vector<ValueType> &vecVelocity, scai::lama::Vector<ValueType> &vecDensity, scai::lama::Vector<ValueType> &vectorModulus)
{
    vectorModulus = vecDensity;
    vectorModulus *= vecVelocity;
    vectorModulus *= vecVelocity;
};

/*! \brief Calculate velocities from a modulus
 *
 *  Calculates Velocity = sqrt( Modulu / Density )
 *
 \param vectorModulus Modulus-Vector which will be used in the calculation
 \param vecDensity Density-Vector which will be used in the calculation
 \param vecVelocity Velocity-Vector which is calculated
 */
template <typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::calcVelocityFromModulus(scai::lama::DenseVector<ValueType> &vectorModulus, scai::lama::DenseVector<ValueType> &vecDensity, scai::lama::DenseVector<ValueType> &vecVelocity)
{
    /* Modulus = pow(velocity,2) * Density */
    /* Velocity = sqrt( Modulus / Density )  */
    vecVelocity = vectorModulus / vecDensity; /* = Modulus / Density */
    vecVelocity = lama::sqrt(vecVelocity);    /* = sqrt( Modulus / Density ) */
    Common::replaceInvalid<ValueType>(vecVelocity, 0.0); // in case of density = 0 in vacuum.
};

/*! \brief Get const reference to inverseDensity model parameter
 *
 * If inverseDensity is dirty eg. because the density was modified, inverseDensity will be calculated from density.
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Modelparameter<ValueType>::getInverseDensity()
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
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Modelparameter<ValueType>::getInverseDensity() const
{
    SCAI_ASSERT(dirtyFlagInverseDensity == false, "Inverse density has to be recalculated, prepareForModelling before run forward simulation! ");
    return (inverseDensity);
}

/*! \brief Get const reference to density model parameter
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Modelparameter<ValueType>::getDensity() const
{
    return (density);
}

/*! \brief Set density model parameter
 */
template <typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::setDensity(scai::lama::Vector<ValueType> const &setDensity)
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
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Modelparameter<ValueType>::getPWaveModulus()
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
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Modelparameter<ValueType>::getPWaveModulus() const
{
    SCAI_ASSERT(dirtyFlagPWaveModulus == false, "P-Wave Modulus has to be recalculated! ");
    return (pWaveModulus);
}

/*! \brief Get const reference to second Lame Parameter sWaveModulus
 *
 * If S-Wave modulus is dirty eg. because the S-Wave Velocity was modified, S-Wave modulus will be calculated from density and S-Wave velocity.
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Modelparameter<ValueType>::getSWaveModulus()
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
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Modelparameter<ValueType>::getSWaveModulus() const
{
    SCAI_ASSERT(dirtyFlagSWaveModulus == false, "Modulus has to be recalculated! ");
    return (sWaveModulus);
}

/*! \brief Get const reference to P-wave velocity
 *
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Modelparameter<ValueType>::getVelocityP() const
{
    return (velocityP);
}

/*! \brief Set P-wave velocity model parameter
 */
template <typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::setVelocityP(scai::lama::Vector<ValueType> const &setVelocityP)
{
    dirtyFlagPWaveModulus = true; // the modulus vector is now dirty
    velocityP = setVelocityP;
}

/*! \brief Get const reference to S-wave velocity
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Modelparameter<ValueType>::getVelocityS() const
{
    return (velocityS);
}

/*! \brief Set S-wave velocity model parameter
 */
template <typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::setVelocityS(scai::lama::Vector<ValueType> const &setVelocityS)
{
    dirtyFlagSWaveModulus = true; // the modulus vector is now dirty
    dirtyFlagAveraging = true;    // If S-Wave velocity will be changed, averaging needs to be redone
    velocityS = setVelocityS;
}

/*! \brief Get const reference to S-wave velocity
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Modelparameter<ValueType>::getBulkModulusRockMatrix() const
{
    return (bulkModulusRockMatrix);
}

/*! \brief Set S-wave velocity model parameter
 */
template <typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::setBulkModulusRockMatrix(scai::lama::Vector<ValueType> const &setBulkModulusRockMatrix)
{
    bulkModulusRockMatrix = setBulkModulusRockMatrix;
}

/*! \brief Get const reference to S-wave velocity
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Modelparameter<ValueType>::getShearModulusRockMatrix() const
{
    return (shearModulusRockMatrix);
}

/*! \brief Set S-wave velocity model parameter
 */
template <typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::setShearModulusRockMatrix(scai::lama::Vector<ValueType> const &setShearModulusRockMatrix)
{
    shearModulusRockMatrix = setShearModulusRockMatrix;
}

/*! \brief Get const reference to S-wave velocity
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Modelparameter<ValueType>::getDensityRockMatrix() const
{
    return (densityRockMatrix);
}

/*! \brief Set S-wave velocity model parameter
 */
template <typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::setDensityRockMatrix(scai::lama::Vector<ValueType> const &setDensityRockMatrix)
{
    densityRockMatrix = setDensityRockMatrix;
}

/*! \brief Get const reference to P-wave velocity
 *
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Modelparameter<ValueType>::getPorosity() const
{
    return (porosity);
}

/*! \brief Set P-wave velocity model parameter
 */
template <typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::setPorosity(scai::lama::Vector<ValueType> const &setPorosity)
{
    porosity = setPorosity;
}

/*! \brief Get const reference to S-wave velocity
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Modelparameter<ValueType>::getSaturation() const
{
    return (saturation);
}

/*! \brief Set S-wave velocity model parameter
 */
template <typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::setSaturation(scai::lama::Vector<ValueType> const &setSaturation)
{
    saturation = setSaturation;
}

/*! \brief Get const reference to tauP
 *
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Modelparameter<ValueType>::getTauP() const
{
    return (tauP);
}

/*! \brief Set tauP velocity model parameter
 */
template <typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::setTauP(scai::lama::Vector<ValueType> const &setTauP)
{
    dirtyFlagPWaveModulus = true; // the modulus vector is now dirty
    tauP = setTauP;
}

/*! \brief Get const reference to tauS
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Modelparameter<ValueType>::getTauS() const
{
    return (tauS);
}

/*! \brief Set tauS velocity model parameter
 */
template <typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::setTauS(scai::lama::Vector<ValueType> const &setTauS)
{
    tauS = setTauS;
    dirtyFlagSWaveModulus = true; // the modulus vector is now dirty
    dirtyFlagAveraging = true;    // If S-Wave velocity will be changed, averaging needs to be redone
}

/*! \brief Getter method for DensityWater */
template <typename ValueType>
ValueType const KITGPI::Modelparameter::Modelparameter<ValueType>::getDensityWater() const
{
    return (DensityWater);
}

/*! \brief Getter method for DensityAir */
template <typename ValueType>
ValueType const KITGPI::Modelparameter::Modelparameter<ValueType>::getDensityAir() const
{
    return (DensityAir);
}

/*! \brief Getter method for BulkModulusWater */
template <typename ValueType>
ValueType const KITGPI::Modelparameter::Modelparameter<ValueType>::getBulkModulusWater() const
{
    ValueType BulkModulusWater = DensityWater * VelocityPWater * VelocityPWater; 
    return (BulkModulusWater);
}

/*! \brief Getter method for BulkModulusAir */
template <typename ValueType>
ValueType const KITGPI::Modelparameter::Modelparameter<ValueType>::getBulkModulusAir() const
{         
    ValueType BulkModulusAir = DensityAir * VelocityPAir * VelocityPAir; 
    return (BulkModulusAir);
}

/*! \brief Getter method for CriticalPorosity */
template <typename ValueType>
ValueType const KITGPI::Modelparameter::Modelparameter<ValueType>::getCriticalPorosity() const
{
    return (CriticalPorosity);
}

//! \brief Calculate averaging matrix in x-direction
/*!
 *
 \param dist Distribution
 */
template <typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::calcAverageMatrixX(Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::DistributionPtr dist)
{
    SCAI_REGION("Modelparameter.calcAverageMatrixX")
    hmemo::HArray<IndexType> ownedIndexes; // all (global) points owned by this process
    dist->getOwnedIndexes(ownedIndexes);

    lama::MatrixAssembly<ValueType> assembly;
    assembly.reserve(ownedIndexes.size() * 2);

    for (IndexType ownedIndex : hmemo::hostReadAccess(ownedIndexes)) {

        Acquisition::coordinate3D coordinate = modelCoordinates.index2coordinate(ownedIndex);
        scai::IndexType dhFactor = modelCoordinates.getDHFactor(coordinate);
        scai::IndexType X = coordinate.x + dhFactor;

        if (X < modelCoordinates.getNX()) {
            //point 1
            assembly.push(ownedIndex, ownedIndex, 1.0 / 2.0);
            //point 2
            IndexType columnIndex = modelCoordinates.coordinate2index(X, coordinate.y, coordinate.z);
            assembly.push(ownedIndex, columnIndex, 1.0 / 2.0);
        } else {
            //point 1
            assembly.push(ownedIndex, ownedIndex, 1.0);
        }
    }

    averageMatrixX = lama::zero<SparseFormat>(dist, dist);
    averageMatrixX.fillFromAssembly(assembly);
}

//! \brief Calculate averaging matrix in y-direction
/*!
 *
 \param dist Distribution
 */
template <typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::calcAverageMatrixY(Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::DistributionPtr dist)
{
    SCAI_REGION("Modelparameter.calcAverageMatrixY")
    hmemo::HArray<IndexType> ownedIndexes; // all (global) points owned by this process
    dist->getOwnedIndexes(ownedIndexes);

    lama::MatrixAssembly<ValueType> assembly;
    assembly.reserve(ownedIndexes.size() * 2);
    IndexType layer = 0;

    for (IndexType ownedIndex : hmemo::hostReadAccess(ownedIndexes)) {

        Acquisition::coordinate3D coordinate = modelCoordinates.index2coordinate(ownedIndex);

        layer = modelCoordinates.getLayer(coordinate);

        scai::IndexType Y = coordinate.y;

        if ((modelCoordinates.locatedOnInterface(coordinate)) && (modelCoordinates.getTransition(coordinate) == 0)) {
            Y += modelCoordinates.getDHFactor(layer + 1);
        } else {
            Y += modelCoordinates.getDHFactor(layer);
        }

        if (Y < modelCoordinates.getNY()) {
            assembly.push(ownedIndex, ownedIndex, 1.0 / 2.0);
            IndexType columnIndex = modelCoordinates.coordinate2index(coordinate.x, Y, coordinate.z);
            assembly.push(ownedIndex, columnIndex, 1.0 / 2.0);
        } else {
            assembly.push(ownedIndex, ownedIndex, 1.0);
        }
    }

    averageMatrixY = lama::zero<SparseFormat>(dist, dist);
    averageMatrixY.fillFromAssembly(assembly);
}

//! \brief Calculate averaging matrix in z-direction
/*!
 *
 \param dist Distribution
 */
template <typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::calcAverageMatrixZ(Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::DistributionPtr dist)
{
    SCAI_REGION("Modelparameter.calcAverageMatrixZ")
    hmemo::HArray<IndexType> ownedIndexes; // all (global) points owned by this process
    dist->getOwnedIndexes(ownedIndexes);

    lama::MatrixAssembly<ValueType> assembly;
    assembly.reserve(ownedIndexes.size() * 2);

    for (IndexType ownedIndex : hmemo::hostReadAccess(ownedIndexes)) {

        Acquisition::coordinate3D coordinate = modelCoordinates.index2coordinate(ownedIndex);
        scai::IndexType dhFactor = modelCoordinates.getDHFactor(coordinate);
        scai::IndexType Z = coordinate.z + dhFactor;
        ;

        if (Z < modelCoordinates.getNZ()) {
            //point 1
            assembly.push(ownedIndex, ownedIndex, 1.0 / 2.0);
            //point 2
            IndexType columnIndex = modelCoordinates.coordinate2index(coordinate.x, coordinate.y, Z);
            assembly.push(ownedIndex, columnIndex, 1.0 / 2.0);
        } else {
            //point 1
            assembly.push(ownedIndex, ownedIndex, 1.0);
        }
    }

    averageMatrixZ = lama::zero<SparseFormat>(dist, dist);
    averageMatrixZ.fillFromAssembly(assembly);
}

template <typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::calc4PointAverageMatrixRow(scai::IndexType rowIndex, scai::IndexType pX[], scai::IndexType pY[], scai::IndexType pZ[], scai::lama::MatrixAssembly<ValueType> &assembly, Acquisition::Coordinates<ValueType> const &modelCoordinates)
{
    IndexType columnIndex;
    /*Points
            1      2
               av
            3      4
        */

    IndexType maxX = *std::max_element(pX, pX + 5);
    IndexType maxY = *std::max_element(pY, pY + 5);
    IndexType maxZ = *std::max_element(pZ, pZ + 5);

    if ((maxX < modelCoordinates.getNX()) && (maxY < modelCoordinates.getNY()) && (maxZ < modelCoordinates.getNZ())) {
        // Point 1 (diagonal element)
        assembly.push(rowIndex, rowIndex, 1.0 / 4.0);

        //Point 2
        columnIndex = modelCoordinates.coordinate2index(pX[2], pY[2], pZ[2]);
        assembly.push(rowIndex, columnIndex, 1.0 / 4.0);

        //Point 3
        columnIndex = modelCoordinates.coordinate2index(pX[3], pY[3], pZ[3]);
        assembly.push(rowIndex, columnIndex, 1.0 / 4.0);

        //Point 4
        columnIndex = modelCoordinates.coordinate2index(pX[4], pY[4], pZ[4]);
        assembly.push(rowIndex, columnIndex, 1.0 / 4.0);
    }

    //bottom side
    if ((maxX < modelCoordinates.getNX()) && (maxY >= modelCoordinates.getNY()) && (maxZ < modelCoordinates.getNZ())) {
        columnIndex = modelCoordinates.coordinate2index(maxX, pY[1], maxZ);
        // Point 2
        assembly.push(rowIndex, columnIndex, 1.0 / 2.0);
        // Point 1
        assembly.push(rowIndex, rowIndex, 1.0 / 2.0);
    }

    // right side
    if ((maxX >= modelCoordinates.getNX()) && (maxY < modelCoordinates.getNY()) && (maxZ < modelCoordinates.getNZ())) {
        columnIndex = modelCoordinates.coordinate2index(pX[1], maxY, maxZ);
        //Point 3
        assembly.push(rowIndex, columnIndex, 1.0 / 2.0);
        // Point 1
        assembly.push(rowIndex, rowIndex, 1.0 / 2.0);
    }

    // back side
    if ((maxX < modelCoordinates.getNX()) && (maxY < modelCoordinates.getNY()) && (maxZ >= modelCoordinates.getNZ())) {
        columnIndex = modelCoordinates.coordinate2index(maxX, maxY, pZ[1]);
        // Point 3
        assembly.push(rowIndex, columnIndex, 1.0 / 2.0);
        // Point 1
        assembly.push(rowIndex, rowIndex, 1.0 / 2.0);
    }

    //bottom-right edge

    if ((maxX >= modelCoordinates.getNX()) && (maxY >= modelCoordinates.getNY()) && (maxZ < modelCoordinates.getNZ())) {
        // Point 1
        assembly.push(rowIndex, rowIndex, 1.0);
    }

    //bottom-back edge

    if ((maxX < modelCoordinates.getNX()) && (maxY >= modelCoordinates.getNY()) && (maxZ >= modelCoordinates.getNZ())) {
        // Point 1
        assembly.push(rowIndex, rowIndex, 1.0);
    }

    //right-back edge

    if ((maxX >= modelCoordinates.getNX()) && (maxY < modelCoordinates.getNY()) && (maxZ >= modelCoordinates.getNZ())) {
        // Point 1
        assembly.push(rowIndex, rowIndex, 1.0);
    }
}

//! \brief Calculate averaging matrix in x and y-direction
/*!
 *
 \param dist Distribution
 */
template <typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::calcAverageMatrixXY(Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::DistributionPtr dist)
{
    SCAI_REGION("Modelparameter.calcAverageMatrixXY")
    hmemo::HArray<IndexType> ownedIndexes; // all (global) points owned by this process
    dist->getOwnedIndexes(ownedIndexes);

    lama::MatrixAssembly<ValueType> assembly;
    assembly.reserve(ownedIndexes.size() * 4);

    for (IndexType ownedIndex : hmemo::hostReadAccess(ownedIndexes)) {

        Acquisition::coordinate3D coordinate = modelCoordinates.index2coordinate(ownedIndex);

        IndexType dhFactor = modelCoordinates.getDHFactor(coordinate);

        IndexType pX[] = {0, coordinate.x, coordinate.x + dhFactor, coordinate.x, coordinate.x + dhFactor};
        IndexType pY[] = {0, coordinate.y, coordinate.y, coordinate.y + dhFactor, coordinate.y + dhFactor};
        IndexType pZ[] = {0, coordinate.z, coordinate.z, coordinate.z, coordinate.z};

        calc4PointAverageMatrixRow(ownedIndex, pX, pY, pZ, assembly, modelCoordinates);
    }

    averageMatrixXY = lama::zero<SparseFormat>(dist, dist);
    averageMatrixXY.fillFromAssembly(assembly);
}

//! \brief Calculate averaging matrix in x and y-direction
/*!
 *
 \param dist Distribution
 */
template <typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::calcAverageMatrixXZ(Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::DistributionPtr dist)
{
    SCAI_REGION("Modelparameter.calcAverageMatrixXZ")
    //     calcAverageMatrix(sWaveModulusAverageMatrixXZ, &Modelparameter<ValueType>::calcNumberRowElements_SWaveModulusAverageMatrixXZ, &Modelparameter<ValueType>::setRowElements_SWaveModulusAverageMatrixXZ, modelCoordinates, dist);

    hmemo::HArray<IndexType> ownedIndexes; // all (global) points owned by this process
    dist->getOwnedIndexes(ownedIndexes);

    lama::MatrixAssembly<ValueType> assembly;
    assembly.reserve(ownedIndexes.size() * 4);

    for (IndexType ownedIndex : hmemo::hostReadAccess(ownedIndexes)) {

        Acquisition::coordinate3D coordinate = modelCoordinates.index2coordinate(ownedIndex);

        IndexType dhFactor = modelCoordinates.getDHFactor(coordinate);

        IndexType pX[] = {0, coordinate.x, coordinate.x + dhFactor, coordinate.x, coordinate.x + dhFactor};
        IndexType pY[] = {0, coordinate.y, coordinate.y, coordinate.y, coordinate.y};
        IndexType pZ[] = {0, coordinate.z, coordinate.z, coordinate.z + dhFactor, coordinate.z + dhFactor};

        calc4PointAverageMatrixRow(ownedIndex, pX, pY, pZ, assembly, modelCoordinates);
    }

    averageMatrixXZ = lama::zero<SparseFormat>(dist, dist);
    averageMatrixXZ.fillFromAssembly(assembly);
}

//! \brief Calculate averaging matrix in y and z-direction
/*!
 *
 \param dist Distribution
 */
template <typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::calcAverageMatrixYZ(Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::DistributionPtr dist)
{
    SCAI_REGION("Modelparameter.calcAverageMatrixYZ")
    hmemo::HArray<IndexType> ownedIndexes; // all (global) points owned by this process
    dist->getOwnedIndexes(ownedIndexes);

    lama::MatrixAssembly<ValueType> assembly;
    assembly.reserve(ownedIndexes.size() * 4);

    for (IndexType ownedIndex : hmemo::hostReadAccess(ownedIndexes)) {

        Acquisition::coordinate3D coordinate = modelCoordinates.index2coordinate(ownedIndex);

        IndexType dhFactor = modelCoordinates.getDHFactor(coordinate);

        IndexType pX[] = {0, coordinate.x, coordinate.x, coordinate.x, coordinate.x};
        IndexType pY[] = {0, coordinate.y, coordinate.y + dhFactor, coordinate.y, coordinate.y + dhFactor};
        IndexType pZ[] = {0, coordinate.z, coordinate.z, coordinate.z + dhFactor, coordinate.z + dhFactor};

        calc4PointAverageMatrixRow(ownedIndex, pX, pY, pZ, assembly, modelCoordinates);
    }

    averageMatrixYZ = lama::zero<SparseFormat>(dist, dist);
    averageMatrixYZ.fillFromAssembly(assembly);
}

/*! \brief calculate averaged inverse density modulus
 *
 \param vecDensity Density vector.
 \param vecInverseAvDensity Averaged inverse density vector which is calculated
 \param avDensityMatrix Averaging matrix which is used to calculate averaged vector
 */
template <typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::calculateInverseAveragedDensity(scai::lama::DenseVector<ValueType> &vecDensity, scai::lama::DenseVector<ValueType> &vecInverseAvDensity, scai::lama::Matrix<ValueType> &avDensityMatrix)
{
    vecInverseAvDensity = avDensityMatrix * vecDensity;
    vecInverseAvDensity = 1 / vecInverseAvDensity;

    Common::replaceInvalid<ValueType>(vecInverseAvDensity, 0.0);
}

/*! \brief calculate averaged s-wave modulus
 *
 \param vecSWaveModulus s-wave modulus vector
 \param vecAvSWaveModulus Averaged s-wave modulus vector which is calculated
 \param avSWaveModulusMatrix Averaging matrix which is used to calculate averaged vector
 */
template <typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::calculateAveragedSWaveModulus(scai::lama::DenseVector<ValueType> &vecSWaveModulus, scai::lama::DenseVector<ValueType> &vecAvSWaveModulus, scai::lama::Matrix<ValueType> &avSWaveModulusMatrix)
{
    //replace values smaller than 1.0 with 1.0 (to avoid infs)
    Common::searchAndReplace<ValueType>(vecSWaveModulus, 1.0, 1.0, 1);

    vecAvSWaveModulus = 1 / vecSWaveModulus;
    auto temp = lama::eval<lama::DenseVector<ValueType>>(avSWaveModulusMatrix * vecAvSWaveModulus);
    vecAvSWaveModulus = 1 / temp;

    // replace values smaller than 4.0 with 0.0 (improved vacuum formulation)
    Common::searchAndReplace<ValueType>(vecAvSWaveModulus, 4.0, 0.0, 1);
}

/*! \brief calculate averaged tauS
 *
 \param vecTauS TauS vector
 \param vecAvTauS Averaged tauS vector which is calculated
 \param avTauSMatrix Averaging matrix which is used to calculate averaged vector
 */
template <typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::calculateAveragedTauS(scai::lama::Vector<ValueType> &vecTauS, scai::lama::Vector<ValueType> &vecAvTauS, scai::lama::Matrix<ValueType> &avTauSMatrix)
{
    vecAvTauS = vecTauS;
    vecAvTauS = avTauSMatrix * vecAvTauS;
}

/*! \brief Get const reference to averaged density in x-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Modelparameter<ValueType>::getInverseDensityAverageX()
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
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Modelparameter<ValueType>::getInverseDensityAverageY()
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
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Modelparameter<ValueType>::getInverseDensityAverageZ()
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
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Modelparameter<ValueType>::getSWaveModulusAverageXY()
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
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Modelparameter<ValueType>::getSWaveModulusAverageXZ()
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
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Modelparameter<ValueType>::getSWaveModulusAverageYZ()
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
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Modelparameter<ValueType>::getTauSAverageXY()
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
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Modelparameter<ValueType>::getTauSAverageXZ()
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
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Modelparameter<ValueType>::getTauSAverageYZ()
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
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Modelparameter<ValueType>::getInverseDensityAverageX() const
{
    SCAI_ASSERT(dirtyFlagAveraging == false, "Averaging has to be recalculated! ");
    return (inverseDensityAverageX);
}

/*! \brief Get const reference to averaged density in y-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Modelparameter<ValueType>::getInverseDensityAverageY() const
{
    SCAI_ASSERT(dirtyFlagAveraging == false, "Averaging has to be recalculated! ");
    return (inverseDensityAverageY);
}

/*! \brief Get const reference to averaged density in z-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Modelparameter<ValueType>::getInverseDensityAverageZ() const
{
    SCAI_ASSERT(dirtyFlagAveraging == false, "Averaging has to be recalculated! ");
    return (inverseDensityAverageZ);
}

/*! \brief Get const reference to averaged s-wave modulus in xy-plane
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Modelparameter<ValueType>::getSWaveModulusAverageXY() const
{
    SCAI_ASSERT(dirtyFlagAveraging == false, "Averaging has to be recalculated! ");
    return (sWaveModulusAverageXY);
}

/*! \brief Get const reference to averaged s-wave modulus in xz-plane
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Modelparameter<ValueType>::getSWaveModulusAverageXZ() const
{
    SCAI_ASSERT(dirtyFlagAveraging == false, "Averaging has to be recalculated! ");
    return (sWaveModulusAverageXZ);
}

/*! \brief Get const reference to averaged s-wave modulus in yz-plane
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Modelparameter<ValueType>::getSWaveModulusAverageYZ() const
{
    SCAI_ASSERT(dirtyFlagAveraging == false, "Averaging has to be recalculated! ");
    return (sWaveModulusAverageYZ);
}

/*! \brief Get const reference to averaged tauS in xy-plane
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Modelparameter<ValueType>::getTauSAverageXY() const
{
    SCAI_ASSERT(dirtyFlagAveraging == false, "Averaging has to be recalculated! ");
    return (tauSAverageXY);
}

/*! \brief Get const reference to averaged tauS in xz-plane
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Modelparameter<ValueType>::getTauSAverageXZ() const
{
    SCAI_ASSERT(dirtyFlagAveraging == false, "Averaging has to be recalculated! ");
    return (tauSAverageXZ);
}

/*! \brief Get const reference to averaged tauS in yz-plane
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Modelparameter<ValueType>::getTauSAverageYZ() const
{
    SCAI_ASSERT(dirtyFlagAveraging == false, "Averaging has to be recalculated! ");
    return (tauSAverageYZ);
}

/*! \brief Getter method for dirtyFlagPWaveModulus */
template <typename ValueType>
bool KITGPI::Modelparameter::Modelparameter<ValueType>::getDirtyFlagPWaveModulus() const
{
    return (dirtyFlagPWaveModulus);
}

/*! \brief Getter method for dirtyFlagSWaveModulus */
template <typename ValueType>
bool KITGPI::Modelparameter::Modelparameter<ValueType>::getDirtyFlagSWaveModulus() const
{
    return (dirtyFlagSWaveModulus);
}

/*! \brief Getter method for dirtyFlagInverseDensity */
template <typename ValueType>
bool KITGPI::Modelparameter::Modelparameter<ValueType>::getDirtyFlagInverseDensity() const
{
    return (dirtyFlagInverseDensity);
}

/*! \brief Getter method for dirtyFlagAveraging */
template <typename ValueType>
bool KITGPI::Modelparameter::Modelparameter<ValueType>::getDirtyFlagAveraging() const
{
    return (dirtyFlagAveraging);
}

/*! \brief Overloading = Operation
 *
 \param rhs Model which is copied.
 */
template <typename ValueType>
KITGPI::Modelparameter::Modelparameter<ValueType> &KITGPI::Modelparameter::Modelparameter<ValueType>::operator=(KITGPI::Modelparameter::Modelparameter<ValueType> const &rhs)
{
    assign(rhs);
    return *this;
}

/*! \brief Overloading += Operation
 *
 \param rhs Model which is added.
 */
template <typename ValueType>
KITGPI::Modelparameter::Modelparameter<ValueType> &KITGPI::Modelparameter::Modelparameter<ValueType>::operator+=(KITGPI::Modelparameter::Modelparameter<ValueType> const &rhs)
{
    plusAssign(rhs);
    return *this;
}

/*! \brief Overloading -= Operation
 *
 \param rhs Model which is substracted.
 */
template <typename ValueType>
KITGPI::Modelparameter::Modelparameter<ValueType> &KITGPI::Modelparameter::Modelparameter<ValueType>::operator-=(KITGPI::Modelparameter::Modelparameter<ValueType> const &rhs)
{
    minusAssign(rhs);
    return *this;
}

template class KITGPI::Modelparameter::Modelparameter<float>;
template class KITGPI::Modelparameter::Modelparameter<double>;
