#include "Modelparameter.hpp"
#include "../IO/IO.hpp"

using namespace scai;
using namespace KITGPI;

template <typename ValueType>
ValueType KITGPI::Modelparameter::ModelparameterEM<ValueType>::getMemoryUsage(scai::dmemo::DistributionPtr dist, scai::IndexType numParameter)
{
    ValueType size = getMemoryModel(dist) / 1024 / 1024 * numParameter;
    return size;
}

//! \brief calculate and return memory usage the of a single ModelParameter
/*!
 */
template <typename ValueType>
ValueType KITGPI::Modelparameter::ModelparameterEM<ValueType>::getMemoryModel(scai::dmemo::DistributionPtr dist)
{
    /* size of a wavefield is the size of a densevector = numGridpoints*size of Valuetype*/
    return (dist->getGlobalSize() * sizeof(ValueType));
}

/*! \brief Getter method for parameterisation */
template <typename ValueType>
IndexType KITGPI::Modelparameter::ModelparameterEM<ValueType>::getParameterisation() const
{
    return (parameterisation);
}

/*! \brief Set method for parameterisation */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterEM<ValueType>::setParameterisation(scai::IndexType const setParameterisation)
{
    parameterisation = setParameterisation;
}

/*! \brief Getter method for inversionType */
template <typename ValueType>
IndexType KITGPI::Modelparameter::ModelparameterEM<ValueType>::getInversionType() const
{
    return (inversionType);
}

/*! \brief Set method for inversionType */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterEM<ValueType>::setInversionType(scai::IndexType const setInversionType)
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
scai::lama::CSRSparseMatrix<ValueType> KITGPI::Modelparameter::ModelparameterEM<ValueType>::getShrinkMatrix(scai::dmemo::DistributionPtr dist, scai::dmemo::DistributionPtr distBig, Acquisition::Coordinates<ValueType> const &modelCoordinates, Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, Acquisition::coordinate3D const cutCoordinate)
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
scai::lama::SparseVector<ValueType> KITGPI::Modelparameter::ModelparameterEM<ValueType>::getEraseVector(scai::dmemo::DistributionPtr dist, scai::dmemo::DistributionPtr distBig, Acquisition::Coordinates<ValueType> const &modelCoordinates, Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, Acquisition::coordinate3D const cutCoordinate, scai::IndexType boundaryWidth)
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
void KITGPI::Modelparameter::ModelparameterEM<ValueType>::prepareForInversion(Configuration::Configuration const &config, scai::dmemo::CommunicatorPtr comm)
{
    HOST_PRINT(comm, "", "Preparation of the model parameters inversion\n");
    setParameterisation(config.get<IndexType>("parameterisation"));
    setInversionType(config.get<IndexType>("inversionType"));
    calcConductivityReference(config.get<ValueType>("CenterFrequencyCPML"));
    setArchieFactors(config.get<ValueType>("aArchie"), config.get<ValueType>("mArchie"), config.get<ValueType>("nArchie"));
    HOST_PRINT(comm, "", "Model ready!\n\n");
}

/*! \brief Getter method for relaxation frequency */
template <typename ValueType>
ValueType KITGPI::Modelparameter::ModelparameterEM<ValueType>::getRelaxationFrequency() const
{
    return (relaxationFrequency);
}

/*! \brief Set relaxation frequency
 */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterEM<ValueType>::setRelaxationFrequency(ValueType const setRelaxationFrequency)
{
    dirtyFlagAveraging = true;    // If modulus will be changed, averaging needs to be redone
    dirtyFlagVelocivityEM = true; // the velocity vector is now dirty
    dirtyFlagDielectricPermittivityEMoptical = true; // the visco Modulus vector is now dirty
    dirtyFlagConductivityEMoptical = true; // the visco Modulus vector is now dirty
    relaxationFrequency = setRelaxationFrequency;
}

/*! \brief Getter method for number of relaxation mechanisms */
template <typename ValueType>
IndexType KITGPI::Modelparameter::ModelparameterEM<ValueType>::getNumRelaxationMechanisms() const
{
    return (numRelaxationMechanisms);
}

/*! \brief Set number of Relaxation Mechanisms
 */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterEM<ValueType>::setNumRelaxationMechanisms(IndexType const setNumRelaxationMechanisms)
{
    dirtyFlagAveraging = true;    // If modulus will be changed, averaging needs to be redone
    dirtyFlagVelocivityEM = true; // the velocity vector is now dirty
    dirtyFlagDielectricPermittivityEMoptical = true; // the visco Modulus vector is now dirty
    dirtyFlagConductivityEMoptical = true; // the visco Modulus vector is now dirty
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

/*! \brief Getter method for ConductivityReference */
template <typename ValueType>
ValueType const KITGPI::Modelparameter::ModelparameterEM<ValueType>::getConductivityReference() const
{        
    return (ConductivityReference);
}

/*! \brief Getter method for TauDielectricPermittivityReference */
template <typename ValueType>
ValueType const KITGPI::Modelparameter::ModelparameterEM<ValueType>::getTauDielectricPermittivityReference() const
{        
    return (TauDielectricPermittivityReference);
}

/*! \brief Getter method for TauConductivityReference */
template <typename ValueType>
ValueType const KITGPI::Modelparameter::ModelparameterEM<ValueType>::getTauConductivityReference() const
{        
    ValueType TauConductivityReference = 10 * TauDielectricPermittivityReference * this->getTauDisplacementEM(); 
    return (TauConductivityReference);
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
void KITGPI::Modelparameter::ModelparameterEM<ValueType>::writeRockMatrixParameter(std::string filename, scai::IndexType fileFormat) const
{
    IO::writeVector(conductivityEMWater, filename + ".sigmaEMw", fileFormat);
    IO::writeVector(relativeDieletricPeimittivityRockMatrix, filename + ".epsilonEMrma", fileFormat);
};

/*! \brief calculate conductivityEMWater and relativeDielectricPermittivityRockMatrix from porosity and saturation */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterEM<ValueType>::calcRockMatrixParameter(Configuration::Configuration const &config)
{            
    scai::lama::DenseVector<ValueType> temp1;
    scai::lama::DenseVector<ValueType> temp2;
    scai::lama::DenseVector<ValueType> conductivityEMWatertemp;
    scai::lama::DenseVector<ValueType> relativeDielectricPermittivityRockMatrixtemp;  
    
    conductivityEMWatertemp = this->getConductivityEM();  
    relativeDielectricPermittivityRockMatrixtemp = this->getDielectricPermittivityEM();      
    relativeDielectricPermittivityRockMatrixtemp /= DielectricPermittivityVacuum;    
    relativeDielectricPermittivityRockMatrixtemp = scai::lama::sqrt(relativeDielectricPermittivityRockMatrixtemp); 
    
    // Based on Archie equation
    temp1 = scai::lama::pow(porosity, -mArchie);
    temp2 = scai::lama::pow(saturation, -nArchie);
    conductivityEMWatertemp *= aArchie;
    conductivityEMWatertemp *= temp1;
    conductivityEMWatertemp *= temp2;
    
    Common::replaceInvalid<ValueType>(conductivityEMWatertemp, 0.0); // in case of air
    Common::searchAndReplace<ValueType>(conductivityEMWatertemp, 0, 0, 1);
    
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
    
    this->setConductivityEMWater(conductivityEMWatertemp);
    this->setRelativeDieletricPeimittivityRockMatrix(relativeDielectricPermittivityRockMatrixtemp);
}

/*! \brief calculate ConductivityReference from the source center frequency and DielectricPermittivityVacuum */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterEM<ValueType>::calcConductivityReference(ValueType const CenterFrequencyCPML)
{
    // please note here ConductivityReference is calculated by f_min=0.5*f_c, which is corresponding to the scaling factor = 0.5 in Lavoue et al (2014).
    ConductivityReference = 0.5 * 2 * M_PI * CenterFrequencyCPML * DielectricPermittivityVacuum; 
}

/*! \brief transform porosity and saturation to EM parameters */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterEM<ValueType>::calcWaveModulusFromPetrophysics()
{
    scai::lama::DenseVector<ValueType> temp1;
    scai::lama::DenseVector<ValueType> temp2;
    scai::lama::DenseVector<ValueType> conductivitytemp;
    scai::lama::DenseVector<ValueType> dielectricPermittivityEMtemp;  
    
    conductivitytemp = this->getConductivityEMWater();  
    dielectricPermittivityEMtemp = this->getRelativeDieletricPeimittivityRockMatrix();   
    
    // Based on Archie equation
    temp1 = scai::lama::pow(porosity, mArchie);
    temp2 = scai::lama::pow(saturation, nArchie);
    conductivitytemp /= aArchie;
    conductivitytemp *= temp1;
    conductivitytemp *= temp2;
    
    // Based on complex refractive index model (CRIM)    
    temp1 = saturation;
    temp1 *= sqrt(RelativeDielectricPermittivityWater);
    temp2 = 1 - saturation;
    temp2 *= sqrt(RelativeDielectricPermittivityVacuum);
    temp1 += temp2;
    temp1 *= porosity;  
    dielectricPermittivityEMtemp = scai::lama::sqrt(dielectricPermittivityEMtemp);
    temp2 = 1 - porosity;
    dielectricPermittivityEMtemp *= temp2;
    dielectricPermittivityEMtemp += temp1;    
    dielectricPermittivityEMtemp = scai::lama::pow(dielectricPermittivityEMtemp, 2);
    
    Common::searchAndReplace<ValueType>(conductivitytemp, 0, 0, 1);
    Common::searchAndReplace<ValueType>(dielectricPermittivityEMtemp, RelativeDielectricPermittivityVacuum, RelativeDielectricPermittivityVacuum, 1);
    Common::searchAndReplace<ValueType>(dielectricPermittivityEMtemp, RelativeDielectricPermittivityWater, RelativeDielectricPermittivityWater, 2);
    
    dielectricPermittivityEMtemp *= DielectricPermittivityVacuum;    
    
    this->setConductivityEM(conductivitytemp);
    this->setDielectricPermittivityEM(dielectricPermittivityEMtemp);
}

/*! \brief transform EM parameters to porosity and saturation  */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterEM<ValueType>::calcPetrophysicsFromWaveModulus()
{
    scai::lama::DenseVector<ValueType> temp1;
    scai::lama::DenseVector<ValueType> temp2;
    scai::lama::DenseVector<ValueType> porositytemp;
    scai::lama::DenseVector<ValueType> dielectricPermittivityEMtemp;  
    scai::lama::DenseVector<ValueType> relativeDieletricPeimittivityRockMatrix;  
    
    relativeDieletricPeimittivityRockMatrix = this->getRelativeDieletricPeimittivityRockMatrix();  
    dielectricPermittivityEMtemp = this->getDielectricPermittivityEM();   
    dielectricPermittivityEMtemp /= DielectricPermittivityVacuum;
    porositytemp = this->getPorosity();   
        
    // Based on complex refractive index model (CRIM)        
    temp1 = 1 - porositytemp;     
    temp2 = scai::lama::sqrt(relativeDieletricPeimittivityRockMatrix);
    temp2 *= temp1;
    temp1 = scai::lama::sqrt(dielectricPermittivityEMtemp);
    temp1 -= temp2;
    temp1 /= porositytemp;
    temp1 -= sqrt(RelativeDielectricPermittivityVacuum);
    temp1 /= (sqrt(RelativeDielectricPermittivityWater) - sqrt(RelativeDielectricPermittivityVacuum));
    Common::replaceInvalid<ValueType>(temp1, 0); // in case of air
    
    this->setSaturation(temp1);      
}
/*! \brief Init a single modelparameter by a constant value
 *
 \param vector Singel modelparameter which will be initialized
 \param ctx Context
 \param dist Distribution
 \param value Value which will be used to initialize the single modelparameter to a homogenoeus model
 */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterEM<ValueType>::initModelparameter(scai::lama::Vector<ValueType> &vector, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, ValueType value)
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
void KITGPI::Modelparameter::ModelparameterEM<ValueType>::initModelparameter(scai::lama::Vector<ValueType> &vector, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, std::string filename, IndexType fileFormat)
{
    HOST_PRINT(dist->getCommunicatorPtr(), "", "initModelParameter from file " << filename << "\n")
    allocateModelparameter(vector, ctx, dist);

    IO::readVector(vector, filename, fileFormat);
}

/*! \brief Allocate a single modelparameter
 */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterEM<ValueType>::allocateModelparameter(scai::lama::Vector<ValueType> &vector, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist)
{
    vector.setContextPtr(ctx);
    vector.allocate(dist);
};

/*! \brief Calculate a modulus from velocity
 *
 *  Calculates Modulus = 1 / (pow(Velocity,2) *  MagneticPermeabilityEM)
 *
 \param vecVelocity Velocity-Vector which will be used in the calculation
 \param vecMagneticPermeabilityEM MagneticPermeabilityEM-Vector which will be used in the calculation
 \param vecDielectricPermittivityEM Modulus-Vector which is calculated
 */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterEM<ValueType>::calcModulusFromVelocity(scai::lama::Vector<ValueType> &vecVelocity, scai::lama::Vector<ValueType> &vecMagneticPermeabilityEM, scai::lama::Vector<ValueType> &vecDielectricPermittivityEM)
{
    vecDielectricPermittivityEM = 1 / vecMagneticPermeabilityEM;
    vecDielectricPermittivityEM /= vecVelocity;
    vecDielectricPermittivityEM /= vecVelocity;
};

/*! \brief Calculate velocities from a modulus
 *
 *  Calculates Velocity = 1 / sqrt( Modulu * MagneticPermeabilityEM )
 *
 \param vecDielectricPermittivityEM Modulus-Vector which will be used in the calculation
 \param vecMagneticPermeabilityEM MagneticPermeabilityEM-Vector which will be used in the calculation
 \param vecVelocity Velocity-Vector which is calculated
 */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterEM<ValueType>::calcVelocityFromModulus(scai::lama::Vector<ValueType> &vecDielectricPermittivityEM, scai::lama::Vector<ValueType> &vecMagneticPermeabilityEM, scai::lama::Vector<ValueType> &vecVelocity)
{
    vecVelocity = vecDielectricPermittivityEM * vecMagneticPermeabilityEM; 
    vecVelocity = lama::sqrt(vecVelocity);    
    vecVelocity = 1 / vecVelocity;
};

/*! \brief Get const reference to EM-wave velocity
 * 
 * If EM-Wave velocity is dirty eg. because the EM-Wave velocity was modified, EM-Wave velocity will be calculated from magneticPermeabilityEM and dielectricPermittivityEM.
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getVelocityEM()
{
    // If the modulus is dirty, then recalculate
    if (dirtyFlagVelocivityEM) {
        HOST_PRINT(dielectricPermittivityEM.getDistributionPtr()->getCommunicatorPtr(), "", "EM-Wave velocity will be calculated from magneticPermeabilityEM and dielectricPermittivityEM\n");
        calcVelocityFromModulus(dielectricPermittivityEM, magneticPermeabilityEM, velocivityEM);
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

/*! \brief Get const reference to magneticPermeabilityEM model parameter
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getMagneticPermeabilityEM() const
{
    return (magneticPermeabilityEM);
}

/*! \brief Set magneticPermeabilityEM model parameter
 */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterEM<ValueType>::setMagneticPermeabilityEM(scai::lama::Vector<ValueType> const &setMagneticPermeabilityEM)
{
    dirtyFlagAveraging = true;      // If magneticPermeabilityEM will be changed, averaging needs to be redone
    dirtyFlagVelocivityEM = true;
    magneticPermeabilityEM = setMagneticPermeabilityEM;
}

/*! \brief Get const reference to conductivityEM
 *
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getConductivityEM() const
{
    return (conductivityEM);
}

/*! \brief Set conductivityEM
 */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterEM<ValueType>::setConductivityEM(scai::lama::Vector<ValueType> const &setConductivityEM)
{
    dirtyFlagAveraging = true;    // If modulus will be changed, averaging needs to be redone
    dirtyFlagVelocivityEM = true;   // the velocity vector is now dirty
    dirtyFlagConductivityEMoptical = true; // the visco Modulus vector is now dirty
    dirtyFlagDielectricPermittivityEMoptical = true; // the visco Modulus vector is now dirty
    conductivityEM = setConductivityEM;
}

/*! \brief Get const reference to dielectricPermittivityEM
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getDielectricPermittivityEM() const
{
    return (dielectricPermittivityEM);
}

/*! \brief Set dielectricPermittivityEM
 */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterEM<ValueType>::setDielectricPermittivityEM(scai::lama::Vector<ValueType> const &setDielectricPermittivityEM)
{
    dirtyFlagAveraging = true;    // If modulus will be changed, averaging needs to be redone
    dirtyFlagVelocivityEM = true; // the velocity vector is now dirty
    dirtyFlagConductivityEMoptical = true; // the visco Modulus vector is now dirty
    dirtyFlagDielectricPermittivityEMoptical = true; // the visco Modulus vector is now dirty
    dielectricPermittivityEM = setDielectricPermittivityEM;
}

/*! \brief Get const reference to conductivityEMoptical
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getConductivityEMoptical()
{
    // If the modulus is dirty, then recalculate
    if (dirtyFlagConductivityEMoptical) {
        HOST_PRINT(conductivityEM.getDistributionPtr()->getCommunicatorPtr(), "conductivityEMoptical will be calculated from\n dielectricPermittivityEM, conductivityEM, tauDielectricPermittivityEM\n, tauDisplacementEM and numRelaxationMechanisms \n");              
        conductivityEMoptical = tauDielectricPermittivityEM / tauDisplacementEM;    
        conductivityEMoptical *= numRelaxationMechanisms; // numRelaxationMechanisms = 1
        conductivityEMoptical /= numRelaxationMechanisms;
        conductivityEMoptical *= dielectricPermittivityEM;
        conductivityEMoptical += conductivityEM;
        dirtyFlagConductivityEMoptical = false;
    }
    
    return (conductivityEMoptical);
}
/*! \brief Get const reference to conductivityEMoptical
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getConductivityEMoptical() const
{
    SCAI_ASSERT(dirtyFlagConductivityEMoptical == false, "conductivityEMoptical has to be recalculated! ");
    
    return (conductivityEMoptical);
}

/*! \brief Set conductivityEMoptical model parameter
 */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterEM<ValueType>::setConductivityEMoptical(scai::lama::Vector<ValueType> const &setConductivityEMoptical)
{
    dirtyFlagAveraging = true;    // If modulus will be changed, averaging needs to be redone
    conductivityEMoptical = setConductivityEMoptical;
}

/*! \brief Get const reference to dielectricPermittivityEMoptical
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getDielectricPermittivityEMoptical()
{
    // If the modulus is dirty, then recalculate
    if (dirtyFlagDielectricPermittivityEMoptical) {
        HOST_PRINT(dielectricPermittivityEM.getDistributionPtr()->getCommunicatorPtr(), "dielectricPermittivityEMoptical will be calculated from\n dielectricPermittivityEM, conductivityEM, tauDielectricPermittivityEM\n, tauConductivityEM and numRelaxationMechanisms \n");
        scai::lama::DenseVector<ValueType> temp;
        temp = numRelaxationMechanisms * tauDielectricPermittivityEM;  // numRelaxationMechanisms = 1
        temp /= numRelaxationMechanisms;
        dielectricPermittivityEMoptical = 1 - temp;      
        dielectricPermittivityEMoptical *= dielectricPermittivityEM;
        temp = conductivityEM * tauConductivityEM;
        dielectricPermittivityEMoptical += temp;
        
        // optical dielectric permittivity should be larger than the dielectric permittivity of vacuum
        Common::searchAndReplace<ValueType>(dielectricPermittivityEMoptical, DielectricPermittivityVacuum, DielectricPermittivityVacuum, 1); 
        
        dirtyFlagDielectricPermittivityEMoptical = false;
    }
    
    return (dielectricPermittivityEMoptical);
}

/*! \brief Get const reference to dielectricPermittivityEMoptical
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getDielectricPermittivityEMoptical() const
{
    SCAI_ASSERT(dirtyFlagDielectricPermittivityEMoptical == false, "dielectricPermittivityEMoptical has to be recalculated! ");
    return (dielectricPermittivityEMoptical);
}

/*! \brief Set dielectricPermittivityEMoptical model parameter
 */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterEM<ValueType>::setDielectricPermittivityEMoptical(scai::lama::Vector<ValueType> const &setDielectricPermittivityEMoptical)
{
    dirtyFlagAveraging = true;    // If modulus will be changed, averaging needs to be redone
    dielectricPermittivityEMoptical = setDielectricPermittivityEMoptical;
}

/*! \brief Get const reference to tauConductivityEM */
template <typename ValueType>
ValueType const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getTauDisplacementEM()
{
    tauDisplacementEM = 1.0 / (2.0 * M_PI * this->getRelaxationFrequency()); 
    return (tauDisplacementEM);
}

/*! \brief Get const reference to tauConductivityEM */
template <typename ValueType>
ValueType const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getTauDisplacementEM() const
{
    return (tauDisplacementEM);
}

/*! \brief Get const reference to tauConductivityEM */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getTauConductivityEM() const
{
    return (tauConductivityEM);
}

/*! \brief Set tauConductivityEM model parameter */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterEM<ValueType>::setTauConductivityEM(scai::lama::Vector<ValueType> const &setTauConductivityEM)
{
    dirtyFlagAveraging = true;    // If modulus will be changed, averaging needs to be redone
    dirtyFlagVelocivityEM = true; // the velocity vector is now dirty
    dirtyFlagDielectricPermittivityEMoptical = true; // the visco Modulus vector is now dirty
    tauConductivityEM = setTauConductivityEM;
}

/*! \brief Get const reference to tauDielectricPermittivityEM */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getTauDielectricPermittivityEM() const
{
    return (tauDielectricPermittivityEM);
}

/*! \brief Set tauDielectricPermittivityEM model parameter */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterEM<ValueType>::setTauDielectricPermittivityEM(scai::lama::Vector<ValueType> const &setTauDielectricPermittivityEM)
{
    dirtyFlagAveraging = true;    // If modulus will be changed, averaging needs to be redone
    dirtyFlagVelocivityEM = true; // the velocity vector is now dirty
    dirtyFlagConductivityEMoptical = true; // the visco Modulus vector is now dirty
    dirtyFlagDielectricPermittivityEMoptical = true; // the visco Modulus vector is now dirty
    tauDielectricPermittivityEM = setTauDielectricPermittivityEM;
}

/*! \brief Get const reference to porosity */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getPorosity() const
{
    return (porosity);
}

/*! \brief Set porosity */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterEM<ValueType>::setPorosity(scai::lama::Vector<ValueType> const &setPorosity)
{
    porosity = setPorosity;
}

/*! \brief Get const reference to saturation */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getSaturation() const
{
    return (saturation);
}

/*! \brief Set saturation
 */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterEM<ValueType>::setSaturation(scai::lama::Vector<ValueType> const &setSaturation)
{
    saturation = setSaturation;
}

/*! \brief Get const reference to conductivityEMWater
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getConductivityEMWater() const
{
    return (conductivityEMWater);
}

/*! \brief Set conductivityEMWater
 */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterEM<ValueType>::setConductivityEMWater(scai::lama::Vector<ValueType> const &setConductivityEMWater)
{
    conductivityEMWater = setConductivityEMWater;
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

//! \brief Calculate averaging matrix in x-direction
/*!
 *
 \param dist Distribution
 */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterEM<ValueType>::calcAverageMatrixX(Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::DistributionPtr dist)
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
void KITGPI::Modelparameter::ModelparameterEM<ValueType>::calcAverageMatrixY(Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::DistributionPtr dist)
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
void KITGPI::Modelparameter::ModelparameterEM<ValueType>::calcAverageMatrixZ(Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::DistributionPtr dist)
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
void KITGPI::Modelparameter::ModelparameterEM<ValueType>::calc4PointAverageMatrixRow(scai::IndexType rowIndex, scai::IndexType pX[], scai::IndexType pY[], scai::IndexType pZ[], scai::lama::MatrixAssembly<ValueType> &assembly, Acquisition::Coordinates<ValueType> const &modelCoordinates)
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
void KITGPI::Modelparameter::ModelparameterEM<ValueType>::calcAverageMatrixXY(Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::DistributionPtr dist)
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
void KITGPI::Modelparameter::ModelparameterEM<ValueType>::calcAverageMatrixXZ(Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::DistributionPtr dist)
{
    SCAI_REGION("Modelparameter.calcAverageMatrixXZ")
    //     calcAverageMatrix(dielectricPermittivityEMAverageMatrixXZ, &ModelparameterEM<ValueType>::calcNumberRowElements_InverseDielectricPermittivityEMAverageMatrixXZ, &ModelparameterEM<ValueType>::setRowElements_InverseDielectricPermittivityEMAverageMatrixXZ, modelCoordinates, dist);

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
void KITGPI::Modelparameter::ModelparameterEM<ValueType>::calcAverageMatrixYZ(Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::DistributionPtr dist)
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

/*! \brief calculate averaged EMparameter
 *
 \param vecEMparameter EM parameter vector
 \param vecAveragedParameter Averaged EM parameter vector which is calculated
 \param averagedMatrix Averaging matrix which is used to calculate averaged vector
 */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterEM<ValueType>::calculateAveragedEMparameter(scai::lama::Vector<ValueType> const &vecEMparameter, scai::lama::DenseVector<ValueType> &vecAveragedParameter, scai::lama::Matrix<ValueType> const &averagedMatrix)
{
    vecAveragedParameter = averagedMatrix * vecEMparameter;
    
    Common::replaceInvalid<ValueType>(vecAveragedParameter, 0.0);
}

/*! \brief calculate averaged inverse EM parameter
 *
 \param vecEMparameter EM parameter vector.
 \param vecInverseAveragedEMparameter Averaged inverse EM parameter vector which is calculated
 \param averagedMatrix Averaging matrix which is used to calculate averaged vector
 */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterEM<ValueType>::calculateInverseAveragedEMparameter(scai::lama::Vector<ValueType> const &vecEMparameter, scai::lama::DenseVector<ValueType> &vecInverseAveragedEMparameter, scai::lama::Matrix<ValueType> const &averagedMatrix)
{
    vecInverseAveragedEMparameter = averagedMatrix * vecEMparameter;
    vecInverseAveragedEMparameter = 1 / vecInverseAveragedEMparameter;

    Common::replaceInvalid<ValueType>(vecInverseAveragedEMparameter, 0.0);
}

/*! \brief Get const reference to averaged magneticPermeabilityEM in yz-plane
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getInverseMagneticPermeabilityEMAverageYZ()
{
    // If Averaging is outdated or has to be calculated for the first time, then recalculate averaging
    if (dirtyFlagAveraging == true) {
        calculateAveraging();
    }
    return (inverseMagneticPermeabilityEMAverageYZ);
}

/*! \brief Get const reference to averaged magneticPermeabilityEM in yz-plane
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getInverseMagneticPermeabilityEMAverageYZ() const
{
    SCAI_ASSERT(dirtyFlagAveraging == false, "Averaging has to be recalculated! ");
    return (inverseMagneticPermeabilityEMAverageYZ);
}

/*! \brief Get const reference to averaged magneticPermeabilityEM in xz-plane
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getInverseMagneticPermeabilityEMAverageXZ()
{
    // If Averaging is outdated or has to be calculated for the first time, then recalculate averaging
    if (dirtyFlagAveraging == true) {
        calculateAveraging();
    }
    return (inverseMagneticPermeabilityEMAverageXZ);
}

/*! \brief Get const reference to averaged magneticPermeabilityEM in xz-plane
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getInverseMagneticPermeabilityEMAverageXZ() const
{
    SCAI_ASSERT(dirtyFlagAveraging == false, "Averaging has to be recalculated! ");
    return (inverseMagneticPermeabilityEMAverageXZ);
}

/*! \brief Get const reference to averaged magneticPermeabilityEM in xy-plane
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getInverseMagneticPermeabilityEMAverageXY()
{
    // If Averaging is outdated or has to be calculated for the first time, then recalculate averaging
    if (dirtyFlagAveraging == true) {
        calculateAveraging();
    }
    return (inverseMagneticPermeabilityEMAverageXY);
}

/*! \brief Get const reference to averaged magneticPermeabilityEM in xy-plane
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getInverseMagneticPermeabilityEMAverageXY() const
{
    SCAI_ASSERT(dirtyFlagAveraging == false, "Averaging has to be recalculated! ");
    return (inverseMagneticPermeabilityEMAverageXY);
}

/*! \brief Get const reference to averaged conductivityEM in x-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getConductivityEMAverageX()
{
    // If Averaging is outdated or has to be calculated for the first time, then recalculate averaging
    if (dirtyFlagAveraging == true) {
        calculateAveraging();
    }
    return (conductivityEMAverageX);
}

/*! \brief Get const reference to averaged conductivityEM in x-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getConductivityEMAverageX() const
{
    SCAI_ASSERT(dirtyFlagAveraging == false, "Averaging has to be recalculated! ");
    return (conductivityEMAverageX);
}

/*! \brief Get const reference to averaged conductivityEM in y-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getConductivityEMAverageY()
{
    // If Averaging is outdated or has to be calculated for the first time, then recalculate averaging
    if (dirtyFlagAveraging == true) {
        calculateAveraging();
    }
    return (conductivityEMAverageY);
}

/*! \brief Get const reference to averaged conductivityEM in y-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getConductivityEMAverageY() const
{
    SCAI_ASSERT(dirtyFlagAveraging == false, "Averaging has to be recalculated! ");
    return (conductivityEMAverageY);
}

/*! \brief Get const reference to averaged conductivityEM in z-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getConductivityEMAverageZ()
{
    // If Averaging is outdated or has to be calculated for the first time, then recalculate averaging
    if (dirtyFlagAveraging == true) {
        calculateAveraging();
    }
    return (conductivityEMAverageZ);
}

/*! \brief Get const reference to averaged conductivityEM in z-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getConductivityEMAverageZ() const
{
    SCAI_ASSERT(dirtyFlagAveraging == false, "Averaging has to be recalculated! ");
    return (conductivityEMAverageZ);
}

/*! \brief Get const reference to averaged dielectricPermittivityEM in x-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getDielectricPermittivityEMAverageX()
{
    // If Averaging is outdated or has to be calculated for the first time, then recalculate averaging
    if (dirtyFlagAveraging == true) {
        calculateAveraging();
    }
    return (dielectricPermittivityEMAverageX);
}

/*! \brief Get const reference to averaged dielectricPermittivityEM in x-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getDielectricPermittivityEMAverageX() const
{
    SCAI_ASSERT(dirtyFlagAveraging == false, "Averaging has to be recalculated! ");
    return (dielectricPermittivityEMAverageX);
}

/*! \brief Get const reference to averaged dielectricPermittivityEM in y-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getDielectricPermittivityEMAverageY()
{
    // If Averaging is outdated or has to be calculated for the first time, then recalculate averaging
    if (dirtyFlagAveraging == true) {
        calculateAveraging();
    }
    return (dielectricPermittivityEMAverageY);
}

/*! \brief Get const reference to averaged dielectricPermittivityEM in y-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getDielectricPermittivityEMAverageY() const
{
    SCAI_ASSERT(dirtyFlagAveraging == false, "Averaging has to be recalculated! ");
    return (dielectricPermittivityEMAverageY);
}

/*! \brief Get const reference to averaged dielectricPermittivityEM in z-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getDielectricPermittivityEMAverageZ()
{
    // If Averaging is outdated or has to be calculated for the first time, then recalculate averaging
    if (dirtyFlagAveraging == true) {
        calculateAveraging();
    }
    return (dielectricPermittivityEMAverageZ);
}

/*! \brief Get const reference to averaged dielectricPermittivityEM in z-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getDielectricPermittivityEMAverageZ() const
{
    SCAI_ASSERT(dirtyFlagAveraging == false, "Averaging has to be recalculated! ");
    return (dielectricPermittivityEMAverageZ);
}

/*! \brief Get const reference to averaged conductivityEMoptical in x-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getConductivityEMopticalAverageX()
{
    // If Averaging is outdated or has to be calculated for the first time, then recalculate averaging
    if (dirtyFlagAveraging == true) {
        calculateAveraging();
    }
    return (conductivityEMopticalAverageX);
}

/*! \brief Get const reference to averaged conductivityEMoptical in x-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getConductivityEMopticalAverageX() const
{
    SCAI_ASSERT(dirtyFlagAveraging == false, "Averaging has to be recalculated! ");
    return (conductivityEMopticalAverageX);
}

/*! \brief Get const reference to averaged conductivityEMoptical in y-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getConductivityEMopticalAverageY()
{
    // If Averaging is outdated or has to be calculated for the first time, then recalculate averaging
    if (dirtyFlagAveraging == true) {
        calculateAveraging();
    }
    return (conductivityEMopticalAverageY);
}

/*! \brief Get const reference to averaged conductivityEMoptical in y-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getConductivityEMopticalAverageY() const
{
    SCAI_ASSERT(dirtyFlagAveraging == false, "Averaging has to be recalculated! ");
    return (conductivityEMopticalAverageY);
}

/*! \brief Get const reference to averaged conductivityEMoptical in z-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getConductivityEMopticalAverageZ()
{
    // If Averaging is outdated or has to be calculated for the first time, then recalculate averaging
    if (dirtyFlagAveraging == true) {
        calculateAveraging();
    }
    return (conductivityEMopticalAverageZ);
}

/*! \brief Get const reference to averaged conductivityEMoptical in z-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getConductivityEMopticalAverageZ() const
{
    SCAI_ASSERT(dirtyFlagAveraging == false, "Averaging has to be recalculated! ");
    return (conductivityEMopticalAverageZ);
}

/*! \brief Get const reference to averaged dielectricPermittivityEMoptical in x-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getDielectricPermittivityEMopticalAverageX()
{
    // If Averaging is outdated or has to be calculated for the first time, then recalculate averaging
    if (dirtyFlagAveraging == true) {
        calculateAveraging();
    }
    return (dielectricPermittivityEMopticalAverageX);
}

/*! \brief Get const reference to averaged dielectricPermittivityEMoptical in x-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getDielectricPermittivityEMopticalAverageX() const
{
    SCAI_ASSERT(dirtyFlagAveraging == false, "Averaging has to be recalculated! ");
    return (dielectricPermittivityEMopticalAverageX);
}

/*! \brief Get const reference to averaged dielectricPermittivityEMoptical in y-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getDielectricPermittivityEMopticalAverageY()
{
    // If Averaging is outdated or has to be calculated for the first time, then recalculate averaging
    if (dirtyFlagAveraging == true) {
        calculateAveraging();
    }
    return (dielectricPermittivityEMopticalAverageY);
}

/*! \brief Get const reference to averaged dielectricPermittivityEMoptical in y-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getDielectricPermittivityEMopticalAverageY() const
{
    SCAI_ASSERT(dirtyFlagAveraging == false, "Averaging has to be recalculated! ");
    return (dielectricPermittivityEMopticalAverageY);
}

/*! \brief Get const reference to averaged dielectricPermittivityEMoptical in z-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getDielectricPermittivityEMopticalAverageZ()
{
    // If Averaging is outdated or has to be calculated for the first time, then recalculate averaging
    if (dirtyFlagAveraging == true) {
        calculateAveraging();
    }
    return (dielectricPermittivityEMopticalAverageZ);
}

/*! \brief Get const reference to averaged dielectricPermittivityEMoptical in z-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getDielectricPermittivityEMopticalAverageZ() const
{
    SCAI_ASSERT(dirtyFlagAveraging == false, "Averaging has to be recalculated! ");
    return (dielectricPermittivityEMopticalAverageZ);
}

/*! \brief Get const reference to averaged tauDielectricPermittivityEM in x-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getTauDielectricPermittivityEMAverageX()
{
    // If Averaging is outdated or has to be calculated for the first time, then recalculate averaging
    if (dirtyFlagAveraging == true) {
        calculateAveraging();
    }
    return (tauDielectricPermittivityEMAverageX);
}

/*! \brief Get const reference to averaged tauDielectricPermittivityEM in x-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getTauDielectricPermittivityEMAverageX() const
{
    SCAI_ASSERT(dirtyFlagAveraging == false, "Averaging has to be recalculated! ");
    return (tauDielectricPermittivityEMAverageX);
}

/*! \brief Get const reference to averaged tauDielectricPermittivityEM in y-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getTauDielectricPermittivityEMAverageY()
{
    // If Averaging is outdated or has to be calculated for the first time, then recalculate averaging
    if (dirtyFlagAveraging == true) {
        calculateAveraging();
    }
    return (tauDielectricPermittivityEMAverageY);
}

/*! \brief Get const reference to averaged tauDielectricPermittivityEM in y-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getTauDielectricPermittivityEMAverageY() const
{
    SCAI_ASSERT(dirtyFlagAveraging == false, "Averaging has to be recalculated! ");
    return (tauDielectricPermittivityEMAverageY);
}

/*! \brief Get const reference to averaged tauDielectricPermittivityEM in z-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getTauDielectricPermittivityEMAverageZ()
{
    // If Averaging is outdated or has to be calculated for the first time, then recalculate averaging
    if (dirtyFlagAveraging == true) {
        calculateAveraging();
    }
    return (tauDielectricPermittivityEMAverageZ);
}

/*! \brief Get const reference to averaged tauDielectricPermittivityEM in z-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getTauDielectricPermittivityEMAverageZ() const
{
    SCAI_ASSERT(dirtyFlagAveraging == false, "Averaging has to be recalculated! ");
    return (tauDielectricPermittivityEMAverageZ);
}

/*! \brief Overloading = Operation
 *
 \param rhs Model which is copied.
 */
template <typename ValueType>
KITGPI::Modelparameter::ModelparameterEM<ValueType> &KITGPI::Modelparameter::ModelparameterEM<ValueType>::operator=(KITGPI::Modelparameter::ModelparameterEM<ValueType> const &rhs)
{
    assign(rhs);
    return *this;
}

/*! \brief Overloading += Operation
 *
 \param rhs Model which is added.
 */
template <typename ValueType>
KITGPI::Modelparameter::ModelparameterEM<ValueType> &KITGPI::Modelparameter::ModelparameterEM<ValueType>::operator+=(KITGPI::Modelparameter::ModelparameterEM<ValueType> const &rhs)
{
    plusAssign(rhs);
    return *this;
}

/*! \brief Overloading -= Operation
 *
 \param rhs Model which is substracted.
 */
template <typename ValueType>
KITGPI::Modelparameter::ModelparameterEM<ValueType> &KITGPI::Modelparameter::ModelparameterEM<ValueType>::operator-=(KITGPI::Modelparameter::ModelparameterEM<ValueType> const &rhs)
{
    minusAssign(rhs);
    return *this;
}

template class KITGPI::Modelparameter::ModelparameterEM<float>;
template class KITGPI::Modelparameter::ModelparameterEM<double>;
