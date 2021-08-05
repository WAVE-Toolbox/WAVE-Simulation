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

/*! \brief Getter method for gradientType */
template <typename ValueType>
IndexType KITGPI::Modelparameter::ModelparameterEM<ValueType>::getGradientType() const
{
    return (gradientType);
}

/*! \brief Set method for gradientType */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterEM<ValueType>::setGradientType(scai::IndexType const setGradientType)
{
    gradientType = setGradientType;
}

/*! \brief Getter method for decomposeType */
template <typename ValueType>
IndexType KITGPI::Modelparameter::ModelparameterEM<ValueType>::getDecomposeType() const
{
    return (decomposeType);
}

/*! \brief Set method for decomposeType */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterEM<ValueType>::setDecomposeType(scai::IndexType const setDecomposeType)
{
    decomposeType = setDecomposeType;
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
    setParameterisation(config.getAndCatch("parameterisation", 0));
    setInversionType(config.getAndCatch("inversionType", 1));
    setGradientType(config.getAndCatch("gradientType", 0));
    setDecomposeType(config.getAndCatch("decomposeType", 0));
    calcElectricConductivityReference(config.get<ValueType>("CenterFrequencyCPML"));
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
    dirtyFlagDielectricPermittivityOptical = true; // the visco Modulus vector is now dirty
    dirtyFlagElectricConductivityOptical = true; // the visco Modulus vector is now dirty
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
    dirtyFlagDielectricPermittivityOptical = true; // the visco Modulus vector is now dirty
    dirtyFlagElectricConductivityOptical = true; // the visco Modulus vector is now dirty
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
void KITGPI::Modelparameter::ModelparameterEM<ValueType>::writeRockMatrixParameter(std::string filename, scai::IndexType fileFormat) const
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
    dirtyFlagElectricConductivityOptical = true; // the visco Modulus vector is now dirty
    dirtyFlagDielectricPermittivityOptical = true; // the visco Modulus vector is now dirty
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
    dirtyFlagElectricConductivityOptical = true; // the visco Modulus vector is now dirty
    dirtyFlagDielectricPermittivityOptical = true; // the visco Modulus vector is now dirty
    dielectricPermittivity = setDielectricPermittivity;
}

/*! \brief Get const reference to electricConductivityOptical
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getElectricConductivityOptical()
{
    // If the modulus is dirty, then recalculate
    if (dirtyFlagElectricConductivityOptical) {
        HOST_PRINT(electricConductivity.getDistributionPtr()->getCommunicatorPtr(), "electricConductivityOptical will be calculated from\n dielectricPermittivity, electricConductivity, tauDielectricPermittivity\n, tauElectricDisplacement and numRelaxationMechanisms \n");              
        electricConductivityOptical = tauDielectricPermittivity / tauElectricDisplacement;    
        electricConductivityOptical *= numRelaxationMechanisms; // numRelaxationMechanisms = 1
        electricConductivityOptical /= numRelaxationMechanisms;
        electricConductivityOptical *= dielectricPermittivity;
        electricConductivityOptical += electricConductivity;
        dirtyFlagElectricConductivityOptical = false;
    }
    
    return (electricConductivityOptical);
}
/*! \brief Get const reference to electricConductivityOptical
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getElectricConductivityOptical() const
{
    SCAI_ASSERT(dirtyFlagElectricConductivityOptical == false, "electricConductivityOptical has to be recalculated! ");
    
    return (electricConductivityOptical);
}

/*! \brief Set electricConductivityOptical model parameter
 */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterEM<ValueType>::setElectricConductivityOptical(scai::lama::Vector<ValueType> const &setElectricConductivityOptical)
{
    dirtyFlagAveraging = true;    // If modulus will be changed, averaging needs to be redone
    electricConductivityOptical = setElectricConductivityOptical;
}

/*! \brief Get const reference to dielectricPermittivityOptical
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getDielectricPermittivityOptical()
{
    // If the modulus is dirty, then recalculate
    if (dirtyFlagDielectricPermittivityOptical) {
        HOST_PRINT(dielectricPermittivity.getDistributionPtr()->getCommunicatorPtr(), "dielectricPermittivityOptical will be calculated from\n dielectricPermittivity, electricConductivity, tauDielectricPermittivity\n, tauElectricConductivity and numRelaxationMechanisms \n");
        scai::lama::DenseVector<ValueType> temp;
        temp = numRelaxationMechanisms * tauDielectricPermittivity;  // numRelaxationMechanisms = 1
        temp /= numRelaxationMechanisms;
        dielectricPermittivityOptical = 1 - temp;      
        dielectricPermittivityOptical *= dielectricPermittivity;
        temp = electricConductivity * tauElectricConductivity;
        dielectricPermittivityOptical += temp;
        
        // optical dielectric permittivity should be larger than the dielectric permittivity of vacuum
        Common::searchAndReplace<ValueType>(dielectricPermittivityOptical, DielectricPermittivityVacuum, DielectricPermittivityVacuum, 1); 
        
        dirtyFlagDielectricPermittivityOptical = false;
    }
    
    return (dielectricPermittivityOptical);
}

/*! \brief Get const reference to dielectricPermittivityOptical
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getDielectricPermittivityOptical() const
{
    SCAI_ASSERT(dirtyFlagDielectricPermittivityOptical == false, "dielectricPermittivityOptical has to be recalculated! ");
    return (dielectricPermittivityOptical);
}

/*! \brief Set dielectricPermittivityOptical model parameter
 */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterEM<ValueType>::setDielectricPermittivityOptical(scai::lama::Vector<ValueType> const &setDielectricPermittivityOptical)
{
    dirtyFlagAveraging = true;    // If modulus will be changed, averaging needs to be redone
    dielectricPermittivityOptical = setDielectricPermittivityOptical;
}

/*! \brief Get const reference to tauElectricConductivity */
template <typename ValueType>
ValueType const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getTauElectricDisplacement()
{
    tauElectricDisplacement = 1.0 / (2.0 * M_PI * this->getRelaxationFrequency()); 
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
    dirtyFlagDielectricPermittivityOptical = true; // the visco Modulus vector is now dirty
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
    dirtyFlagElectricConductivityOptical = true; // the visco Modulus vector is now dirty
    dirtyFlagDielectricPermittivityOptical = true; // the visco Modulus vector is now dirty
    tauDielectricPermittivity = setTauDielectricPermittivity;
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

/*! \brief Get const reference to reflectivity
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getReflectivity() const
{
    return (reflectivity);
}

/*! \brief Set reflectivity
 */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterEM<ValueType>::setReflectivity(scai::lama::Vector<ValueType> const &setReflectivity)
{
    reflectivity = setReflectivity;
}

/*! \brief reset reflectivity
 */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterEM<ValueType>::resetReflectivity()
{
    reflectivity = 0.0;
}

/*! \brief calculate reflectivity from permittivity
 */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterEM<ValueType>::calcReflectivity(Acquisition::Coordinates<ValueType> const &modelCoordinates, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives, ValueType DT)
{
    scai::lama::DenseVector<ValueType> dielectricPermittivitytemp;
    scai::lama::DenseVector<ValueType> dielectricPermittivityAverageYtemp;
    auto dist = dielectricPermittivity.getDistributionPtr();
    auto ctx = dielectricPermittivity.getContextPtr();
    auto const &Dyf = derivatives.getDyf();
    this->calcAverageMatrixY(modelCoordinates, dist);
    averageMatrixY.setContextPtr(ctx);
    dielectricPermittivitytemp = scai::lama::sqrt(dielectricPermittivity);
    this->calculateAveragedParameter(dielectricPermittivitytemp, dielectricPermittivityAverageYtemp, averageMatrixY);
    averageMatrixY.purge();
    dielectricPermittivityAverageYtemp *= 2;
    reflectivity = Dyf * dielectricPermittivitytemp;
    reflectivity *= -modelCoordinates.getDH() / DT;
    reflectivity /= dielectricPermittivityAverageYtemp;
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
    //     calcAverageMatrix(dielectricPermittivityAverageMatrixXZ, &ModelparameterEM<ValueType>::calcNumberRowElements_InverseDielectricPermittivityAverageMatrixXZ, &ModelparameterEM<ValueType>::setRowElements_InverseDielectricPermittivityAverageMatrixXZ, modelCoordinates, dist);

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

/*! \brief calculate averaged Parameter
 *
 \param vecParameter EM parameter vector
 \param vecAveragedParameter Averaged EM parameter vector which is calculated
 \param averagedMatrix Averaging matrix which is used to calculate averaged vector
 */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterEM<ValueType>::calculateAveragedParameter(scai::lama::Vector<ValueType> const &vecParameter, scai::lama::DenseVector<ValueType> &vecAveragedParameter, scai::lama::Matrix<ValueType> const &averagedMatrix)
{
    vecAveragedParameter = averagedMatrix * vecParameter;
    
    Common::replaceInvalid<ValueType>(vecAveragedParameter, 0.0);
}

/*! \brief calculate averaged inverse EM parameter
 *
 \param vecParameter EM parameter vector.
 \param vecInverseAveragedParameter Averaged inverse EM parameter vector which is calculated
 \param averagedMatrix Averaging matrix which is used to calculate averaged vector
 */
template <typename ValueType>
void KITGPI::Modelparameter::ModelparameterEM<ValueType>::calculateInverseAveragedParameter(scai::lama::Vector<ValueType> const &vecParameter, scai::lama::DenseVector<ValueType> &vecInverseAveragedParameter, scai::lama::Matrix<ValueType> const &averagedMatrix)
{
    vecInverseAveragedParameter = averagedMatrix * vecParameter;
    vecInverseAveragedParameter = 1 / vecInverseAveragedParameter;

    Common::replaceInvalid<ValueType>(vecInverseAveragedParameter, 0.0);
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

/*! \brief Get const reference to averaged electricConductivityOptical in x-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getElectricConductivityOpticalAverageX()
{
    // If Averaging is outdated or has to be calculated for the first time, then recalculate averaging
    if (dirtyFlagAveraging == true) {
        calculateAveraging();
    }
    return (electricConductivityOpticalAverageX);
}

/*! \brief Get const reference to averaged electricConductivityOptical in x-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getElectricConductivityOpticalAverageX() const
{
    SCAI_ASSERT(dirtyFlagAveraging == false, "Averaging has to be recalculated! ");
    return (electricConductivityOpticalAverageX);
}

/*! \brief Get const reference to averaged electricConductivityOptical in y-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getElectricConductivityOpticalAverageY()
{
    // If Averaging is outdated or has to be calculated for the first time, then recalculate averaging
    if (dirtyFlagAveraging == true) {
        calculateAveraging();
    }
    return (electricConductivityOpticalAverageY);
}

/*! \brief Get const reference to averaged electricConductivityOptical in y-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getElectricConductivityOpticalAverageY() const
{
    SCAI_ASSERT(dirtyFlagAveraging == false, "Averaging has to be recalculated! ");
    return (electricConductivityOpticalAverageY);
}

/*! \brief Get const reference to averaged electricConductivityOptical in z-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getElectricConductivityOpticalAverageZ()
{
    // If Averaging is outdated or has to be calculated for the first time, then recalculate averaging
    if (dirtyFlagAveraging == true) {
        calculateAveraging();
    }
    return (electricConductivityOpticalAverageZ);
}

/*! \brief Get const reference to averaged electricConductivityOptical in z-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getElectricConductivityOpticalAverageZ() const
{
    SCAI_ASSERT(dirtyFlagAveraging == false, "Averaging has to be recalculated! ");
    return (electricConductivityOpticalAverageZ);
}

/*! \brief Get const reference to averaged dielectricPermittivityOptical in x-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getDielectricPermittivityOpticalAverageX()
{
    // If Averaging is outdated or has to be calculated for the first time, then recalculate averaging
    if (dirtyFlagAveraging == true) {
        calculateAveraging();
    }
    return (dielectricPermittivityOpticalAverageX);
}

/*! \brief Get const reference to averaged dielectricPermittivityOptical in x-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getDielectricPermittivityOpticalAverageX() const
{
    SCAI_ASSERT(dirtyFlagAveraging == false, "Averaging has to be recalculated! ");
    return (dielectricPermittivityOpticalAverageX);
}

/*! \brief Get const reference to averaged dielectricPermittivityOptical in y-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getDielectricPermittivityOpticalAverageY()
{
    // If Averaging is outdated or has to be calculated for the first time, then recalculate averaging
    if (dirtyFlagAveraging == true) {
        calculateAveraging();
    }
    return (dielectricPermittivityOpticalAverageY);
}

/*! \brief Get const reference to averaged dielectricPermittivityOptical in y-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getDielectricPermittivityOpticalAverageY() const
{
    SCAI_ASSERT(dirtyFlagAveraging == false, "Averaging has to be recalculated! ");
    return (dielectricPermittivityOpticalAverageY);
}

/*! \brief Get const reference to averaged dielectricPermittivityOptical in z-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getDielectricPermittivityOpticalAverageZ()
{
    // If Averaging is outdated or has to be calculated for the first time, then recalculate averaging
    if (dirtyFlagAveraging == true) {
        calculateAveraging();
    }
    return (dielectricPermittivityOpticalAverageZ);
}

/*! \brief Get const reference to averaged dielectricPermittivityOptical in z-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::ModelparameterEM<ValueType>::getDielectricPermittivityOpticalAverageZ() const
{
    SCAI_ASSERT(dirtyFlagAveraging == false, "Averaging has to be recalculated! ");
    return (dielectricPermittivityOpticalAverageZ);
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
