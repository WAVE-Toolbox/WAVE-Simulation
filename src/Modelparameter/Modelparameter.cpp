#include "Modelparameter.hpp"
using namespace scai;
using namespace KITGPI;

/*! \brief Getter method for partitionedIn */
template <typename ValueType>
IndexType KITGPI::Modelparameter::Modelparameter<ValueType>::getPartitionedIn()
{
    return (PartitionedIn);
}

/*! \brief Getter method for partitionedOut */
template <typename ValueType>
IndexType KITGPI::Modelparameter::Modelparameter<ValueType>::getPartitionedOut()
{
    return (PartitionedOut);
}

/*! \brief Getter method for parametrisation */
template <typename ValueType>
IndexType KITGPI::Modelparameter::Modelparameter<ValueType>::getParametrisation()
{
    return (parametrisation);
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
 *  Reads a single model from an external mtx file.
 \param vector Singel modelparameter which will be initialized
 \param ctx Context
 \param dist Distribution
 \param filename Location of external file which will be read in
 \param partitionedIn Partitioned input
 */
template <typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::initModelparameter(scai::lama::Vector<ValueType> &vector, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, std::string filename, IndexType partitionedIn)
{
    HOST_PRINT(dist->getCommunicatorPtr(), "", "initModelParameter from file " << filename << "\n")
    allocateModelparameter(vector, ctx, dist);

    readModelparameter(vector, filename, dist, partitionedIn);

    vector.redistribute(dist);
}

/*! \brief Write singe modelparameter to an external file
 *
 *  Write a single model to an external file block.
 \param vector Single modelparameter which will be written to filename
 \param filename Name of file in which modelparameter will be written
 \param partitionedOut Partitioned output
 */
template <typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::writeModelparameter(scai::lama::Vector<ValueType> const &vector, std::string filename, IndexType partitionedOut) const
{
    PartitionedInOut::PartitionedInOut<ValueType> partitionOut;

    switch (partitionedOut) {
    case false:
        vector.writeToFile(filename);
        HOST_PRINT(vector.getDistributionPtr()->getCommunicatorPtr(), "writing " << filename << "\n");
        break;

    case true:
        partitionOut.writeToDistributedFiles(vector, filename);
        break;

    default:
        COMMON_THROWEXCEPTION("Unexpected output option!")
        break;
    }
};

/*! \brief Read a modelparameter from file
 */
template <typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::readModelparameter(scai::lama::Vector<ValueType> &vector, std::string filename, scai::dmemo::DistributionPtr dist, IndexType partitionedIn)
{

    PartitionedInOut::PartitionedInOut<ValueType> partitionIn;

    switch (partitionedIn) {
    case false:
        partitionIn.readFromOneFile(vector, filename, dist);
        break;

    case true:
        partitionIn.readFromDistributedFiles(vector, filename, dist);
        break;

    default:
        COMMON_THROWEXCEPTION("Unexpected input option!")
        break;
    }
};

/*! \brief Allocate a single modelparameter
 */
template <typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::allocateModelparameter(scai::lama::Vector<ValueType> &vector, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist)
{
    vector.setContextPtr(ctx);
    vector.allocate(dist);
};

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
void KITGPI::Modelparameter::Modelparameter<ValueType>::calcVelocityFromModulus(scai::lama::Vector<ValueType> &vectorModulus, scai::lama::Vector<ValueType> &vecDensity, scai::lama::Vector<ValueType> &vecVelocity)
{
    /* Modulus = pow(velocity,2) * Density */
    /* Velocity = sqrt( Modulus / Density )  */
    vecVelocity = vectorModulus / vecDensity; /* = Modulus / Density */
    vecVelocity = lama::sqrt(vecVelocity);    /* = sqrt( Modulus / Density ) */
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
    SCAI_ASSERT(dirtyFlagInverseDensity == false, "Inverse density has to be recalculated! ");
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
        dirtyFlagPWaveModulus = false;
        calcModulusFromVelocity(velocityP, density, pWaveModulus);
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

//! \brief Calculate density averaging matrix in x-direction
/*!
 *
 \param dist Distribution
 */
template <typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::calcDensityAverageMatrixX(Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::DistributionPtr dist)
{
    hmemo::HArray<IndexType> ownedIndexes; // all (global) points owned by this process
    dist->getOwnedIndexes(ownedIndexes);

    lama::MatrixAssembly<ValueType> assembly;
    assembly.reserve(ownedIndexes.size() * 2);

    for (IndexType ownedIndex : hmemo::hostReadAccess(ownedIndexes)) {

        Acquisition::coordinate3D coordinate = modelCoordinates.index2coordinate(ownedIndex);
        scai::IndexType X = coordinate.x;
        if ((X + modelCoordinates.getDHFactor(coordinate)*0.5) == modelCoordinates.getNX() - 0.5)
            assembly.push(ownedIndex, ownedIndex, 1.0);
        else
            assembly.push(ownedIndex, ownedIndex, 1.0 / 2.0);

        X+=modelCoordinates.getDHFactor(coordinate);
        if (X < modelCoordinates.getNX()) {
            IndexType columnIndex = modelCoordinates.coordinate2index(X, coordinate.y, coordinate.z);
            assembly.push(ownedIndex, columnIndex, 1.0 / 2.0);
        }
    }

    DensityAverageMatrixX = lama::zero<SparseFormat>(dist, dist);
    DensityAverageMatrixX.fillFromAssembly(assembly);
    DensityAverageMatrixX.writeToFile("test.mtx");
}

//! \brief Calculate density averaging matrix in y-direction
/*!
 *
 \param dist Distribution
 */
template <typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::calcDensityAverageMatrixY(Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::DistributionPtr dist)
{

    hmemo::HArray<IndexType> ownedIndexes; // all (global) points owned by this process
    dist->getOwnedIndexes(ownedIndexes);

    lama::MatrixAssembly<ValueType> assembly;
    assembly.reserve(ownedIndexes.size() * 2);
    IndexType layer=0;
    
    for (IndexType ownedIndex : hmemo::hostReadAccess(ownedIndexes)) {

        
        Acquisition::coordinate3D coordinate = modelCoordinates.index2coordinate(ownedIndex);
        
        layer=modelCoordinates.getLayer(coordinate);
        
        scai::IndexType Y = coordinate.y;
        if ((Y + modelCoordinates.getDHFactor(layer)*0.5) == modelCoordinates.getNY() - 0.5)
            assembly.push(ownedIndex, ownedIndex, 1.0);
        else
            assembly.push(ownedIndex, ownedIndex, 1.0 / 2.0);

        Y+=modelCoordinates.getDHFactor(layer);
        if ((modelCoordinates.locatedOnInterface(coordinate)) && (!modelCoordinates.getTransition(coordinate)))
        Y+=modelCoordinates.getDHFactor(layer-1)-modelCoordinates.getDHFactor(layer);  
        
        if (Y < modelCoordinates.getNY()) {
            IndexType columnIndex = modelCoordinates.coordinate2index(coordinate.x, Y, coordinate.z);
            assembly.push(ownedIndex, columnIndex, 1.0 / 2.0);
        }
    }

    DensityAverageMatrixY = lama::zero<SparseFormat>(dist, dist);
    DensityAverageMatrixY.fillFromAssembly(assembly);

}

//! \brief Calculate density averaging matrix in z-direction
/*!
 *
 \param dist Distribution
 */
template <typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::calcDensityAverageMatrixZ(Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::DistributionPtr dist)
{
    hmemo::HArray<IndexType> ownedIndexes; // all (global) points owned by this process
    dist->getOwnedIndexes(ownedIndexes);

    lama::MatrixAssembly<ValueType> assembly;
    assembly.reserve(ownedIndexes.size() * 2);

    for (IndexType ownedIndex : hmemo::hostReadAccess(ownedIndexes)) {

        Acquisition::coordinate3D coordinate = modelCoordinates.index2coordinate(ownedIndex);
        scai::IndexType Z = coordinate.z;
        if ((Z + modelCoordinates.getDHFactor(coordinate)*0.5) == modelCoordinates.getNZ() - 0.5)
            assembly.push(ownedIndex, ownedIndex, 1.0);
        else
            assembly.push(ownedIndex, ownedIndex, 1.0 / 2.0);

        Z+= modelCoordinates.getDHFactor(coordinate);
        if (Z < modelCoordinates.getNZ()) {
            IndexType columnIndex = modelCoordinates.coordinate2index(coordinate.x, coordinate.y, Z);
            assembly.push(ownedIndex, columnIndex, 1.0 / 2.0);
        }
    }

    DensityAverageMatrixZ = lama::zero<SparseFormat>(dist, dist);
    DensityAverageMatrixZ.fillFromAssembly(assembly);
}

//! \brief Calculate s-wave modulus averaging matrix in x-direction
/*!
 *
 \param dist Distribution
 */
template <typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::calcSWaveModulusAverageMatrixXY(Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::DistributionPtr dist)
{
    hmemo::HArray<IndexType> ownedIndexes; // all (global) points owned by this process
    dist->getOwnedIndexes(ownedIndexes);

    lama::MatrixAssembly<ValueType> assembly;
    assembly.reserve(ownedIndexes.size() * 4);

    for (IndexType ownedIndex : hmemo::hostReadAccess(ownedIndexes)) {

        Acquisition::coordinate3D coordinate = modelCoordinates.index2coordinate(ownedIndex);

        scai::IndexType X = coordinate.x;
        scai::IndexType Y = coordinate.y;
        IndexType columnIndex;

        /*Points
            1      2
               av
            3      4
        */
        // Point 1 (diagonal element)
        if ((X == modelCoordinates.getNX() - 1) && (Y == modelCoordinates.getNY() - 1))
            assembly.push(ownedIndex, ownedIndex, 1.0);
        else if (!(X == modelCoordinates.getNX() - 1) != !(Y == modelCoordinates.getNY() - 1))
            assembly.push(ownedIndex, ownedIndex, 1.0 / 2.0);
        else
            assembly.push(ownedIndex, ownedIndex, 1.0 / 4.0);

        //Point 2
        X = coordinate.x + 1;
        Y = coordinate.y;
        if (X < modelCoordinates.getNX()) {
            columnIndex = modelCoordinates.coordinate2index(X, Y, coordinate.z);
            if (Y == modelCoordinates.getNY() - 1)
                assembly.push(ownedIndex, columnIndex, 1.0 / 2.0);
            else
                assembly.push(ownedIndex, columnIndex, 1.0 / 4.0);
        }

        //Point 3
        X = coordinate.x;
        Y = coordinate.y + 1;
        if (Y < modelCoordinates.getNY()) {
            columnIndex = modelCoordinates.coordinate2index(X, Y, coordinate.z);
            if (X == modelCoordinates.getNX() - 1)
                assembly.push(ownedIndex, columnIndex, 1.0 / 2.0);
            else
                assembly.push(ownedIndex, columnIndex, 1.0 / 4.0);
        }

        //Point 4
        X = coordinate.x + 1;
        Y = coordinate.y + 1;
        if ((X < modelCoordinates.getNX()) && (Y < modelCoordinates.getNY())) {
            columnIndex = modelCoordinates.coordinate2index(X, Y, coordinate.z);
            assembly.push(ownedIndex, columnIndex, 1.0 / 4.0);
        }
    }

    sWaveModulusAverageMatrixXY = lama::zero<SparseFormat>(dist, dist);
    sWaveModulusAverageMatrixXY.fillFromAssembly(assembly);
}

//! \brief Calculate s-wave modulus averaging matrix in y-direction
/*!
 *
 \param dist Distribution
 */
template <typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::calcSWaveModulusAverageMatrixXZ(Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::DistributionPtr dist)
{
    //     calcAverageMatrix(sWaveModulusAverageMatrixXZ, &Modelparameter<ValueType>::calcNumberRowElements_SWaveModulusAverageMatrixXZ, &Modelparameter<ValueType>::setRowElements_SWaveModulusAverageMatrixXZ, modelCoordinates, dist);

    hmemo::HArray<IndexType> ownedIndexes; // all (global) points owned by this process
    dist->getOwnedIndexes(ownedIndexes);

    lama::MatrixAssembly<ValueType> assembly;
    assembly.reserve(ownedIndexes.size() * 4);

    for (IndexType ownedIndex : hmemo::hostReadAccess(ownedIndexes)) {

        Acquisition::coordinate3D coordinate = modelCoordinates.index2coordinate(ownedIndex);

        scai::IndexType X = coordinate.x;
        scai::IndexType Z = coordinate.z;
        IndexType columnIndex;

        /*Points
            1      2
               av
            3      4
        */
        // Point 1 (diagonal element)
        if ((X == modelCoordinates.getNX() - 1) && (Z == modelCoordinates.getNZ() - 1))
            assembly.push(ownedIndex, ownedIndex, 1.0);
        else if (!(X == modelCoordinates.getNX() - 1) != !(Z == modelCoordinates.getNZ() - 1))
            assembly.push(ownedIndex, ownedIndex, 1.0 / 2.0);
        else
            assembly.push(ownedIndex, ownedIndex, 1.0 / 4.0);

        //Point 2
        X = coordinate.x + 1;
        Z = coordinate.z;
        if (X < modelCoordinates.getNX()) {
            columnIndex = modelCoordinates.coordinate2index(X, coordinate.y, Z);
            if (Z == modelCoordinates.getNZ() - 1)
                assembly.push(ownedIndex, columnIndex, 1.0 / 2.0);
            else
                assembly.push(ownedIndex, columnIndex, 1.0 / 4.0);
        }

        //Point 3
        X = coordinate.x;
        Z = coordinate.z + 1;
        if (Z < modelCoordinates.getNZ()) {
            columnIndex = modelCoordinates.coordinate2index(X, coordinate.y, Z);
            if (X == modelCoordinates.getNX() - 1)
                assembly.push(ownedIndex, columnIndex, 1.0 / 2.0);
            else
                assembly.push(ownedIndex, columnIndex, 1.0 / 4.0);
        }

        //Point 4
        X = coordinate.x + 1;
        Z = coordinate.z + 1;
        if ((X < modelCoordinates.getNX()) && (Z < modelCoordinates.getNZ())) {
            columnIndex = modelCoordinates.coordinate2index(X, coordinate.y, Z);
            assembly.push(ownedIndex, columnIndex, 1.0 / 4.0);
        }
    }

    sWaveModulusAverageMatrixXZ = lama::zero<SparseFormat>(dist, dist);
    sWaveModulusAverageMatrixXZ.fillFromAssembly(assembly);
}

//! \brief Calculate s-wave modulus averaging matrix in z-direction
/*!
 *
 \param dist Distribution
 */
template <typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::calcSWaveModulusAverageMatrixYZ(Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::DistributionPtr dist)
{
    hmemo::HArray<IndexType> ownedIndexes; // all (global) points owned by this process
    dist->getOwnedIndexes(ownedIndexes);

    lama::MatrixAssembly<ValueType> assembly;
    assembly.reserve(ownedIndexes.size() * 4);

    for (IndexType ownedIndex : hmemo::hostReadAccess(ownedIndexes)) {

        Acquisition::coordinate3D coordinate = modelCoordinates.index2coordinate(ownedIndex);

        scai::IndexType Y = coordinate.y;
        scai::IndexType Z = coordinate.z;
        IndexType columnIndex;

        /*Points
            1      2
               av
            3      4
        */
        // Point 1 (diagonal element)
        if ((Y == modelCoordinates.getNY() - 1) && (Z == modelCoordinates.getNZ() - 1))
            assembly.push(ownedIndex, ownedIndex, 1.0);
        else if (!(Y == modelCoordinates.getNY() - 1) != !(Z == modelCoordinates.getNZ() - 1))
            assembly.push(ownedIndex, ownedIndex, 1.0 / 2.0);
        else
            assembly.push(ownedIndex, ownedIndex, 1.0 / 4.0);

        //Point 2
        Y = coordinate.y + 1;
        Z = coordinate.z;
        if (Y < modelCoordinates.getNY()) {
            columnIndex = modelCoordinates.coordinate2index(coordinate.x, Y, Z);
            if (Z == modelCoordinates.getNZ() - 1)
                assembly.push(ownedIndex, columnIndex, 1.0 / 2.0);
            else
                assembly.push(ownedIndex, columnIndex, 1.0 / 4.0);
        }

        //Point 3
        Y = coordinate.y;
        Z = coordinate.z + 1;
        if (Z < modelCoordinates.getNZ()) {
            columnIndex = modelCoordinates.coordinate2index(coordinate.x, Y, Z);
            if (Y == modelCoordinates.getNY() - 1)
                assembly.push(ownedIndex, columnIndex, 1.0 / 2.0);
            else
                assembly.push(ownedIndex, columnIndex, 1.0 / 4.0);
        }

        //Point 4
        Y = coordinate.y + 1;
        Z = coordinate.z + 1;
        if ((Y < modelCoordinates.getNY()) && (Z < modelCoordinates.getNZ())) {
            columnIndex = modelCoordinates.coordinate2index(coordinate.x, Y, Z);
            assembly.push(ownedIndex, columnIndex, 1.0 / 4.0);
        }
    }

    sWaveModulusAverageMatrixYZ = lama::zero<SparseFormat>(dist, dist);
    sWaveModulusAverageMatrixYZ.fillFromAssembly(assembly);
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
    Common::searchAndReplace<ValueType>(vecSWaveModulus, 1.0, 1.0, 1);

    vecAvSWaveModulus = 1 / vecSWaveModulus;
    auto temp = lama::eval<lama::DenseVector<ValueType>>(avSWaveModulusMatrix * vecAvSWaveModulus);
    vecAvSWaveModulus = 1 / temp;

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
