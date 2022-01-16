#include "Modelparameter.hpp"
#include "../IO/IO.hpp"

using namespace scai;
using namespace KITGPI;

/*! \brief Getter method for memory usage */
template <typename ValueType>
ValueType KITGPI::Modelparameter::Modelparameter<ValueType>::getMemoryUsage(scai::dmemo::DistributionPtr dist, scai::IndexType numParameter)
{
    ValueType size = getMemoryModel(dist) / 1024 / 1024 * numParameter;
    return size;
}

/*! \brief Getter method for center frequency of CPML */
template <typename ValueType>
ValueType const KITGPI::Modelparameter::Modelparameter<ValueType>::getCenterFrequencyCPML() const
{
    return (centerFrequencyCPML);
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

/*! \brief Getter method for parameterisation */
template <typename ValueType>
bool KITGPI::Modelparameter::Modelparameter<ValueType>::getEffectiveParameterisation() const
{
    return (effectiveParameterisation);
}

/*! \brief Set method for parameterisation */
template <typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::setEffectiveParameterisation(bool const setEffectiveParameterisation)
{
    effectiveParameterisation = setEffectiveParameterisation;
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

/*! \brief Getter method for gradientType */
template <typename ValueType>
IndexType KITGPI::Modelparameter::Modelparameter<ValueType>::getGradientType() const
{
    return (gradientType);
}

/*! \brief Set method for gradientType */
template <typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::setGradientType(scai::IndexType const setGradientType)
{
    gradientType = setGradientType;
}

/*! \brief Getter method for decomposition */
template <typename ValueType>
IndexType KITGPI::Modelparameter::Modelparameter<ValueType>::getDecomposeType() const
{
    return (decomposition);
}

/*! \brief Set method for decomposition */
template <typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::setDecomposeType(scai::IndexType const setDecomposeType)
{
    decomposition = setDecomposeType;
}

/*! \brief Getter method for vmim */
template <typename ValueType>
ValueType KITGPI::Modelparameter::Modelparameter<ValueType>::getVmin()
{    
    scai::lama::DenseVector<ValueType> vMinVec;
    bool isSeismic = Common::checkEquationType<ValueType>(equationType);
    if (isSeismic) {
        if (equationType == "sh" || equationType == "viscosh" || equationType == "elastic" || equationType == "viscoelastic") {
            vMinVec = this->getVelocityS();
        } else if (equationType == "acoustic") {
            vMinVec = this->getVelocityP();
        }
    } else {
        vMinVec = this->getVelocityEM();
    }
    ValueType vMin = vMinVec.min();
    if (vMin <= 1.0) {        
        Common::searchAndReplace<ValueType>(vMinVec, 1.0, vMinVec.max(), 3);
        vMin = vMinVec.min();
    }
    
    return (vMin);
}

/*! \brief Getter method for compensation */
template <typename ValueType>
scai::lama::DenseVector<ValueType> KITGPI::Modelparameter::Modelparameter<ValueType>::getCompensation(ValueType DT, scai::IndexType tStep) const
{    
    scai::lama::DenseVector<ValueType> compensation;
    bool isSeismic = Common::checkEquationType<ValueType>(equationType);
    if (isSeismic) {
        COMMON_THROWEXCEPTION("There is no compensation in an Seismic modelling")
    } else {
        compensation = this->getElectricConductivity();
        compensation /= this->getDielectricPermittivity();
        compensation *= tStep * DT;
        compensation.unaryOp(compensation, common::UnaryOp::EXP);
    }
    
    return (compensation);
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

/*! \brief Get shrink-vector that shrinks the old values in the big model
 \param dist Distribution of the pershot
 \param distBig Distribution of the big model
 \param modelCoordinates coordinate class object of the pershot
 \param modelCoordinatesBig coordinate class object of the big model
 \param cutCoordinate coordinate where to cut the pershot
 */
template <typename ValueType>
scai::lama::SparseVector<ValueType> KITGPI::Modelparameter::Modelparameter<ValueType>::getShrinkVector(scai::dmemo::DistributionPtr dist, scai::dmemo::DistributionPtr distBig, Acquisition::Coordinates<ValueType> const &modelCoordinates, Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, Acquisition::coordinate3D const cutCoordinate, scai::IndexType boundaryWidthLeft, scai::IndexType boundaryWidthRight)
{
    scai::lama::SparseVector<ValueType> shrinkVector(distBig, 0.0); //!< Shrink Multiplication matrix
      
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
        assembly.push(indexBig, 1.0);
    }
    shrinkVector.fillFromAssembly(assembly);
    
    // damp the boundary boarders
    for (IndexType y = 0; y < modelCoordinatesBig.getNY(); y++) {
        for (IndexType i = 0; i < boundaryWidthLeft; i++) {
            ValueType tmp = sin(i * M_PI / ((ValueType)boundaryWidthLeft * 2.0));
            tmp *= tmp;
            shrinkVector[modelCoordinatesBig.coordinate2index(cutCoordinate.x+i, y, 0)] = tmp;
        }
        for (IndexType i = 0; i < boundaryWidthRight; i++) {
            ValueType tmp = sin(i * M_PI / ((ValueType)boundaryWidthRight * 2.0));
            tmp *= tmp;
            shrinkVector[modelCoordinatesBig.coordinate2index(cutCoordinate.x+modelCoordinates.getNX()-1-i, y, 0)] = tmp;
        }
    }
    return shrinkVector;
}

/*! \brief Getter method for relaxation frequency */
template <typename ValueType>
std::vector<ValueType> KITGPI::Modelparameter::Modelparameter<ValueType>::getRelaxationFrequency() const
{
    return (relaxationFrequency);
}

/*! \brief Getter method for number of relaxation mechanisms */
template <typename ValueType>
IndexType KITGPI::Modelparameter::Modelparameter<ValueType>::getNumRelaxationMechanisms() const
{
    return (numRelaxationMechanisms);
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

/*! \brief Get const reference to reflectivity
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Modelparameter<ValueType>::getReflectivity() const
{
    return (reflectivity);
}

/*! \brief Set reflectivity
 */
template <typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::setReflectivity(scai::lama::Vector<ValueType> const &setReflectivity)
{
    reflectivity = setReflectivity;
}

/*! \brief reset reflectivity
 */
template <typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::resetReflectivity()
{
    reflectivity = 0.0;
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
void KITGPI::Modelparameter::Modelparameter<ValueType>::calcInverseAveragedParameter(scai::lama::Vector<ValueType> const &vecParameter, scai::lama::DenseVector<ValueType> &vecInverseAveragedParameter, scai::lama::Matrix<ValueType> const &averagedMatrix)
{
    vecInverseAveragedParameter = averagedMatrix * vecParameter;
    vecInverseAveragedParameter = 1 / vecInverseAveragedParameter;

    Common::replaceInvalid<ValueType>(vecInverseAveragedParameter, 0.0);
}

/*! \brief calculate averaged tauS
 *
 \param vecTauS TauS vector
 \param vecAvTauS Averaged tauS vector which is calculated
 \param avTauSMatrix Averaging matrix which is used to calculate averaged vector
 */
template <typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::calcAveragedParameter(scai::lama::Vector<ValueType> const &vecParameter, scai::lama::DenseVector<ValueType> &vecAveragedParameter, scai::lama::Matrix<ValueType> const &averagedMatrix)
{
    vecAveragedParameter = averagedMatrix * vecParameter;
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
