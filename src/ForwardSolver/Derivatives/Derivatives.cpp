#include "Derivatives.hpp"
#include "../Common/Common.hpp"

using namespace scai;

//! \brief Setup configuration of the derivative object
/*!
 *
 \param config Configuration
 */
template <typename ValueType>
void KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::setup(Configuration::Configuration const &config)
{
    useFreeSurface = config.get<IndexType>("FreeSurface");
    useVarGrid = config.get<bool>("useVariableGrid");

    std::string type = config.get<std::string>("equationType");
    std::transform(type.begin(), type.end(), type.begin(), ::tolower);
    if ((type.compare("sh") == 0) || (type.compare("viscosh") == 0) || (type.compare("elastic") == 0) || (type.compare("viscoelastic") == 0)) {
        isElastic = true;
    }

    DT = config.get<ValueType>("DT");
    setFDCoef();

    try {
        useStencilMatrix = config.get<bool>("useStencilMatrix");
    } catch (...) {
        useStencilMatrix = false;
    }

    if (config.get<IndexType>("partitioning") != 1 && useStencilMatrix == true) {
        COMMON_THROWEXCEPTION("It is not possible to use the stencil matrix without grid distribution (partitioning=1)!")
    }

    try {
        useHybridFreeSurface = config.get<bool>("useHybridFreeSurface");
    } catch (...) {
        useHybridFreeSurface = false;
    }

    if (useStencilMatrix == false && useHybridFreeSurface == true) {
        COMMON_THROWEXCEPTION("It is not possible to use the hybrid matrix without stencil matrix!")
    }

    if ((!useStencilMatrix) && (config.get<bool>("useVariableFDoperators"))) {
        useVarFDorder = true;
        setFDOrder(config.get<std::string>("gridConfigurationFilename"));
    } else {
        SCAI_ASSERT(!config.get<bool>("useVariableFDoperators"), "Variable FD operators are not available for stencil matrices")
        // Set FD-order to class member
        setFDOrder(config.get<IndexType>("spatialFDorder"));
    }
    isSetup = true;
}

//! \brief Setup configuration of the derivative object
/*!
 *
 \param config Configuration
 \param FDorder std::vector with FD orders per layer
 */
template <typename ValueType>
void KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::setup(Configuration::Configuration const &config, std::vector<IndexType> &FDorder)
{
    useFreeSurface = config.get<IndexType>("FreeSurface");
    useVarGrid = config.get<bool>("useVariableGrid");

    std::string type = config.get<std::string>("equationType");
    std::transform(type.begin(), type.end(), type.begin(), ::tolower);
    if ((type.compare("sh") == 0) || (type.compare("viscosh") == 0) || (type.compare("elastic") == 0) || (type.compare("viscoelastic") == 0)) {
        isElastic = true;
    }

    DT = config.get<ValueType>("DT");

    useStencilMatrix = false;
    useHybridFreeSurface = false;

    SCAI_ASSERT(config.get<IndexType>("partitioning") != 1, "grid partition is not available for variable FDorders")

    useVarFDorder = true;
    setFDCoef();
    setFDOrder(FDorder);
    isSetup = true;
}

template <typename ValueType>
ValueType KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::getMemoryUsage(scai::dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> const &modelCoordinates, IndexType numDMatrices, IndexType numInterpMatrices)
{
    ValueType size = 0;
    ValueType mega = 1024 * 1024;
    ValueType sizeInterp = 0;
    if (useStencilMatrix) {
        size = getMemoryStencilMatrix(dist) / mega * numDMatrices;
    } else {
        size = getMemorySparseMatrix(dist, modelCoordinates) / mega * numDMatrices;
    }

    if (useVarGrid) {
        sizeInterp = getMemoryInterpolationMatrix(dist) / mega * numInterpMatrices;
    }
    return (size + sizeInterp);
}

//! \brief Calculate Dxf matrix
/*!
 *
 \param dist Distribution
 */
template <typename ValueType>
void KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::calcDxf(scai::dmemo::DistributionPtr dist)
{
    SCAI_REGION("Derivatives.calcDxfStencil")
    // Attention: keep in mind topology NZ x NY x NX

    common::Stencil1D<ValueType> stencilId(1);
    common::Stencil3D<ValueType> stencil(stencilId, stencilId, stencilFDmap[spatialFDorderVec.at(0)]);
    // use dist for distribution
    Dxf.define(dist, stencil);
}
//! \brief Calculate Dxf sparse matrix
/*!
 *
 \param modelCoordinates Coordinate class, which eg. maps 3D coordinates to 1D model indices
 \param dist Distribution
 */
template <typename ValueType>
void KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::calcDxf(Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::DistributionPtr dist)
{
    SCAI_REGION("Derivatives.calcDxfSparse")
    hmemo::HArray<IndexType> ownedIndexes; // all (global) points owned by this process
    dist->getOwnedIndexes(ownedIndexes);

    lama::MatrixAssembly<ValueType> assembly;
    //assembly.reserve(ownedIndexes.size() * 6);
    IndexType X = 0;
    IndexType columnIndex = 0;
    IndexType j = 0;

    for (IndexType ownedIndex : hmemo::hostReadAccess(ownedIndexes)) {

        Acquisition::coordinate3D coordinate = modelCoordinates.index2coordinate(ownedIndex);
        const IndexType &dhFactor = modelCoordinates.getDHFactor(coordinate);

        auto spatialFDorder = spatialFDorderVec.at(0);
        if (useVarFDorder) {
            const auto layer = modelCoordinates.getLayer(coordinate);
            spatialFDorder = spatialFDorderVec[layer];
        }

        for (j = 0; j < spatialFDorder; j++) {

            if (!modelCoordinates.locatedOnInterface(coordinate))
                X = coordinate.x + modelCoordinates.getDHFactor(coordinate) * (j - spatialFDorder / 2 + 1);
            else
                X = coordinate.x - dhFactor / 3 + modelCoordinates.getDHFactor(coordinate) * (j - spatialFDorder / 2 + 1);

            if ((X >= 0) && (X < modelCoordinates.getNX())) {
                columnIndex = modelCoordinates.coordinate2index(X, coordinate.y, coordinate.z);
                assembly.push(ownedIndex, columnIndex, stencilFDmap[spatialFDorder].values()[j] / modelCoordinates.getDH(coordinate));
            }
        }
    }

    DxfSparse = lama::zero<SparseFormat>(dist, dist);
    DxfSparse.fillFromAssembly(assembly);
}

//! \brief Calculate Dyf matrix
/*!
 *
 \param dist Distribution
 */
template <typename ValueType>
void KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::calcDyf(scai::dmemo::DistributionPtr dist)
{
    SCAI_REGION("Derivatives.calcDyfStencil")
    common::Stencil1D<ValueType> stencilId(1);
    common::Stencil3D<ValueType> stencil(stencilFDmap[spatialFDorderVec.at(0)], stencilId, stencilId);
    // use dist for distribution
    Dyf.define(dist, stencil);
}

//! \brief Calculate Dyf sparse matrix
/*!
 *
 \param modelCoordinates Coordinate class, which eg. maps 3D coordinates to 1D model indices
 \param dist Distribution
 */
template <typename ValueType>
void KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::calcDyf(Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::DistributionPtr dist)
{
    SCAI_REGION("Derivatives.calcDyfSparse")
    hmemo::HArray<IndexType> ownedIndexes; // all (global) points owned by this process
    dist->getOwnedIndexes(ownedIndexes);

    lama::MatrixAssembly<ValueType> assembly;
    //assembly.reserve(ownedIndexes.size() * 6);
    IndexType Y = 0;
    IndexType columnIndex = 0;
    IndexType j = 0;
    ValueType DH = 0;
    IndexType dhFactor = 0;
    IndexType layer = 0;

    for (IndexType ownedIndex : hmemo::hostReadAccess(ownedIndexes)) {

        Acquisition::coordinate3D coordinate = modelCoordinates.index2coordinate(ownedIndex);

        layer = modelCoordinates.getLayer(coordinate);
        auto spatialFDorder = spatialFDorderVec.at(0);
        if (useVarFDorder) {
            spatialFDorder = spatialFDorderVec[layer];
        }

        /* reduce FDorder at the variable grid interfaces*/
        if (useVarGrid) {
            auto distance = modelCoordinates.distToInterface(coordinate.y) / modelCoordinates.getDHFactor(layer);
            if (distance == 0)
                spatialFDorder = 2;
            else if (spatialFDorder > distance * 2)
                spatialFDorder = distance * 2;
        }

        //Transition from coarse (layer) to fine grid (layer+1) uses fine operator for Dyf
        if ((modelCoordinates.locatedOnInterface(coordinate)) && (modelCoordinates.getTransition(coordinate) == -1)) {
            dhFactor = modelCoordinates.getDHFactor(layer + 1);
            DH = modelCoordinates.getDH(layer + 1);
        } else {
            dhFactor = modelCoordinates.getDHFactor(layer);
            DH = modelCoordinates.getDH(layer);
        }

        for (j = 0; j < spatialFDorder; j++) {
            Y = coordinate.y + dhFactor * (j - spatialFDorder / 2 + 1);

            if ((Y >= 0) && (Y < modelCoordinates.getNY())) {
                columnIndex = modelCoordinates.coordinate2index(coordinate.x, Y, coordinate.z);
                assembly.push(ownedIndex, columnIndex, stencilFDmap[spatialFDorder].values()[j] / DH);
            }
        }
    }

    DyfSparse = lama::zero<SparseFormat>(dist, dist);
    DyfSparse.fillFromAssembly(assembly);
}

//! \brief Calculate Dzf matrix
/*!
 *
 \param dist Distribution
 */
template <typename ValueType>
void KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::calcDzf(scai::dmemo::DistributionPtr dist)
{
    SCAI_REGION("Derivatives.calcDzfStencil")
    common::Stencil1D<ValueType> stencilId(1);
    common::Stencil3D<ValueType> stencil(stencilId, stencilFDmap[spatialFDorderVec.at(0)], stencilId);
    // use dist for distribution
    Dzf.define(dist, stencil);
}

//! \brief Calculate Dzf sparse matrix
/*!
 *
\param modelCoordinates Coordinate class, which eg. maps 3D coordinates to 1D model indices
 \param dist Distribution
 */
template <typename ValueType>
void KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::calcDzf(Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::DistributionPtr dist)
{
    SCAI_REGION("Derivatives.calcDzfSparse")
    hmemo::HArray<IndexType> ownedIndexes; // all (global) points owned by this process
    dist->getOwnedIndexes(ownedIndexes);

    lama::MatrixAssembly<ValueType> assembly;
    // assembly.reserve(ownedIndexes.size() * 6);
    IndexType Z = 0;
    IndexType columnIndex = 0;
    IndexType j = 0;

    for (IndexType ownedIndex : hmemo::hostReadAccess(ownedIndexes)) {

        Acquisition::coordinate3D coordinate = modelCoordinates.index2coordinate(ownedIndex);

        const IndexType &dhFactor = modelCoordinates.getDHFactor(coordinate);
        auto spatialFDorder = spatialFDorderVec.at(0);
        if (useVarFDorder) {
            const auto layer = modelCoordinates.getLayer(coordinate);
            spatialFDorder = spatialFDorderVec[layer];
        }

        for (j = 0; j < spatialFDorder; j++) {

            if (!modelCoordinates.locatedOnInterface(coordinate))
                Z = coordinate.z + modelCoordinates.getDHFactor(coordinate) * (j - spatialFDorder / 2 + 1);
            else
                Z = coordinate.z - dhFactor / 3 + modelCoordinates.getDHFactor(coordinate) * (j - spatialFDorder / 2 + 1);

            if ((Z >= 0) && (Z < modelCoordinates.getNZ())) {
                columnIndex = modelCoordinates.coordinate2index(coordinate.x, coordinate.y, Z);
                assembly.push(ownedIndex, columnIndex, stencilFDmap[spatialFDorder].values()[j] / modelCoordinates.getDH(coordinate));
            }
        }
    }
    DzfSparse = lama::zero<SparseFormat>(dist, dist);
    DzfSparse.fillFromAssembly(assembly);
}
//! \brief Calculate DyfFreeSurface matrix
/*!
 *
 \param modelCoordinates Coordinate class, which eg. maps 3D coordinates to 1D model indices
 \param dist Distribution
 */
template <typename ValueType>
void KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::calcDyfFreeSurface(Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::DistributionPtr dist)
{
    SCAI_REGION("Derivatives.calcDyfFreeSurface")
    hmemo::HArray<IndexType> ownedIndexes; // all (global) points owned by this process
    dist->getOwnedIndexes(ownedIndexes);

    lama::MatrixAssembly<ValueType> assembly;

    const ValueType ZERO = 0;

    ValueType DH = ZERO;
    IndexType dhFactor = ZERO;

    for (IndexType ownedIndex : hmemo::hostReadAccess(ownedIndexes)) {

        Acquisition::coordinate3D coordinate = modelCoordinates.index2coordinate(ownedIndex);
        IndexType layer = modelCoordinates.getLayer(coordinate);

        //Transition from coarse (layer) to fine grid (layer+1) uses fine operator for Dyf
        if ((modelCoordinates.locatedOnInterface(coordinate)) && (modelCoordinates.getTransition(coordinate) == -1)) {
            dhFactor = modelCoordinates.getDHFactor(layer + 1);
            DH = modelCoordinates.getDH(layer + 1);
        } else {
            dhFactor = modelCoordinates.getDHFactor(layer);
            DH = modelCoordinates.getDH(layer);
        }
        auto spatialFDorder = spatialFDorderVec.at(0);
        if (useVarFDorder) {
            spatialFDorder = spatialFDorderVec[layer];
        }

        if (useVarGrid) {
            /* reduce FDorder at the variable grid interfaces*/
            auto distance = modelCoordinates.distToInterface(coordinate.y) / modelCoordinates.getDHFactor(layer);
            if (distance == 0)
                spatialFDorder = 2;
            else if (spatialFDorder > distance * 2)
                spatialFDorder = distance * 2;
        }

        for (IndexType j = 0; j < spatialFDorder; j++) {
            IndexType Y = coordinate.y + dhFactor * (j - spatialFDorder / 2 + 1);

            ValueType fdCoeff = stencilFDmap[spatialFDorder].values()[j];
            ValueType diffCoeff = ZERO;

            if (spatialFDorder >= (2 + 2 * coordinate.y / dhFactor + j)) {
                IndexType ImageIndex = spatialFDorder - 2 - 2 * coordinate.y / dhFactor - j;
                diffCoeff = stencilFDmap[spatialFDorder].values()[ImageIndex];
            }

            if ((Y >= 0) && (Y < modelCoordinates.getNY())) {
                IndexType columnIndex = modelCoordinates.coordinate2index(coordinate.x, Y, coordinate.z);
                if (!useHybridFreeSurface)
                    assembly.push(ownedIndex, columnIndex, (fdCoeff - diffCoeff) / DH); // push all coefficients
                else if (ZERO != diffCoeff)
                    assembly.push(ownedIndex, columnIndex, -diffCoeff / DH); // push only diffs to stencil matrix
            }
        }
    }

    DyfFreeSurfaceSparse = lama::zero<SparseFormat>(dist, dist);
    DyfFreeSurfaceSparse.fillFromAssembly(assembly);

    if (useHybridFreeSurface) {
        // define the stencil matrix for hybrid matrix
        // ToDo: why not simply use the stencil matrix Dyb
        common::Stencil1D<ValueType> stencilId(1);
        common::Stencil3D<ValueType> stencil(stencilFDmap[spatialFDorderVec.at(0)], stencilId, stencilId);
        DyfFreeSurfaceStencil.define(dist, stencil);
        DyfFreeSurfaceStencil *= 1 / modelCoordinates.getDH();
    }
}

//! \brief Calculate DybFreeSurface matrix
/*!
  \param modelCoordinates Coordinate class, which eg. maps 3D coordinates to 1D model indices
 \param dist Distribution
 */
template <typename ValueType>
void KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::calcDybFreeSurface(Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::DistributionPtr dist)
{
    SCAI_REGION("Derivatives.calcDybFreeSurface")

    hmemo::HArray<IndexType> ownedIndexes; // all (global) points owned by this process
    dist->getOwnedIndexes(ownedIndexes);

    lama::MatrixAssembly<ValueType> assembly;

    const ValueType ZERO = 0;

    for (IndexType ownedIndex : hmemo::hostReadAccess(ownedIndexes)) {

        Acquisition::coordinate3D coordinate = modelCoordinates.index2coordinate(ownedIndex);
        const IndexType &layer = modelCoordinates.getLayer(coordinate);
        const IndexType &dhFactor = modelCoordinates.getDHFactor(coordinate);

        auto spatialFDorder = spatialFDorderVec.at(0);
        if (useVarFDorder) {
            spatialFDorder = spatialFDorderVec[layer];
        }

        if (useVarGrid) {
            auto distance = modelCoordinates.distToInterface(coordinate.y) / modelCoordinates.getDHFactor(layer);

            /* reduce FDorder at the variable grid interfaces*/
            if (distance == 0)
                spatialFDorder = 2;
            else if (spatialFDorder > distance * 2)
                spatialFDorder = distance * 2;
        }

        for (IndexType j = 0; j < spatialFDorder; j++) {
            IndexType Y = coordinate.y + dhFactor * (j - spatialFDorder / 2);

            // apply coordinate correction in the fine staggered grid (staggered in y-direction, coordinates are only correct for full grid points)
            if (modelCoordinates.locatedOnInterface(coordinate)) {
                // fG<->cG transotion=1 -> layer-1, cG<->fG transotion=-1  -> layer+1
                if ((j == 0) && (modelCoordinates.getTransition(coordinate) == 1))
                    Y += modelCoordinates.getDHFactor(layer - 1);
                if ((j == 1) && (modelCoordinates.getTransition(coordinate) == -1))
                    Y += modelCoordinates.getDHFactor(layer + 1);
            }

            ValueType fdCoeff = stencilFDmap[spatialFDorder].values()[j];
            ValueType diffCoeff = ZERO;

            if (spatialFDorder >= (1 + 2 * coordinate.y / dhFactor + j)) {
                IndexType ImageIndex = spatialFDorder - 1 - 2 * coordinate.y / dhFactor - j;
                diffCoeff = stencilFDmap[spatialFDorder].values()[ImageIndex];
            }

            if ((Y >= 0) && (Y < modelCoordinates.getNY())) {
                IndexType columnIndex = modelCoordinates.coordinate2index(coordinate.x, Y, coordinate.z); // push all coefficients
                if (!useHybridFreeSurface)
                    assembly.push(ownedIndex, columnIndex, (fdCoeff - diffCoeff) / modelCoordinates.getDH(coordinate));
                else if (ZERO != diffCoeff)
                    assembly.push(ownedIndex, columnIndex, -diffCoeff / modelCoordinates.getDH(coordinate)); // push only diffs to stencil matrix
            }
        }
    }

    DybFreeSurfaceSparse = lama::zero<SparseFormat>(dist, dist);
    DybFreeSurfaceSparse.fillFromAssembly(assembly);

    if (useHybridFreeSurface) {
        // define the stencil matrix for hybrid matrix
        // ToDo: why not simply use the stencil matrix Dyb

        common::Stencil1D<ValueType> stencilId(1);
        common::Stencil1D<ValueType> stencilBD;
        stencilBD.transpose(stencilFDmap[spatialFDorderVec.at(0)]);
        stencilBD.scale(-1);
        common::Stencil3D<ValueType> stencil(stencilBD, stencilId, stencilId);

        DybFreeSurfaceStencil.define(dist, stencil);
        DybFreeSurfaceStencil *= 1 / modelCoordinates.getDH();
    }
}

//! \brief Calculate DybFreeSurface matrix
/*!
  \param modelCoordinates Coordinate class, which eg. maps 3D coordinates to 1D model indices
 \param dist Distribution
 */
template <typename ValueType>
void KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::calcDybStaggeredXFreeSurface(Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::DistributionPtr dist)
{
    SCAI_REGION("Derivatives.calcDybStaggeredXFreeSurface")

    hmemo::HArray<IndexType> ownedIndexes; // all (global) points owned by this process
    dist->getOwnedIndexes(ownedIndexes);

    lama::MatrixAssembly<ValueType> assembly;
    // assembly.reserve(ownedIndexes.size() * 6);
    IndexType Y = 0;
    ValueType fdCoeff = 0;
    IndexType ImageIndex = 0;
    IndexType columnIndex = 0;
    IndexType j = 0;

    for (IndexType ownedIndex : hmemo::hostReadAccess(ownedIndexes)) {

        Acquisition::coordinate3D coordinate = modelCoordinates.index2coordinate(ownedIndex);
        const IndexType &layer = modelCoordinates.getLayer(coordinate);
        const IndexType &dhFactor = modelCoordinates.getDHFactor(coordinate);
        auto spatialFDorder = spatialFDorderVec.at(0);
        if (useVarFDorder) {
            spatialFDorder = spatialFDorderVec[layer];
        }

        if (useVarGrid) {
            auto distance = modelCoordinates.distToInterface(coordinate.y) / modelCoordinates.getDHFactor(layer);

            /* reduce FDorder at the variable grid interfaces*/
            if (distance == 0)
                spatialFDorder = 2;
            else if (spatialFDorder > distance * 2)
                spatialFDorder = distance * 2;
        }

        std::vector<IndexType> X(spatialFDorder, coordinate.x);

        if ((modelCoordinates.locatedOnInterface(coordinate.y - (spatialFDorder / 2 * dhFactor))) && (modelCoordinates.getTransition(coordinate.y - (spatialFDorder / 2 * dhFactor)) == 1)) {
            if (X[1] < modelCoordinates.getNX() - dhFactor / 3) {
                X[0] += dhFactor / 3;
            }
        }

        if ((modelCoordinates.locatedOnInterface(coordinate)) && (modelCoordinates.getTransition(coordinate) == -1)) {
            if (X[1] >= dhFactor / 3) {
                X[0] -= dhFactor / 3;
            }
        }

        for (j = 0; j < spatialFDorder; j++) {
            Y = coordinate.y + dhFactor * (j - spatialFDorder / 2);

            // apply coordinate correction in the fine staggered grid (staggered in y-direction, coordinates are only correct for full grid points)
            if (modelCoordinates.locatedOnInterface(coordinate)) {
                // fG<->cG transotion=1 -> layer-1, cG<->fG transotion=-1  -> layer+1
                if ((j == 0) && (modelCoordinates.getTransition(coordinate) == 1))
                    Y += modelCoordinates.getDHFactor(layer - 1);
                if ((j == 1) && (modelCoordinates.getTransition(coordinate) == -1))
                    Y += modelCoordinates.getDHFactor(layer + 1);
            }

            fdCoeff = stencilFDmap[spatialFDorder].values()[j];

            if (spatialFDorder >= (1 + 2 * coordinate.y / dhFactor + j)) {
                ImageIndex = spatialFDorder - 1 - 2 * coordinate.y / dhFactor - j;
                fdCoeff -= stencilFDmap[spatialFDorder].values()[ImageIndex];
            }

            if ((Y >= 0) && (Y < modelCoordinates.getNY())) {
                columnIndex = modelCoordinates.coordinate2index(X[j], Y, coordinate.z);
                assembly.push(ownedIndex, columnIndex, fdCoeff / modelCoordinates.getDH(coordinate));
            }
        }
    }

    DybStaggeredXFreeSurface = lama::zero<SparseFormat>(dist, dist);
    DybStaggeredXFreeSurface.fillFromAssembly(assembly);
}

//! \brief Calculate DybFreeSurface matrix
/*!
  \param modelCoordinates Coordinate class, which eg. maps 3D coordinates to 1D model indices
 \param dist Distribution
 */
template <typename ValueType>
void KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::calcDybStaggeredZFreeSurface(Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::DistributionPtr dist)
{
    SCAI_REGION("Derivatives.calcDybStaggeredZFreeSurface")

    hmemo::HArray<IndexType> ownedIndexes; // all (global) points owned by this process
    dist->getOwnedIndexes(ownedIndexes);

    lama::MatrixAssembly<ValueType> assembly;
    // assembly.reserve(ownedIndexes.size() * 6);
    IndexType Y = 0;
    ValueType fdCoeff = 0;
    IndexType ImageIndex = 0;
    IndexType columnIndex = 0;
    IndexType j = 0;

    for (IndexType ownedIndex : hmemo::hostReadAccess(ownedIndexes)) {

        Acquisition::coordinate3D coordinate = modelCoordinates.index2coordinate(ownedIndex);
        const IndexType &layer = modelCoordinates.getLayer(coordinate);
        const IndexType &dhFactor = modelCoordinates.getDHFactor(coordinate);

        auto spatialFDorder = spatialFDorderVec.at(0);
        if (useVarFDorder) {
            spatialFDorder = spatialFDorderVec[layer];
        }

        if (useVarGrid) {
            auto distance = modelCoordinates.distToInterface(coordinate.y) / modelCoordinates.getDHFactor(layer);

            /* reduce FDorder at the variable grid interfaces*/
            if (distance == 0)
                spatialFDorder = 2;
            else if (spatialFDorder > distance * 2)
                spatialFDorder = distance * 2;
        }

        std::vector<IndexType> Z(spatialFDorder, coordinate.z);

        if ((modelCoordinates.locatedOnInterface(coordinate.y - (spatialFDorder / 2 * dhFactor))) && (modelCoordinates.getTransition(coordinate.y - (spatialFDorder / 2 * dhFactor)) == 1)) {
            if (Z[1] < modelCoordinates.getNZ() - dhFactor / 3) {
                Z[0] += dhFactor / 3;
            }
        }

        if ((modelCoordinates.locatedOnInterface(coordinate)) && (modelCoordinates.getTransition(coordinate) == -1)) {
            if (Z[1] >= dhFactor / 3) {
                Z[0] -= dhFactor / 3;
            }
        }

        for (j = 0; j < spatialFDorder; j++) {
            Y = coordinate.y + dhFactor * (j - spatialFDorder / 2);

            // apply coordinate correction in the fine staggered grid (staggered in y-direction, coordinates are only correct for full grid points)
            if (modelCoordinates.locatedOnInterface(coordinate)) {
                // fG<->cG transotion=1 -> layer-1, cG<->fG transotion=-1  -> layer+1
                if ((j == 0) && (modelCoordinates.getTransition(coordinate) == 1))
                    Y += modelCoordinates.getDHFactor(layer - 1);
                if ((j == 1) && (modelCoordinates.getTransition(coordinate) == -1))
                    Y += modelCoordinates.getDHFactor(layer + 1);
            }

            fdCoeff = stencilFDmap[spatialFDorder].values()[j];

            if (spatialFDorder >= (1 + 2 * coordinate.y / dhFactor + j)) {
                ImageIndex = spatialFDorder - 1 - 2 * coordinate.y / dhFactor - j;
                fdCoeff -= stencilFDmap[spatialFDorder].values()[ImageIndex];
            }

            if ((Y >= 0) && (Y < modelCoordinates.getNY())) {
                columnIndex = modelCoordinates.coordinate2index(coordinate.x, Y, Z[j]);
                assembly.push(ownedIndex, columnIndex, fdCoeff / modelCoordinates.getDH(coordinate));
            }
        }
    }

    DybStaggeredZFreeSurface = lama::zero<SparseFormat>(dist, dist);
    DybStaggeredZFreeSurface.fillFromAssembly(assembly);
}

//! \brief Calculate Dxb sparse matrix
/*!
 *
 \param modelCoordinates Coordinate class, which eg. maps 3D coordinates to 1D model indices
 \param dist Distribution
 */
template <typename ValueType>
void KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::calcDxb(Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::DistributionPtr dist)
{
    SCAI_REGION("Derivatives.calcDxbSparse")
    hmemo::HArray<IndexType> ownedIndexes; // all (global) points owned by this process
    dist->getOwnedIndexes(ownedIndexes);

    lama::MatrixAssembly<ValueType> assembly;
    //assembly.reserve(ownedIndexes.size() * 6);
    IndexType X = 0;
    IndexType columnIndex = 0;
    IndexType j = 0;

    for (IndexType ownedIndex : hmemo::hostReadAccess(ownedIndexes)) {

        Acquisition::coordinate3D const coordinate = modelCoordinates.index2coordinate(ownedIndex);
        const IndexType &dhFactor = modelCoordinates.getDHFactor(coordinate);

        auto spatialFDorder = spatialFDorderVec.at(0);
        if (useVarFDorder) {
            const auto layer = modelCoordinates.getLayer(coordinate);
            spatialFDorder = spatialFDorderVec[layer];
        }

        for (j = 0; j < spatialFDorder; j++) {

            if (!modelCoordinates.locatedOnInterface(coordinate)) {
                X = coordinate.x + modelCoordinates.getDHFactor(coordinate) * (j - spatialFDorder / 2);
            } else {
                X = coordinate.x + (dhFactor / 3) + modelCoordinates.getDHFactor(coordinate) * (j - spatialFDorder / 2);
            }
            if ((X >= 0) && (X < modelCoordinates.getNX())) {
                columnIndex = modelCoordinates.coordinate2index(X, coordinate.y, coordinate.z);
                assembly.push(ownedIndex, columnIndex, stencilFDmap[spatialFDorder].values()[j] / modelCoordinates.getDH(coordinate));
            }
        }
    }

    DxbSparse = lama::zero<SparseFormat>(dist, dist);
    DxbSparse.fillFromAssembly(assembly);
}

//! \brief Calculate Dyb matrix
/*!
  \param modelCoordinates Coordinate class, which eg. maps 3D coordinates to 1D model indices
 \param dist Distribution
 */
template <typename ValueType>
void KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::calcDyb(Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::DistributionPtr dist)
{
    SCAI_REGION("Derivatives.calcDybSparse")

    hmemo::HArray<IndexType> ownedIndexes; // all (global) points owned by this process
    dist->getOwnedIndexes(ownedIndexes);

    lama::MatrixAssembly<ValueType> assembly;
    //assembly.reserve(ownedIndexes.size() * 6);
    IndexType Y = 0;
    IndexType columnIndex = 0;
    IndexType j = 0;

    for (IndexType ownedIndex : hmemo::hostReadAccess(ownedIndexes)) {

        Acquisition::coordinate3D coordinate = modelCoordinates.index2coordinate(ownedIndex);

        const IndexType &layer = modelCoordinates.getLayer(coordinate);
        const IndexType &dhFactor = modelCoordinates.getDHFactor(coordinate);

        auto spatialFDorder = spatialFDorderVec[0];
        if (useVarFDorder) {
            spatialFDorder = spatialFDorderVec[layer];
        }

        if (useVarGrid) {
            auto distance = modelCoordinates.distToInterface(coordinate.y) / modelCoordinates.getDHFactor(layer);

            /* reduce FDorder at the variable grid interfaces*/
            if (distance == 0)
                spatialFDorder = 2;
            else if (spatialFDorder > distance * 2)
                spatialFDorder = distance * 2;
        }

        for (j = 0; j < spatialFDorder; j++) {

            Y = coordinate.y + dhFactor * (j - spatialFDorder / 2);

            // apply coordinate correction in the fine staggered grid (staggered in y-direction, coordinates are only correct for full grid points)
            if (modelCoordinates.locatedOnInterface(coordinate)) {
                // fG<->cG transotion=1 -> layer-1, cG<->fG transotion=-1  -> layer+1
                if ((j == 0) && (modelCoordinates.getTransition(coordinate) == 1))
                    Y += modelCoordinates.getDHFactor(layer - 1);
                if ((j == 1) && (modelCoordinates.getTransition(coordinate) == -1))
                    Y += modelCoordinates.getDHFactor(layer + 1);
            }

            if ((Y >= 0) && (Y < modelCoordinates.getNY())) {
                columnIndex = modelCoordinates.coordinate2index(coordinate.x, Y, coordinate.z);
                assembly.push(ownedIndex, columnIndex, stencilFDmap[spatialFDorder].values()[j] / modelCoordinates.getDH(coordinate));
            }
        }
    }

    DybSparse = lama::zero<SparseFormat>(dist, dist);
    DybSparse.fillFromAssembly(assembly);
}

//! \brief Calculate Dyf sparse matrix
/*!
 *
 \param modelCoordinates Coordinate class, which eg. maps 3D coordinates to 1D model indices
 \param dist Distribution
 */
template <typename ValueType>
void KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::calcDyfStaggeredX(Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::DistributionPtr dist)
{
    SCAI_REGION("Derivatives.calcDyfStaggeredX")

    hmemo::HArray<IndexType> ownedIndexes; // all (global) points owned by this process
    dist->getOwnedIndexes(ownedIndexes);

    lama::MatrixAssembly<ValueType> assembly;
    //assembly.reserve(ownedIndexes.size() * 6);

    IndexType Y = 0;
    IndexType columnIndex = 0;
    IndexType j = 0;
    ValueType DH = 0;
    IndexType dhFactor = 0;
    IndexType layer = 0;

    for (IndexType ownedIndex : hmemo::hostReadAccess(ownedIndexes)) {

        Acquisition::coordinate3D coordinate = modelCoordinates.index2coordinate(ownedIndex);
        layer = modelCoordinates.getLayer(coordinate);
        auto spatialFDorder = spatialFDorderVec[0];
        if (useVarFDorder) {
            spatialFDorder = spatialFDorderVec[layer];
        }
        /* reduce FDorder at the variable grid interfaces*/
        auto distance = modelCoordinates.distToInterface(coordinate.y) / modelCoordinates.getDHFactor(layer);
        if (distance == 0)
            spatialFDorder = 2;
        else if (spatialFDorder > distance * 2)
            spatialFDorder = distance * 2;

        std::vector<IndexType> X(spatialFDorder, coordinate.x);

        //Transition from coarse (layer) to fine grid (layer+1) uses fine operator for Dyf
        if ((modelCoordinates.locatedOnInterface(coordinate)) && (modelCoordinates.getTransition(coordinate) == -1)) {
            dhFactor = modelCoordinates.getDHFactor(layer + 1);
            DH = modelCoordinates.getDH(layer + 1);
        } else {
            dhFactor = modelCoordinates.getDHFactor(layer);
            DH = modelCoordinates.getDH(layer);
        }

        if ((modelCoordinates.locatedOnInterface(coordinate)) && (modelCoordinates.getTransition(coordinate) == 1)) {
            if (X[0] >= dhFactor / 3) {
                X[1] -= dhFactor / 3;
            }
            /*corresponding staggered point in coarse layer has a different X coordinate
              0  0  1  1  2  2  3  3
              x  ^  x  ^-|x  ^  x  ^
                         |
                         |
              v        x |      v  
                         |          
                         |       
              x        ^-|      x        ^
              0        0        3        3
              */
        }

        if ((modelCoordinates.locatedOnInterface(coordinate.y + (spatialFDorder / 2 * dhFactor))) && (modelCoordinates.getTransition(coordinate.y + (spatialFDorder / 2 * dhFactor)) == -1)) {
            if (X[0] < modelCoordinates.getNX() - dhFactor / 3) {
                X[spatialFDorder - 1] += dhFactor / 3;
            }
        }

        for (j = 0; j < spatialFDorder; j++) {
            Y = coordinate.y + dhFactor * (j - spatialFDorder / 2 + 1);

            if ((Y >= 0) && (Y < modelCoordinates.getNY())) {
                columnIndex = modelCoordinates.coordinate2index(X[j], Y, coordinate.z);
                assembly.push(ownedIndex, columnIndex, stencilFDmap[spatialFDorder].values()[j] / DH);
            }
        }
    }

    DyfStaggeredXSparse = lama::zero<SparseFormat>(dist, dist);
    DyfStaggeredXSparse.fillFromAssembly(assembly);
}

//! \brief Calculate DybStaggeredX matrix
/*!
  \param modelCoordinates Coordinate class, which eg. maps 3D coordinates to 1D model indices
 \param dist Distribution
 */
template <typename ValueType>
void KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::calcDybStaggeredX(Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::DistributionPtr dist)
{
    SCAI_REGION("Derivatives.calcDybStaggeredX")

    hmemo::HArray<IndexType> ownedIndexes; // all (global) points owned by this process
    dist->getOwnedIndexes(ownedIndexes);

    lama::MatrixAssembly<ValueType> assembly;
    //assembly.reserve(ownedIndexes.size() * 6);
    IndexType Y = 0;
    IndexType columnIndex = 0;
    IndexType j = 0;

    for (IndexType ownedIndex : hmemo::hostReadAccess(ownedIndexes)) {

        Acquisition::coordinate3D const coordinate = modelCoordinates.index2coordinate(ownedIndex);

        const IndexType &layer = modelCoordinates.getLayer(coordinate);
        const IndexType &dhFactor = modelCoordinates.getDHFactor(coordinate);
        auto distance = modelCoordinates.distToInterface(coordinate.y) / modelCoordinates.getDHFactor(layer);
        auto spatialFDorder = spatialFDorderVec[0];
        if (useVarFDorder) {
            spatialFDorder = spatialFDorderVec[layer];
        }
        /* reduce FDorder at the variable grid interfaces*/
        if (distance == 0)
            spatialFDorder = 2;
        else if (spatialFDorder > distance * 2)
            spatialFDorder = distance * 2;

        std::vector<IndexType> X(spatialFDorder, coordinate.x);

        if ((modelCoordinates.locatedOnInterface(coordinate.y - (spatialFDorder / 2 * dhFactor))) && (modelCoordinates.getTransition(coordinate.y - (spatialFDorder / 2 * dhFactor)) == 1)) {
            if (X[1] < modelCoordinates.getNX() - dhFactor / 3) {
                X[0] += dhFactor / 3;
            }
        }

        if ((modelCoordinates.locatedOnInterface(coordinate)) && (modelCoordinates.getTransition(coordinate) == -1)) {
            if (X[1] >= dhFactor / 3) {
                X[0] -= dhFactor / 3;
            }
        }

        for (j = 0; j < spatialFDorder; j++) {

            Y = coordinate.y + dhFactor * (j - spatialFDorder / 2);

            // apply coordinate correction in the fine staggered grid (staggered in y-direction, coordinates are only correct for full grid points)
            if (modelCoordinates.locatedOnInterface(coordinate)) {
                // fG<->cG transotion=1 -> layer-1, cG<->fG transotion=-1  -> layer+1
                if ((j == 0) && (modelCoordinates.getTransition(coordinate) == 1))
                    Y += modelCoordinates.getDHFactor(layer - 1);
                if ((j == 1) && (modelCoordinates.getTransition(coordinate) == -1))
                    Y += modelCoordinates.getDHFactor(layer + 1);
            }

            if ((Y >= 0) && (Y < modelCoordinates.getNY())) {
                columnIndex = modelCoordinates.coordinate2index(X[j], Y, coordinate.z);
                assembly.push(ownedIndex, columnIndex, stencilFDmap[spatialFDorder].values()[j] / modelCoordinates.getDH(coordinate));
            }
        }
    }

    DybStaggeredXSparse = lama::zero<SparseFormat>(dist, dist);
    DybStaggeredXSparse.fillFromAssembly(assembly);
}

//! \brief Calculate Dyf sparse matrix
/*!
 *
 \param modelCoordinates Coordinate class, which eg. maps 3D coordinates to 1D model indices
 \param dist Distribution
 */
template <typename ValueType>
void KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::calcDyfStaggeredZ(Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::DistributionPtr dist)
{
    SCAI_REGION("Derivatives.calcDyfStaggeredZ")

    hmemo::HArray<IndexType> ownedIndexes; // all (global) points owned by this process
    dist->getOwnedIndexes(ownedIndexes);

    lama::MatrixAssembly<ValueType> assembly;
    //assembly.reserve(ownedIndexes.size() * 6);

    IndexType Y = 0;
    IndexType columnIndex = 0;
    IndexType j = 0;
    ValueType DH = 0;
    IndexType dhFactor = 0;
    IndexType layer = 0;

    for (IndexType ownedIndex : hmemo::hostReadAccess(ownedIndexes)) {

        Acquisition::coordinate3D coordinate = modelCoordinates.index2coordinate(ownedIndex);
        layer = modelCoordinates.getLayer(coordinate);

        auto spatialFDorder = spatialFDorderVec[0];
        if (useVarFDorder) {
            spatialFDorder = spatialFDorderVec[layer];
        }
        /* reduce FDorder at the variable grid interfaces*/
        auto distance = modelCoordinates.distToInterface(coordinate.y) / modelCoordinates.getDHFactor(layer);
        if (distance == 0)
            spatialFDorder = 2;
        else if (spatialFDorder > distance * 2)
            spatialFDorder = distance * 2;

        std::vector<IndexType> Z(spatialFDorder, coordinate.z);

        //Transition from coarse (layer) to fine grid (layer+1) uses fine operator for Dyf
        if ((modelCoordinates.locatedOnInterface(coordinate)) && (modelCoordinates.getTransition(coordinate) == -1)) {
            dhFactor = modelCoordinates.getDHFactor(layer + 1);
            DH = modelCoordinates.getDH(layer + 1);
        } else {
            dhFactor = modelCoordinates.getDHFactor(layer);
            DH = modelCoordinates.getDH(layer);
        }

        if ((modelCoordinates.locatedOnInterface(coordinate)) && (modelCoordinates.getTransition(coordinate) == 1)) {
            if (Z[0] >= dhFactor / 3) {
                Z[1] -= dhFactor / 3;
            }
            /*corresponding staggered point in coarse layer has a different X coordinate
              0  0  1  1  2  2  3  3
              x  ^  x  ^-|x  ^  x  ^
                         |
                         |
              v        x |      v  
                         |          
                         |       
              x        ^-|      x        ^
              0        0        3        3
              */
        }

        if ((modelCoordinates.locatedOnInterface(coordinate.y + (spatialFDorder / 2 * dhFactor))) && (modelCoordinates.getTransition(coordinate.y + (spatialFDorder / 2 * dhFactor)) == -1)) {
            if (Z[0] < modelCoordinates.getNZ() - dhFactor / 3) {
                Z[spatialFDorder - 1] += dhFactor / 3;
            }
        }

        for (j = 0; j < spatialFDorder; j++) {
            Y = coordinate.y + dhFactor * (j - spatialFDorder / 2 + 1);

            if ((Y >= 0) && (Y < modelCoordinates.getNY())) {
                columnIndex = modelCoordinates.coordinate2index(coordinate.x, Y, Z[j]);
                assembly.push(ownedIndex, columnIndex, stencilFDmap[spatialFDorder].values()[j] / DH);
            }
        }
    }

    DyfStaggeredZSparse = lama::zero<SparseFormat>(dist, dist);
    DyfStaggeredZSparse.fillFromAssembly(assembly);
}

//! \brief Calculate DybStaggeredX matrix
/*!
  \param modelCoordinates Coordinate class, which eg. maps 3D coordinates to 1D model indices
 \param dist Distribution
 */
template <typename ValueType>
void KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::calcDybStaggeredZ(Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::DistributionPtr dist)
{
    SCAI_REGION("Derivatives.calcDybStaggeredZ")

    hmemo::HArray<IndexType> ownedIndexes; // all (global) points owned by this process
    dist->getOwnedIndexes(ownedIndexes);

    lama::MatrixAssembly<ValueType> assembly;
    //assembly.reserve(ownedIndexes.size() * 6);
    IndexType Y = 0;
    IndexType columnIndex = 0;
    IndexType j = 0;

    for (IndexType ownedIndex : hmemo::hostReadAccess(ownedIndexes)) {

        Acquisition::coordinate3D coordinate = modelCoordinates.index2coordinate(ownedIndex);

        const IndexType &layer = modelCoordinates.getLayer(coordinate);
        const IndexType &dhFactor = modelCoordinates.getDHFactor(coordinate);

        auto distance = modelCoordinates.distToInterface(coordinate.y) / modelCoordinates.getDHFactor(layer);

        auto spatialFDorder = spatialFDorderVec[0];
        if (useVarFDorder) {
            spatialFDorder = spatialFDorderVec[layer];
        }
        /* reduce FDorder at the variable grid interfaces*/
        if (distance == 0)
            spatialFDorder = 2;
        else if (spatialFDorder > distance * 2)
            spatialFDorder = distance * 2;

        std::vector<IndexType> Z(spatialFDorder, coordinate.z);

        if ((modelCoordinates.locatedOnInterface(coordinate.y - (spatialFDorder / 2 * dhFactor))) && (modelCoordinates.getTransition(coordinate.y - (spatialFDorder / 2 * dhFactor)) == 1)) {
            if (Z[1] < modelCoordinates.getNZ() - dhFactor / 3) {
                Z[0] += dhFactor / 3;
            }
        }

        if ((modelCoordinates.locatedOnInterface(coordinate)) && (modelCoordinates.getTransition(coordinate) == -1)) {
            if (Z[1] >= dhFactor / 3) {
                Z[0] -= dhFactor / 3;
            }
        }

        for (j = 0; j < spatialFDorder; j++) {

            Y = coordinate.y + dhFactor * (j - spatialFDorder / 2);

            // apply coordinate correction in the fine staggered grid (staggered in y-direction, coordinates are only correct for full grid points)
            if (modelCoordinates.locatedOnInterface(coordinate)) {
                // fG<->cG transotion=1 -> layer-1, cG<->fG transotion=-1  -> layer+1
                if ((j == 0) && (modelCoordinates.getTransition(coordinate) == 1))
                    Y += modelCoordinates.getDHFactor(layer - 1);
                if ((j == 1) && (modelCoordinates.getTransition(coordinate) == -1))
                    Y += modelCoordinates.getDHFactor(layer + 1);
            }

            if ((Y >= 0) && (Y < modelCoordinates.getNY())) {
                columnIndex = modelCoordinates.coordinate2index(coordinate.x, Y, Z[j]);
                assembly.push(ownedIndex, columnIndex, stencilFDmap[spatialFDorder].values()[j] / modelCoordinates.getDH(coordinate));
            }
        }
    }

    DybStaggeredZSparse = lama::zero<SparseFormat>(dist, dist);
    DybStaggeredZSparse.fillFromAssembly(assembly);
}

//! \brief Calculate Dzb sparse matrix
/*!
 *
\param modelCoordinates Coordinate class, which eg. maps 3D coordinates to 1D model indices
 \param dist Distribution
 */
template <typename ValueType>
void KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::calcDzb(Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::DistributionPtr dist)
{
    SCAI_REGION("Derivatives.calcDzbSparse")

    hmemo::HArray<IndexType> ownedIndexes; // all (global) points owned by this process
    dist->getOwnedIndexes(ownedIndexes);

    lama::MatrixAssembly<ValueType> assembly;
    //assembly.reserve(ownedIndexes.size() * 6);
    IndexType Z = 0;
    IndexType columnIndex = 0;
    IndexType j = 0;

    for (IndexType ownedIndex : hmemo::hostReadAccess(ownedIndexes)) {

        Acquisition::coordinate3D coordinate = modelCoordinates.index2coordinate(ownedIndex);

        const IndexType &dhFactor = modelCoordinates.getDHFactor(coordinate);

        auto spatialFDorder = spatialFDorderVec[0];
        if (useVarFDorder) {
            const auto layer = modelCoordinates.getLayer(coordinate);
            spatialFDorder = spatialFDorderVec[layer];
        }

        for (j = 0; j < spatialFDorder; j++) {

            if (!modelCoordinates.locatedOnInterface(coordinate))
                Z = coordinate.z + modelCoordinates.getDHFactor(coordinate) * (j - spatialFDorder / 2);
            else
                Z = coordinate.z + dhFactor / 3 + modelCoordinates.getDHFactor(coordinate) * (j - spatialFDorder / 2);

            if ((Z >= 0) && (Z < modelCoordinates.getNZ())) {
                columnIndex = modelCoordinates.coordinate2index(coordinate.x, coordinate.y, Z);
                assembly.push(ownedIndex, columnIndex, stencilFDmap[spatialFDorder].values()[j] / modelCoordinates.getDH(coordinate));
            }
        }
    }
    DzbSparse = lama::zero<SparseFormat>(dist, dist);
    DzbSparse.fillFromAssembly(assembly);
}

//! \brief Calculate interpolation Matrix acoustic (2D/3D) variable grid simulations
/*!
 * Bilinear interpolation is used for the interpolation
 * 
 *   21 ---o---o---22
 *    |            |
 *    o    o   o   o
 *    |            |
 * z  o    o   o   o
 * ^  |            |
 * | 11 ---o---o---21
 *   --> x
\param modelCoordinates Coordinate class, which eg. maps 3D coordinates to 1D model indices
 \param dist Distribution
 */
template <typename ValueType>
void KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::calcInterpolationFull(Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::DistributionPtr dist)
{
    SCAI_REGION("Derivatives.calcInterpolationFull")

    hmemo::HArray<IndexType> ownedIndexes; // all (global) points owned by this process

    dist->getOwnedIndexes(ownedIndexes);

    lama::MatrixAssembly<ValueType> assembly;
    //assembly.reserve(ownedIndexes.size() * 2);

    ValueType denom = 0;
    IndexType dhFactorFineGrid = 0;
    IndexType modx = 0;
    IndexType modz = 0;
    ValueType value;
    IndexType index;

    for (IndexType ownedIndex : hmemo::hostReadAccess(ownedIndexes)) {
        Acquisition::coordinate3D coordinate = modelCoordinates.index2coordinate(ownedIndex);

        const IndexType &x = coordinate.x;
        const IndexType &y = coordinate.y;
        const IndexType &z = coordinate.z;

        const IndexType &NX = modelCoordinates.getNX();
        const IndexType &NZ = modelCoordinates.getNZ();

        const IndexType &layer = modelCoordinates.getLayer(coordinate);
        const IndexType &dhFactor = modelCoordinates.getDHFactor(layer);
        dhFactorFineGrid = dhFactor;
        denom = ValueType(1) / ValueType(dhFactor * dhFactor);

        if (modelCoordinates.locatedOnInterface(coordinate)) {
            if (modelCoordinates.getTransition(coordinate) == 1)
                dhFactorFineGrid = modelCoordinates.getDHFactor(layer - 1);
            else if (modelCoordinates.getTransition(coordinate) == -1)
                dhFactorFineGrid = modelCoordinates.getDHFactor(layer + 1);
        }

        if (!modelCoordinates.locatedOnInterface(coordinate)) {
            assembly.push(ownedIndex, ownedIndex, ValueType(1));
        } else {
            modx = x % dhFactor;
            modz = z % dhFactor;

            // Point 11
            value = (dhFactor - modx) * (dhFactor - modz) * denom;
            index = modelCoordinates.coordinate2index(x - modx, y, z - modz);
            assembly.push(ownedIndex, index, value);

            // Point 12
            if (x + 2 * dhFactorFineGrid < NX) {
                value = modx * (dhFactor - modz) * denom;
                index = modelCoordinates.coordinate2index(x + dhFactor - modx, y, z - modz);
                assembly.push(ownedIndex, index, value);
            }

            //Point 21
            if (z + 2 * dhFactorFineGrid < NZ) {
                value = (dhFactor - modx) * modz * denom;
                try {
                    index = modelCoordinates.coordinate2index(x - modx, y, z + dhFactor - modz);
                } catch (scai::common::Exception &e) {
                    COMMON_THROWEXCEPTION("Error: coordinate2index called from calcInterpolationFull in derivatives.cpp" << e.what());
                }

                assembly.push(ownedIndex, index, value);
            }

            //Point 22
            if ((x + 2 * dhFactorFineGrid < NX) && (z + 2 * dhFactorFineGrid < NZ)) {
                value = modx * modz * denom;
                try {
                    index = modelCoordinates.coordinate2index(x + dhFactor - modx, y, z + dhFactor - modz);
                } catch (scai::common::Exception &e) {
                    COMMON_THROWEXCEPTION("Error: coordinate2index called from calcInterpolationFull in derivatives.cpp" << e.what());
                }
                assembly.push(ownedIndex, index, value);
            }
        }
    }
    // auto colDist = std::make_shared<dmemo::NoDistribution( dist.getGlobalSize() );
    InterpolationFull = lama::zero<SparseFormat>(dist, dist);

    InterpolationFull.fillFromAssembly(assembly);
}

//! \brief Calculate interpolation Matrix acoustic (2D/3D) variable grid simulations
/*!
 * Bilinear interpolation is used for the interpolation
 * 
 *   21 ---o---o---22
 *    |            |
 *    o    o   o   o
 *    |            |
 * z  o    o   o   o
 * ^  |            |
 * | 11 ---o---o---21
 *   --> x
\param modelCoordinates Coordinate class, which eg. maps 3D coordinates to 1D model indices
 \param dist Distribution
 */
template <typename ValueType>
void KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::calcInterpolationStaggeredX(Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::DistributionPtr dist)
{
    SCAI_REGION("Derivatives.calcInterpolationStaggeredX")

    hmemo::HArray<IndexType> ownedIndexes; // all (global) points owned by this process

    dist->getOwnedIndexes(ownedIndexes);

    lama::MatrixAssembly<ValueType> assembly;
    // assembly.reserve(ownedIndexes.size() * 2);

    ValueType denom = 0;
    IndexType dhFactorFineGrid = 0;
    IndexType modx = 0;
    IndexType modz = 0;
    ValueType value;
    IndexType index;

    for (IndexType ownedIndex : hmemo::hostReadAccess(ownedIndexes)) {
        Acquisition::coordinate3D coordinate = modelCoordinates.index2coordinate(ownedIndex);

        const IndexType &x = coordinate.x;
        const IndexType &y = coordinate.y;
        const IndexType &z = coordinate.z;

        const IndexType &NX = modelCoordinates.getNX();
        const IndexType &NZ = modelCoordinates.getNZ();

        const IndexType &layer = modelCoordinates.getLayer(coordinate);
        const IndexType &dhFactor = modelCoordinates.getDHFactor(layer);
        dhFactorFineGrid = dhFactor;
        denom = ValueType(1) / ValueType(dhFactor * dhFactor);

        if (modelCoordinates.locatedOnInterface(coordinate)) {
            if (modelCoordinates.getTransition(coordinate) == 1)
                dhFactorFineGrid = modelCoordinates.getDHFactor(layer - 1);
            else if (modelCoordinates.getTransition(coordinate) == -1)
                dhFactorFineGrid = modelCoordinates.getDHFactor(layer + 1);
        }

        if (!modelCoordinates.locatedOnInterface(coordinate)) {
            assembly.push(ownedIndex, ownedIndex, ValueType(1));
        } else {
            modx = (x - int(dhFactor / 2)) % dhFactor;
            modz = z % dhFactor;

            // Point 11
            if (x >= modx) {
                value = (dhFactor - modx) * (dhFactor - modz) * denom;
                index = modelCoordinates.coordinate2index(x - modx, y, z - modz);
                assembly.push(ownedIndex, index, value);
            }

            // Point 12
            if (x + 1 * dhFactorFineGrid < NX) {
                value = modx * (dhFactor - modz) * denom;
                index = modelCoordinates.coordinate2index(x + dhFactor - modx, y, z - modz);
                assembly.push(ownedIndex, index, value);
            }

            //Point 21
            if ((x >= modx) && (z + 2 * dhFactorFineGrid < NZ)) {
                value = (dhFactor - modx) * modz * denom;
                try {
                    index = modelCoordinates.coordinate2index(x - modx, y, z + dhFactor - modz);
                } catch (scai::common::Exception &e) {
                    COMMON_THROWEXCEPTION("Error: coordinate2index called from calcInterpolationStaggeredX in derivatives.cpp" << e.what());
                }

                assembly.push(ownedIndex, index, value);
            }

            //Point 22
            if ((x + 1 * dhFactorFineGrid < NX) && (z + 2 * dhFactorFineGrid < NZ)) {
                value = modx * modz * denom;
                try {
                    index = modelCoordinates.coordinate2index(x + dhFactor - modx, y, z + dhFactor - modz);
                } catch (scai::common::Exception &e) {
                    COMMON_THROWEXCEPTION("Error: coordinate2index called from calcInterpolationStaggeredX in derivatives.cpp" << e.what());
                }
                assembly.push(ownedIndex, index, value);
            }
        }
    }
    // auto colDist = std::make_shared<dmemo::NoDistribution( dist.getGlobalSize() );
    InterpolationStaggeredX = lama::zero<SparseFormat>(dist, dist);

    InterpolationStaggeredX.fillFromAssembly(assembly);
}

//! \brief Calculate interpolation Matrix acoustic (2D/3D) variable grid simulations
/*!
 * Bilinear interpolation is used for the interpolation
 * 
 *   21 ---o---o---22
 *    |            |
 *    o    o   o   o
 *    |            |
 * z  o    o   o   o
 * ^  |            |
 * | 11 ---o---o---21
 *   --> x
\param modelCoordinates Coordinate class, which eg. maps 3D coordinates to 1D model indices
 \param dist Distribution
 */
template <typename ValueType>
void KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::calcInterpolationStaggeredZ(Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::DistributionPtr dist)
{
    SCAI_REGION("Derivatives.calcInterpolationStaggeredZ")

    hmemo::HArray<IndexType> ownedIndexes; // all (global) points owned by this process

    dist->getOwnedIndexes(ownedIndexes);

    lama::MatrixAssembly<ValueType> assembly;
    //assembly.reserve(ownedIndexes.size() * 2);

    ValueType denom = 0;
    IndexType dhFactorFineGrid = 0;
    IndexType modx = 0;
    IndexType modz = 0;
    ValueType value;
    IndexType index;

    for (IndexType ownedIndex : hmemo::hostReadAccess(ownedIndexes)) {
        Acquisition::coordinate3D coordinate = modelCoordinates.index2coordinate(ownedIndex);

        const IndexType &x = coordinate.x;
        const IndexType &y = coordinate.y;
        const IndexType &z = coordinate.z;

        const IndexType &NX = modelCoordinates.getNX();
        const IndexType &NZ = modelCoordinates.getNZ();

        const IndexType &layer = modelCoordinates.getLayer(coordinate);
        const IndexType &dhFactor = modelCoordinates.getDHFactor(layer);
        dhFactorFineGrid = dhFactor;
        denom = ValueType(1) / ValueType(dhFactor * dhFactor);

        if (modelCoordinates.locatedOnInterface(coordinate)) {
            if (modelCoordinates.getTransition(coordinate) == 1)
                dhFactorFineGrid = modelCoordinates.getDHFactor(layer - 1);
            else if (modelCoordinates.getTransition(coordinate) == -1)
                dhFactorFineGrid = modelCoordinates.getDHFactor(layer + 1);
        }

        if (!modelCoordinates.locatedOnInterface(coordinate)) {
            assembly.push(ownedIndex, ownedIndex, ValueType(1));
        } else {
            modx = (x) % dhFactor;
            modz = (z - int(dhFactor / 2)) % dhFactor;

            // Point 11
            if (z >= modz) {
                value = (dhFactor - modx) * (dhFactor - modz) * denom;
                index = modelCoordinates.coordinate2index(x - modx, y, z - modz);
                assembly.push(ownedIndex, index, value);
            }

            // Point 12
            if (x + 2 * dhFactorFineGrid < NX) {
                value = modx * (dhFactor - modz) * denom;
                index = modelCoordinates.coordinate2index(x + dhFactor - modx, y, z - modz);
                assembly.push(ownedIndex, index, value);
            }

            //Point 21
            if ((x >= modx) && (z + 1 * dhFactorFineGrid < NZ)) {
                value = (dhFactor - modx) * modz * denom;
                try {
                    index = modelCoordinates.coordinate2index(x - modx, y, z + dhFactor - modz);
                } catch (scai::common::Exception &e) {
                    COMMON_THROWEXCEPTION("Error: coordinate2index called from calcInterpolationStaggeredX in derivatives.cpp" << e.what());
                }

                assembly.push(ownedIndex, index, value);
            }

            //Point 22
            if ((x + 2 * dhFactorFineGrid < NX) && (z + 1 * dhFactorFineGrid < NZ)) {
                value = modx * modz * denom;
                try {
                    index = modelCoordinates.coordinate2index(x + dhFactor - modx, y, z + dhFactor - modz);
                } catch (scai::common::Exception &e) {
                    COMMON_THROWEXCEPTION("Error: coordinate2index called from calcInterpolationStaggeredX in derivatives.cpp" << e.what());
                }
                assembly.push(ownedIndex, index, value);
            }
        }
    }
    // auto colDist = std::make_shared<dmemo::NoDistribution( dist.getGlobalSize() );
    InterpolationStaggeredZ = lama::zero<SparseFormat>(dist, dist);

    InterpolationStaggeredZ.fillFromAssembly(assembly);
}

//! \brief Calculate interpolation Matrix acoustic (2D/3D) variable grid simulations
/*!
 * Bilinear interpolation is used for the interpolation
 * 
 *   21 ---o---o---22
 *    |            |
 *    o    o   o   o
 *    |            |
 * z  o    o   o   o
 * ^  |            |
 * | 11 ---o---o---12
 *   --> x
\param modelCoordinates Coordinate class, which eg. maps 3D coordinates to 1D model indices
 \param dist Distribution
 */
template <typename ValueType>
void KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::calcInterpolationStaggeredXZ(Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::DistributionPtr dist)
{
    SCAI_REGION("Derivatives.calcInterpolationStaggeredXZ")

    hmemo::HArray<IndexType> ownedIndexes; // all (global) points owned by this process

    dist->getOwnedIndexes(ownedIndexes);

    lama::MatrixAssembly<ValueType> assembly;
    //assembly.reserve(ownedIndexes.size() * 2);

    ValueType denom = 0;
    IndexType dhFactorFineGrid = 0;
    IndexType modx = 0;
    IndexType modz = 0;
    ValueType value;
    IndexType index;

    for (IndexType ownedIndex : hmemo::hostReadAccess(ownedIndexes)) {
        Acquisition::coordinate3D coordinate = modelCoordinates.index2coordinate(ownedIndex);

        const IndexType &x = coordinate.x;
        const IndexType &y = coordinate.y;
        const IndexType &z = coordinate.z;

        const IndexType &NX = modelCoordinates.getNX();
        const IndexType &NZ = modelCoordinates.getNZ();

        const IndexType &layer = modelCoordinates.getLayer(coordinate);
        const IndexType &dhFactor = modelCoordinates.getDHFactor(layer);
        dhFactorFineGrid = dhFactor;
        denom = ValueType(1) / ValueType(dhFactor * dhFactor);

        if (modelCoordinates.locatedOnInterface(coordinate)) {
            if (modelCoordinates.getTransition(coordinate) == 1)
                dhFactorFineGrid = modelCoordinates.getDHFactor(layer - 1);
            else if (modelCoordinates.getTransition(coordinate) == -1)
                dhFactorFineGrid = modelCoordinates.getDHFactor(layer + 1);
        }

        if (!modelCoordinates.locatedOnInterface(coordinate)) {
            assembly.push(ownedIndex, ownedIndex, ValueType(1));
        } else {
            modx = (x - int(dhFactor / 2)) % dhFactor;
            modz = (z - int(dhFactor / 2)) % dhFactor;

            // Point 11
            if ((x >= modx) && (z >= modz)) {
                value = (dhFactor - modx) * (dhFactor - modz) * denom;
                index = modelCoordinates.coordinate2index(x - modx, y, z - modz);
                assembly.push(ownedIndex, index, value);
            }

            // Point 12
            if (x + 1 * dhFactorFineGrid < NX) {
                value = modx * (dhFactor - modz) * denom;
                index = modelCoordinates.coordinate2index(x + dhFactor - modx, y, z - modz);
                assembly.push(ownedIndex, index, value);
            }

            //Point 21
            if ((x >= modx) && (z + 1 * dhFactorFineGrid < NZ)) {
                value = (dhFactor - modx) * modz * denom;
                try {
                    index = modelCoordinates.coordinate2index(x - modx, y, z + dhFactor - modz);
                } catch (scai::common::Exception &e) {
                    COMMON_THROWEXCEPTION("Error: coordinate2index called from calcInterpolationStaggeredX in derivatives.cpp" << e.what());
                }

                assembly.push(ownedIndex, index, value);
            }

            //Point 22
            if ((x + 1 * dhFactorFineGrid < NX) && (z + 1 * dhFactorFineGrid < NZ)) {
                value = modx * modz * denom;
                try {
                    index = modelCoordinates.coordinate2index(x + dhFactor - modx, y, z + dhFactor - modz);
                } catch (scai::common::Exception &e) {
                    COMMON_THROWEXCEPTION("Error: coordinate2index called from calcInterpolationStaggeredX in derivatives.cpp" << e.what());
                }
                assembly.push(ownedIndex, index, value);
            }
        }
    }
    // auto colDist = std::make_shared<dmemo::NoDistribution( dist.getGlobalSize() );
    InterpolationStaggeredXZ = lama::zero<SparseFormat>(dist, dist);

    InterpolationStaggeredXZ.fillFromAssembly(assembly);
}

//! \brief set variable FDorder
/*!
\param FDorderFilename Filename of the FDorder File. This File contains the FDorder for every layer in the variable grid.
 */
template <typename ValueType>
void KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::setFDOrder(std::string const &gridConfigFileName)
{

    unsigned int column = 2;
    Common::readColumnFromFile(gridConfigFileName, spatialFDorderVec, column);

    for (auto i : spatialFDorderVec) {
        if (stencilFDmap.find(i) == stencilFDmap.end())
            COMMON_THROWEXCEPTION("spatialFDorder = " << i << " Unsupported spatialFDorder value.");
    }
}

//! \brief set variable FDorder
/*!
\param FDorder std vector of FDorders per layer
 */
template <typename ValueType>
void KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::setFDOrder(std::vector<scai::IndexType> &FDorder)
{
    for (auto i : FDorder) {
        if (stencilFDmap.find(i) == stencilFDmap.end())
            COMMON_THROWEXCEPTION("spatialFDorder = " << i << " Unsupported spatialFDorder value.");
    }
    spatialFDorderVec = FDorder;
}

//! \brief set constant FDorder
/*!
\param FDorder spatial FD order
 */
template <typename ValueType>
void KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::setFDOrder(IndexType FDorder)
{
    spatialFDorderVec.push_back(FDorder);
    if (stencilFDmap.find(FDorder) == stencilFDmap.end())
        COMMON_THROWEXCEPTION("spatialFDorder = " << FDorder << " Unsupported spatialFDorder value.");
}
//! \brief calculate and return memory usage the of a stencil matrix
/*!
 */
template <typename ValueType>
ValueType KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::getMemoryStencilMatrix(scai::dmemo::DistributionPtr dist)
{
    /* size of a stencil Matrix= size of the halo area which is a CSR Storage of size numGridpoints * IndexType 
     Other contributions (stencil values neighborship relations etc. are neglectable */

    IndexType numGridpoints = dist->getGlobalSize();
    return (numGridpoints * sizeof(IndexType));
}

//! \brief calculate and return memory usage the of a sparse matrix with constant FD order
/*!
 * \param dist distribution
 */
template <typename ValueType>
ValueType KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::getMemorySparseMatrix(scai::dmemo::DistributionPtr dist)
{
    /* size of a sparse Matrix = Number of non Zero Values + Indexes (numRows*spatialFDorder) + num rows (IA -> Indexes per row) +
     numRows * IndexType (halo) */

    IndexType numGridpoints = dist->getGlobalSize();
    return (numGridpoints * spatialFDorderVec.at(0) * (sizeof(ValueType) + sizeof(IndexType)) + 2 * numGridpoints * sizeof(IndexType));
}

//! \brief calculate and return memory usage the of a sparse matrix with constant FD order
/*!
 * \param dist distribution
 */
template <typename ValueType>
ValueType KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::getMemoryInterpolationMatrix(scai::dmemo::DistributionPtr dist)
{
    /* size of a sparse Matrix = Number of non Zero Values + Indexes (numRows*spatialFDorder) + num rows (IA -> Indexes per row) +
     numRows * IndexType (halo). 
     Here only the size of the identity matrix is calculated the actual interpolation planes are neglectable*/

    IndexType numGridpoints = dist->getGlobalSize();
    return (numGridpoints * (sizeof(ValueType) + sizeof(IndexType)) + 2 * numGridpoints * sizeof(IndexType));
}

//! \brief calculate and return memory usage the of a sparse matrix with constant FD order
/*!
 * \param dist distribution
 */
template <typename ValueType>
ValueType KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::getMemorySparseMatrix(scai::dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> const &modelCoordinates)
{
    /* size of a sparse Matrix = Number of non Zero Values + Indexes (numRows*spatialFDorder) + num rows (IA -> Indexes per row) +
     numRows * IndexType (halo) */
    ValueType size = 0;
    auto spatialFDOrder = spatialFDorderVec.at(0);
    for (IndexType layer = 0; layer < modelCoordinates.getNumLayers(); layer++) {
        if (useVarFDorder) {
            spatialFDOrder = spatialFDorderVec.at(layer);
        }
        size += modelCoordinates.getNGridpoints(layer) * spatialFDOrder * (sizeof(ValueType) + sizeof(IndexType)) + 2 * modelCoordinates.getNGridpoints(layer) * sizeof(IndexType);
    }

    return (size);
}

//! \brief Getter method for the spatial FD-order
template <typename ValueType>
IndexType KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::getSpatialFDorder() const
{
    return (spatialFDorderVec.at(0));
}

//! \brief Getter method for derivative matrix DybFreeSurface
template <typename ValueType>
scai::lama::Matrix<ValueType> const &KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::getDybFreeSurface() const
{
    if (!useHybridFreeSurface)
        return DybFreeSurfaceSparse;
    else
        return DybFreeSurfaceHybrid;
}

//! \brief Getter method for derivative matrix DybFreeSurface
template <typename ValueType>
scai::lama::Matrix<ValueType> &KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::getDybFreeSurface()
{
    if (!useHybridFreeSurface)
        return DybFreeSurfaceSparse;
    else
        return DybFreeSurfaceHybrid;
}
//! \brief Getter method for derivative matrix DybFreeSurface
template <typename ValueType>
scai::lama::Matrix<ValueType> const &KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::getDybStaggeredXFreeSurface() const
{
    if ((isElastic) && (useVarGrid)) {
        return (DybStaggeredXFreeSurface);
    } else if (!useHybridFreeSurface) {
        return DybFreeSurfaceSparse;
    } else {
        return DybFreeSurfaceHybrid;
    }
}

//! \brief Getter method for derivative matrix DybFreeSurface
template <typename ValueType>
scai::lama::Matrix<ValueType> const &KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::getDybStaggeredZFreeSurface() const
{
    if ((isElastic) && (useVarGrid)) {
        return (DybStaggeredZFreeSurface);
    } else if (!useHybridFreeSurface) {
        return DybFreeSurfaceSparse;
    } else {
        return DybFreeSurfaceHybrid;
    }
}

//! \brief Getter method for derivative matrix DyfFreeSurface
template <typename ValueType>
scai::lama::Matrix<ValueType> const &KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::getDyfFreeSurface() const
{
    if (!useHybridFreeSurface)
        return DyfFreeSurfaceSparse;
    else
        return DyfFreeSurfaceHybrid;
}

//! \brief Getter method for derivative matrix DyfFreeSurface
template <typename ValueType>
scai::lama::Matrix<ValueType> &KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::getDyfFreeSurface()
{
    if (!useHybridFreeSurface)
        return DyfFreeSurfaceSparse;
    else
        return DyfFreeSurfaceHybrid;
}

//! \brief Getter method for derivative matrix Dxf
template <typename ValueType>
scai::lama::Matrix<ValueType> const &KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::getDxf() const
{
    if (useStencilMatrix)
        return (Dxf);
    else
        return (DxfSparse);
}

//! \brief Getter method for derivative matrix Dyf
template <typename ValueType>
scai::lama::Matrix<ValueType> const &KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::getDyf() const
{
    if (useStencilMatrix)
        return (Dyf);
    else
        return (DyfSparse);
}

//! \brief Getter method for derivative matrix Dyf
template <typename ValueType>
scai::lama::Matrix<ValueType> &KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::getDyf()
{
    if (useStencilMatrix)
        return (Dyf);
    else
        return (DyfSparse);
}

//! \brief Getter method for derivative matrix DyfStaggeredX
template <typename ValueType>
scai::lama::Matrix<ValueType> const &KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::getDyfStaggeredX() const
{
    if ((isElastic) && (useVarGrid))
        return (DyfStaggeredXSparse);
    else if (useStencilMatrix)
        return (Dyf);
    else
        return (DyfSparse);
}

//! \brief Getter method for derivative matrix DybStaggeredX
template <typename ValueType>
scai::lama::Matrix<ValueType> const &KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::getDybStaggeredX() const
{
    if ((isElastic) && (useVarGrid))
        return (DybStaggeredXSparse);
    else if (useStencilMatrix)
        return (Dyb);
    else
        return (DybSparse);
}
//! \brief Getter method for derivative matrix DyfStaggeredZ
template <typename ValueType>
scai::lama::Matrix<ValueType> const &KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::getDyfStaggeredZ() const
{
    if ((isElastic) && (useVarGrid))
        return (DyfStaggeredZSparse);
    else if (useStencilMatrix)
        return (Dyf);
    else
        return (DyfSparse);
}

//! \brief Getter method for derivative matrix DybStaggeredX
template <typename ValueType>
scai::lama::Matrix<ValueType> const &KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::getDybStaggeredZ() const
{
    if ((isElastic) && (useVarGrid))
        return (DybStaggeredZSparse);
    else if (useStencilMatrix)
        return (Dyb);
    else
        return (DybSparse);
}

//! \brief Getter method for derivative matrix Dzf
template <typename ValueType>
scai::lama::Matrix<ValueType> const &KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::getDzf() const
{
    if (useStencilMatrix)
        return (Dzf);
    else
        return (DzfSparse);
}

//! \brief Getter method for derivative matrix Dxb
template <typename ValueType>
scai::lama::Matrix<ValueType> const &KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::getDxb() const
{
    if (useStencilMatrix)
        return (Dxb);
    else
        return (DxbSparse);
}

//! \brief Getter method for derivative matrix Dyb
template <typename ValueType>
scai::lama::Matrix<ValueType> const &KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::getDyb() const
{
    if (useStencilMatrix)
        return (Dyb);
    else
        return (DybSparse);
}

//! \brief Getter method for derivative matrix Dzb
template <typename ValueType>
scai::lama::Matrix<ValueType> const &KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::getDzb() const
{
    if (useStencilMatrix)
        return (Dzb);
    else
        return (DzbSparse);
}

//! \brief Getter method for derivative interpolation matrix of P
template <typename ValueType>
scai::lama::Matrix<ValueType> const *KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::getInterFull() const
{
    if (InterpolationFull.getNumRows() > 0) {
        return &InterpolationFull;
    } else {
        return NULL;
    }
}

//! \brief Getter method for derivative interpolation matrix of points staggered in x direction
template <typename ValueType>
scai::lama::Matrix<ValueType> const *KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::getInterStaggeredX() const
{
    if (InterpolationStaggeredX.getNumRows() > 0) {
        return &InterpolationStaggeredX;
    } else {
        return NULL;
    }
}

//! \brief Getter method for derivative interpolation matrix of points staggered in x direction
template <typename ValueType>
scai::lama::Matrix<ValueType> const *KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::getInterStaggeredZ() const
{
    if (InterpolationStaggeredZ.getNumRows() > 0) {
        return &InterpolationStaggeredZ;
    } else {
        return NULL;
    }
}

//! \brief Getter method for derivative interpolation matrix of points staggered in x direction
template <typename ValueType>
scai::lama::Matrix<ValueType> const *KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::getInterStaggeredXZ() const
{
    if (InterpolationStaggeredXZ.getNumRows() > 0) {
        return &InterpolationStaggeredXZ;
    } else {
        return NULL;
    }
}

//! \brief Set FD coefficients for each order
/*!
 *
 */
template <typename ValueType>
void KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::setFDCoef()
{
    const ValueType FD2[] = {-1.0, 1.0};
    stencilFDmap[2] = common::Stencil1D<ValueType>(2, FD2);

    const ValueType FD4[] = {1.0 / 24.0, -9.0 / 8.0, 9.0 / 8.0, -1.0 / 24.0};
    stencilFDmap[4] = common::Stencil1D<ValueType>(4, FD4);

    const ValueType FD6[] = {-3.0 / 640.0, 25.0 / 384.0, -75.0 / 64.0,
                             75.0 / 64.0, -25.0 / 384.0, 3.0 / 640.0};
    stencilFDmap[6] = common::Stencil1D<ValueType>(6, FD6);

    const ValueType FD8[] = {5.0 / 7168.0, -49.0 / 5120.0, 245.0 / 3072.0, -1225.0 / 1024.0,
                             1225.0 / 1024.0, -245.0 / 3072.0, 49.0 / 5120.0, -5.0 / 7168.0};
    stencilFDmap[8] = common::Stencil1D<ValueType>(8, FD8);

    const ValueType FD10[] = {-8756999275442633.0 / 73786976294838206464.0,
                              8142668969129685.0 / 4611686018427387904.0,
                              -567.0 / 40960.0,
                              735.0 / 8192.0,
                              -19845.0 / 16384.0,
                              19845.0 / 16384.0,
                              -735.0 / 8192.0,
                              567.0 / 40960.0,
                              -8142668969129685.0 / 4611686018427387904.0,
                              8756999275442633.0 / 73786976294838206464.0};
    stencilFDmap[10] = common::Stencil1D<ValueType>(10, FD10);

    const ValueType FD12[] = {6448335830095439.0 / 295147905179352825856.0,
                              -1655620175512543.0 / 4611686018427387904.0,
                              6842103786556949.0 / 2305843009213693952.0,
                              -628618285389933.0 / 36028797018963968.0,
                              436540475965291.0 / 4503599627370496.0,
                              -2750204998582123.0 / 2251799813685248.0,
                              2750204998582123.0 / 2251799813685248.0,
                              -436540475965291.0 / 4503599627370496.0,
                              628618285389933.0 / 36028797018963968.0,
                              -6842103786556949.0 / 2305843009213693952.0,
                              1655620175512543.0 / 4611686018427387904.0,
                              -6448335830095439.0 / 295147905179352825856.0};
    stencilFDmap[12] = common::Stencil1D<ValueType>(12, FD12);
}

template class KITGPI::ForwardSolver::Derivatives::Derivatives<float>;
template class KITGPI::ForwardSolver::Derivatives::Derivatives<double>;
