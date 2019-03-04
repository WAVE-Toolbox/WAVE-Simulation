#include "Coordinates.hpp"
#include <algorithm>
#include <vector>
#include <cmath>
#include <iterator> 
using namespace scai;

/* ------- */
/* Mapping */
/* ------- */

/*! \brief constructor for variable grid
 *
 \param NX Number of grid points in X
 \param NY Number of grid points in Y
 \param NZ Number of grid points in Z
 *
 */
template <typename ValueType>
KITGPI::Acquisition::Coordinates<ValueType>::Coordinates(IndexType nx, IndexType ny, IndexType nz, ValueType dh, std::vector<IndexType> &dhFactors, std::vector<IndexType> &interfaces) : NX(nx), NY(ny), NZ(nz), DH(dh), dhFactor(dhFactors), interface(interfaces)
{
    
    SCAI_ASSERT_ERROR(NX > 0, "NX<=0");
    SCAI_ASSERT_ERROR(NY > 0, "NY<=0");
    SCAI_ASSERT_ERROR(NZ > 0, "NZ<=0");
    
    numLayers = dhFactor.size();

    varNX.resize(numLayers);
    varNY.resize(numLayers);
    varNZ.resize(numLayers);
    varDH.resize(numLayers);

    layerStart.resize(numLayers);
    layerEnd.resize(numLayers);
    transision.resize(numLayers);

    nGridpointsPerLayer.resize(numLayers);
    nGridpoints = 0;

    interface.insert(interface.begin(), -1);
    interface.push_back(ny - 1);

    /* loop over all layers in the variable grid 
     DH will be multiplied with the enlargement factor dhFactor.  dhFactor must be a multiple of 3.
     The begin and end of each layer will be estimated from the layer interfaces. 
     Begin and end doesn't match the interface locations because internally the ghost points are part of the fine grid instead of the coarse grid. This simplyfies the mapping.
     Also the number of gridpoints in each dimension and in total will be estimated. 
     */
    for (IndexType layer = 0; layer < numLayers - 1; layer++) {
        if (dhFactor[layer] < dhFactor[layer + 1]) {
            transision[layer] = 1;
            layerEnd[layer] = interface[layer + 1];
            layerStart[layer + 1] = interface[layer + 1] + dhFactor[layer + 1];
        } else {
            transision[layer] = 0;
            layerEnd[layer] = interface[layer + 1] - dhFactor[layer];
            layerStart[layer + 1] = interface[layer + 1];
        }
    }
    
    // set the last elements
    transision[numLayers - 1] = 0;
    layerEnd[numLayers - 1] = interface[numLayers];

    for (IndexType layer = 0; layer < numLayers; layer++) {

        varNY[layer] = (layerEnd[layer] - layerStart[layer]) / dhFactor[layer] + 1;
        varNX[layer] = floor(nx / dhFactor[layer]);
        varNZ[layer] = floor(nz / dhFactor[layer]);
        
       
        
        if (dhFactor[layer] > 1) {
            varNX[layer] = varNX[layer] + 1;
            varNZ[layer] = varNZ[layer] + 1;
        }

        
         
        nGridpointsPerLayer[layer] = varNX[layer] * varNY[layer] * varNZ[layer];
        nGridpoints += nGridpointsPerLayer[layer];
        
  //       std::cout << "nGridpointsPerLayer [" << layer << "] = "<< nGridpointsPerLayer[layer] << std::endl;
        varDH[layer]=dh*dhFactor[layer];
    }
    std::cout << "nGridpoints " << nGridpoints << std::endl;
}


/*! \brief constructor for regular grid
 *
 \param NX Number of grid points in X
 \param NY Number of grid points in Y
 \param NZ Number of grid points in Z
 *
 */
template <typename ValueType>
KITGPI::Acquisition::Coordinates<ValueType>::Coordinates(scai::IndexType nx, scai::IndexType ny, scai::IndexType nz, ValueType dh) : NX(nx), NY(ny), NZ(nz), DH(dh)
{
    numLayers = 1;

    varNX.push_back(nx);
    varNY.push_back(ny);
    varNZ.push_back(nz);
    varDH.push_back(dh);
    
    layerStart.push_back(0);
    layerEnd.push_back(ny - 1);
    // transision.resize(numLayers);

    nGridpoints = nx * ny * nz;

    nGridpointsPerLayer.push_back(nGridpoints);

    interface.insert(interface.begin(), 0);
    interface.push_back(ny - 1);
    dhFactor.push_back(1);
    transision.push_back(0);

    SCAI_ASSERT_ERROR(NX > 0, "NX<=0");
    SCAI_ASSERT_ERROR(NY > 0, "NY<=0");
    SCAI_ASSERT_ERROR(NZ > 0, "NZ<=0");
}

/*! \brief getter function for DH
 *
 *
 */
template <typename ValueType>
ValueType KITGPI::Acquisition::Coordinates<ValueType>::getDH() const
{
    return (DH);
}

/*! \brief getter function for NX
 *
 *
 */
template <typename ValueType>
scai::IndexType KITGPI::Acquisition::Coordinates<ValueType>::getNX() const
{
    return (NX);
}

/*! \brief getter function for NY
 *
 *
 */
template <typename ValueType>
scai::IndexType KITGPI::Acquisition::Coordinates<ValueType>::getNY() const
{
    return (NY);
}

/*! \brief getter function for NZ
 *
 *
 */
template <typename ValueType>
scai::IndexType KITGPI::Acquisition::Coordinates<ValueType>::getNZ() const
{
    return (NZ);
}

/*! \brief getter function for numberOfGridpoints
 *
 *
 */
template <typename ValueType>
scai::IndexType KITGPI::Acquisition::Coordinates<ValueType>::getNGridpoints() const
{
    return (nGridpoints);
}

/*! \brief getter function for DH
 *
 *
 */
template <typename ValueType>
ValueType KITGPI::Acquisition::Coordinates<ValueType>::getDH(coordinate3D coordinate) const
{
    IndexType layer = 0;

    for (layer = 0; layer < numLayers; ++layer) {
        if ((coordinate.y < interface[layer + 1]) && (coordinate.y > interface[layer])) {
            break;
        } else if (coordinate.y == interface[layer + 1]) {
            layer += transision[layer];
            break;
        } 
    }
//     if (varDH[layer]== 0)
//         std::cout << "error at y = " << coordinate.y << " layer = " << layer << " interface = " << interface[0] << " " << interface[1] << " transition = " << transision[0] << " " << transision[1] << std::endl;

    return (varDH[layer]);
}

/*! \brief getter function for dhFactor
 *
 *
 */
template <typename ValueType>
IndexType KITGPI::Acquisition::Coordinates<ValueType>::getDHFactor(coordinate3D coordinate) const
{
    IndexType layer = 0;

    for (layer = 0; layer < numLayers; layer++) {
        if ((coordinate.y < interface[layer + 1]) && (coordinate.y > interface[layer])) {
            break;
        } else if (coordinate.y == interface[layer + 1]) {
            layer += transision.at(layer);
            break;
        }
    }

    return (dhFactor[layer]);
}
/*! \brief getter function for NX
 *
 *
 */
template <typename ValueType>
scai::IndexType KITGPI::Acquisition::Coordinates<ValueType>::getNX(coordinate3D coordinate) const
{
    IndexType layer;

    for (layer = 0; layer < numLayers; layer++) {
        if ((coordinate.y <= layerEnd[layer]) && (coordinate.y >= layerStart[layer])) {
            break;
        }
    }

    return (varNX[layer]);
}

/*! \brief getter function for NY
 *
 *
 */
template <typename ValueType>
scai::IndexType KITGPI::Acquisition::Coordinates<ValueType>::getNY(coordinate3D coordinate) const
{
    IndexType layer;

    for (layer = 0; layer < numLayers; layer++) {
        if ((coordinate.y <= layerEnd[layer]) && (coordinate.y >= layerStart[layer])) {
            break;
        }
    }

    return (varNY[layer]);
}

/*! \brief getter function for NZ
 *
 *
 */
template <typename ValueType>
scai::IndexType KITGPI::Acquisition::Coordinates<ValueType>::getNZ(coordinate3D coordinate) const
{
    IndexType layer;

    for (layer = 0; layer < numLayers; layer++) {
        if ((coordinate.y <= layerEnd[layer]) && (coordinate.y >= layerStart[layer])) {
            break;
        }
    }

    return (varNZ[layer]);
}

/*! \brief getter function for coordinates
 *
 *
 */
template <typename ValueType>
std::vector<scai::lama::DenseVector<ValueType>> KITGPI::Acquisition::Coordinates<ValueType>::getCoordinates(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx) const
{

    hmemo::HArray<IndexType> ownedIndexes; // all (global) points owned by this process
    dist->getOwnedIndexes(ownedIndexes);

    std::vector<lama::VectorAssembly<ValueType>> assembly(3);
    for (int i = 0; i < 3; i++) {
        assembly[i].reserve(ownedIndexes.size());
    }

    for (IndexType ownedIndex : hmemo::hostReadAccess(ownedIndexes)) {
        coordinate3D coordinate = index2coordinate(ownedIndex);

        assembly[0].push(ownedIndex, ValueType(coordinate.x));
        assembly[1].push(ownedIndex, ValueType(coordinate.y));
        assembly[2].push(ownedIndex, ValueType(coordinate.z));
    }

    std::vector<scai::lama::DenseVector<ValueType>> coordinates;
    coordinates.resize(3);

    for (int i = 0; i < 3; i++) {
        coordinates[i].allocate(dist);
        coordinates[i].setContextPtr(ctx);
        coordinates[i].fillFromAssembly(assembly[i]);
    }

    return (coordinates);
}

/*! \brief Returns bool if given coordinate is located on the surface
 *
 * This method determines if a given coordinate is located on the surface of the modelling domain.
 \param index Model vector Index
 *
 */
template <typename ValueType>
bool KITGPI::Acquisition::Coordinates<ValueType>::locatedOnSurface(IndexType index) const
{
    coordinate3D result;
    result = mapIndex2coordinate(index);
    if (result.y == 0) {
        return (true);
    } else {
        return (false);
    }
}

/*! \brief General mapping from 1-D index to 3-D coordinate
 *
 \param index Model vector Index
 *
 */
template <typename ValueType>
KITGPI::Acquisition::coordinate3D KITGPI::Acquisition::Coordinates<ValueType>::mapIndex2coordinate(IndexType index) const
{
    IndexType layer;
    // find out which in which layer the index is and reduce index to index inside this layer
    for (layer = 0; layer < numLayers; layer++) {
        index -= nGridpointsPerLayer[layer];
        if (index < 0) {
            index += nGridpointsPerLayer[layer];
            break;
        }
    }

    coordinate3D result;

    result.y = IndexType(index / (varNX[layer] * varNZ[layer]));
    index -= result.y * (varNX[layer] * varNZ[layer]);

    result.z = IndexType(index / (varNX[layer]));
    index -= result.z * (varNX[layer]);

    result.x = index;

    // coordinates in reference to the fine grid
    result.x *= dhFactor[layer];
    result.y *= dhFactor[layer];
    result.z *= dhFactor[layer];

    // move the subgrid coordinates to global coordinates
    result.y += layerStart[layer];

    return (result);
}

/*! \brief General mapping from 1-D index to 3-D coordinate
 *
 * Maps a 1-D index into 3-D coordinates.
 * The 3-D grid starts at 0 and runs to (NX-1), (NY-1) or (NZ-1).
 *
 \param index 1-D coordinate
 *
 */
template <typename ValueType>
KITGPI::Acquisition::coordinate3D KITGPI::Acquisition::Coordinates<ValueType>::index2coordinate(IndexType index) const
{
    return (mapIndex2coordinate(index));
}

/*! \brief General mapping from 3-D coordinates to 1-D indeces
 *
 \param X 3-D coordinate in X (Horizontal 1)
 \param Y 3-D coordinate in Y (Depth)
 \param Z 3-D coordinate in Z (Horizontal 2)
 */
template <typename ValueType>
IndexType KITGPI::Acquisition::Coordinates<ValueType>::map3Dcoordinate2index(IndexType X, IndexType Y, IndexType Z) const
{

//     SCAI_ASSERT(Z < getNZ(), "Could not map from coordinate to index!");
//     SCAI_ASSERT(X < getNX(), "Could not map from coordinate to index!");
//     SCAI_ASSERT(Y < getNY(), "Could not map from coordinate to index! " << getNY());
    SCAI_ASSERT(Z >= 0, "Could not map from coordinate to index!");
    SCAI_ASSERT(Y >= 0, "Could not map from coordinate to index!");
    SCAI_ASSERT(X >= 0, "Could not map from coordinate to index!");

    IndexType layer;

    for (layer = 0; layer < numLayers; layer++) {
        if ((Y <= layerEnd[layer]) && (Y >= layerStart[layer])) {
            Y -= layerStart[layer];
            break;
        }
    }

    IndexType index = (X / dhFactor[layer]) + (Z / dhFactor[layer]) * varNX[layer] + (Y / dhFactor[layer]) * varNX[layer] * varNZ[layer];

    for (IndexType l = 1; l <= layer; l++) {
        index += nGridpointsPerLayer[l - 1];
    }

    return (index);
}

/* ---------- */
/* Interfaces */
/* ---------- */

/*! \brief Convert 3-D coordinates to 1-D indeces
 *
 * This method returns the 1-D coordinate of 3-D coordinates.
 \param X 3-D coordinate in X (Horizontal 1)
 \param Y 3-D coordinate in Y (Depth)
 \param Z 3-D coordinate in Z (Horizontal 2)
 \return 1-D index
 */
template <typename ValueType>
IndexType KITGPI::Acquisition::Coordinates<ValueType>::coordinate2index(IndexType X, IndexType Y, IndexType Z) const
{
    return (map3Dcoordinate2index(X, Y, Z));
}

/*! \brief Convert 3-D coordinates to 1-D coordinates
 *
 * This method returns the 1-D coordinate of 3-D coordinates.
 \param coordinate as a coordinate3D struct
 \return 1-D index
 */
template <typename ValueType>
IndexType KITGPI::Acquisition::Coordinates<ValueType>::coordinate2index(coordinate3D coordinate) const
{
    return (map3Dcoordinate2index(coordinate.x, coordinate.y, coordinate.z));
}

/*! \brief Calculation of distance to boundaries of the modelling domain
 *
 * This method calculates the distance of a given coordinate to the boundaries of the modelling domain.
 \param X Coordinate X
 \param Y Coordinate Y
 \param Z Coordinate Z
 */
template <typename ValueType>
KITGPI::Acquisition::coordinate3D KITGPI::Acquisition::Coordinates<ValueType>::estimateDistanceToEdges3D(IndexType X, IndexType Y, IndexType Z) const
{
//     SCAI_ASSERT(Z < NZ[0], "No valid argument!");
//     SCAI_ASSERT(X < NX[0], "No valid argument!");
//     SCAI_ASSERT(Y < NY[0], "No valid argument!");
    SCAI_ASSERT(Z >= 0, "No valid argument!");
    SCAI_ASSERT(Y >= 0, "No valid argument!");
    SCAI_ASSERT(X >= 0, "No valid argument!");

    coordinate3D distance;

//     distance.x = !((NX[0] - 1 - X) < (X)) ? (X) : (NX[0] - 1 - X);
//     distance.y = !((NY[0] - 1 - Y) < (Y)) ? (Y) : (NY[0] - 1 - Y);
//     distance.z = !((NZ[0] - 1 - Z) < (Z)) ? (Z) : (NZ[0] - 1 - Z);

    return (distance);
}

/*! \brief Determination of distance to boundaries
 *
 * This method calculates the distance of a given coordinate3D to the boundaries of the modelling domain.
 \param coordinate 3D-coordinate structs
 */
template <typename ValueType>
KITGPI::Acquisition::coordinate3D KITGPI::Acquisition::Coordinates<ValueType>::edgeDistance(coordinate3D coordinate) const
{
    return (estimateDistanceToEdges3D(coordinate.x, coordinate.y, coordinate.z));
}

template class KITGPI::Acquisition::Coordinates<float>;
template class KITGPI::Acquisition::Coordinates<double>;
