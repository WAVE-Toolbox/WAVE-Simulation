#include "Coordinates.hpp"
#include <algorithm>
#include <cmath>
#include <iterator>
#include <vector>
using namespace scai;

/* ------- */
/* Mapping */
/* ------- */

/*! \brief constructor for variable grid
 *
 \param config Configuration class
 *
 */
template <typename ValueType>
KITGPI::Acquisition::Coordinates<ValueType>::Coordinates(Configuration::Configuration const &config) : NX(config.get<IndexType>("NX")), NY(config.get<IndexType>("NY")), NZ(config.get<IndexType>("NZ")), DH(config.get<ValueType>("DH"))
{

    if (config.get<bool>("useVariableGrid")) {
        VariableGrid = true;
        std::vector<IndexType> dhFactor;
        std::vector<IndexType> interface;
        std::ifstream is(config.get<std::string>("interfaceFilename"));
        std::istream_iterator<IndexType> start(is), end;
        interface.assign(start, end);
        std::ifstream is2(config.get<std::string>("dhFactorFilename"));
        std::istream_iterator<IndexType> start2(is2);
        dhFactor.assign(start2, end);
        init(dhFactor, interface);
    } else {
        VariableGrid = false;
        init();
    }
}
/*! \brief constructor for variable grid
 *
 \param nx Number of grid points in X
 \param ny Number of grid points in Y
 \param nz Number of grid points in Z
 \param dh finest grid spacing
 \param dhFactors grid spacing factor per layer (3^n)
 \param interfaces std::vector with 3D indeces (with respect to a regular grid) of the interfaces
 *
 */
template <typename ValueType>
KITGPI::Acquisition::Coordinates<ValueType>::Coordinates(IndexType nx, IndexType ny, IndexType nz, ValueType dh, std::vector<IndexType> &dhFactors, std::vector<IndexType> &interfaces) : NX(nx), NY(ny), NZ(nz), DH(dh), dhFactor(dhFactors), interface(interfaces)
{
    VariableGrid = true;
    init(dhFactors, interface);
}

/*! \brief constructor for variable grid
 *
 \param dhFactors std::vector with grid stretch factors f=3^n
 \param interfaces std::vector with 3D indeces (with respect to a regular grid) of the interfaces
 *
 */
template <typename ValueType>
void KITGPI::Acquisition::Coordinates<ValueType>::init(std::vector<IndexType> &dhFactors, std::vector<IndexType> &interfaces)
{
    dhFactor = dhFactors;
    interface = interfaces;
    SCAI_ASSERT_ERROR(dhFactor.size() > 0, "vector of fifferent grid spacings: dhFactor is emty");
    SCAI_ASSERT_ERROR(interface.size() == dhFactor.size() - 1, "number of interfaces doesn't match to the number of different grid spacings");

    SCAI_ASSERT_ERROR(NX > 0, "NX<=0");
    SCAI_ASSERT_ERROR(NY > 0, "NY<=0");
    SCAI_ASSERT_ERROR(NZ > 0, "NZ<=0");

    numLayers = dhFactor.size();

    varNX.assign(numLayers, 0);
    varNY.assign(numLayers, 0);
    varNZ.assign(numLayers, 0);
    varDH.assign(numLayers, 0.0);

    layerStart.assign(numLayers, 0);
    layerEnd.assign(numLayers, 0);
    transition.assign(numLayers, 0);

    nGridpointsPerLayer.assign(numLayers, 0);
    nGridpoints = 0;

    interface.push_back(NY - 1);
    interface.insert(interface.begin(), -1);

    IndexType dhMax = 0;
    /*Check for incompatible NX,NZ,NY interfaces and dhFactors and correct them if necessary*/
    for (IndexType layer = 0; layer < numLayers; layer++) {
        SCAI_ASSERT_ERROR(floor(std::log(dhFactor[layer]) / std::log(3)) == std::log(dhFactor[layer]) / std::log(3), "incompatible dhFactor, dhFactor must be 3^n")
        if (dhFactor[layer] > dhMax)
            dhMax = dhFactor[layer];
    }

    IndexType NXmax = NX;
    IndexType NZmax = NZ;

    if (dhMax != 1) {
        while (true) {
            if (NXmax == floor(NXmax / dhMax) * dhMax + 1 + floor(dhMax / 2)) {
                break;
            }
            NXmax--;
        }

        while ((true) && (NZ != 1)) {
            if (NZmax == floor(NZmax / dhMax) * dhMax + 1 + floor(dhMax / 2)) {
                break;
            }
            NZmax--;
        }

        if (NXmax != NX)
            std::cout << "NX for the variable grid has been changed from " << NX << " to " << NXmax << std::endl;
        if (NZmax != NZ)
            std::cout << "NZ for the variable grid has been changed from " << NZ << " to " << NZmax << std::endl;

        /* This has to be fixed: Because interface[0]=-1 (should be changed to 0) interface 1 has to be checked seperately */
        IndexType layer = 0;
        while (layer < 1) {
            layer++;
            ValueType test = ValueType(interface[layer] - interface[layer - 1] - 1) / dhFactor[layer - 1];
            if (floor(test) != test) {
                std::cout << "interface[" << layer << "] has been changed from " << interface[layer] << " to " << interface[layer] - 1 << std::endl;
                interface[layer]--;
                layer--;
            }
        }

        layer = 1;
        while (layer < numLayers) {
            layer++;
            ValueType test = ValueType(interface[layer] - interface[layer - 1]) / dhFactor[layer - 1];
            if (floor(test) != test) {
                std::cout << "interface[" << layer << "] has been changed from " << interface[layer] << " to " << interface[layer] - 1 << std::endl;
                interface[layer]--;
                layer--;
            }
        }
    }
    /* loop over all layers in the variable grid 
     DH will be multiplied with the enlargement factor dhFactor.  dhFactor must be a multiple of 3.
     The begin and end of each layer will be estimated from the layer interfaces. 
     Begin and end doesn't match the interface locations because internally the ghost points are part of the fine grid instead of the coarse grid. This simplyfies the mapping.
     Also the number of gridpoints in each dimension and in total will be estimated. 
     */
    for (IndexType layer = 0; layer < numLayers - 1; layer++) {
        if (dhFactor[layer] < dhFactor[layer + 1]) {
            transition[layer] = 1;
            layerEnd[layer] = interface[layer + 1];
            layerStart[layer + 1] = interface[layer + 1] + dhFactor[layer + 1];
        } else {
            transition[layer] = 0;
            layerEnd[layer] = interface[layer + 1] - dhFactor[layer];
            layerStart[layer + 1] = interface[layer + 1];
        }
    }

    // set the last elements
    transition[numLayers - 1] = 0;
    layerEnd[numLayers - 1] = interface[numLayers];

    //calculate NX,NY,NZ and number of gridpoints per layer
    for (IndexType layer = 0; layer < numLayers; layer++) {

        varNY[layer] = (layerEnd[layer] - layerStart[layer]) / dhFactor[layer] + 1;
        varNX[layer] = floor(NXmax / dhFactor[layer]);
        varNZ[layer] = floor(NZmax / dhFactor[layer]);

        if (dhFactor[layer] > 1) {
            varNX[layer] = varNX[layer] + 1;
            varNZ[layer] = varNZ[layer] + 1;
        }

        nGridpointsPerLayer[layer] = varNX[layer] * varNY[layer] * varNZ[layer];
        nGridpoints += nGridpointsPerLayer[layer];

        varDH[layer] = DH * dhFactor[layer];
    }
}

/*! \brief constructor for regular grid
 *
 \param nx Number of grid points in X
 \param ny Number of grid points in Y
 \param nz Number of grid points in Z
 \param dh grid spacing in meter
 *
 */
template <typename ValueType>
KITGPI::Acquisition::Coordinates<ValueType>::Coordinates(scai::IndexType nx, scai::IndexType ny, scai::IndexType nz, ValueType dh) : NX(nx), NY(ny), NZ(nz), DH(dh)
{
    VariableGrid = false;
    init();
}

/*! \brief constructor for regular grid
 *
 *
 */
template <typename ValueType>
void KITGPI::Acquisition::Coordinates<ValueType>::init()
{
    numLayers = 1;

    varNX.push_back(NX);
    varNY.push_back(NY);
    varNZ.push_back(NZ);
    varDH.push_back(DH);

    layerStart.push_back(0);
    layerEnd.push_back(NY - 1);
    // transition.resize(numLayers);

    nGridpoints = NX * NY * NZ;

    nGridpointsPerLayer.push_back(nGridpoints);

    interface.push_back(NY - 1);
    interface.insert(interface.begin(), -1);

    dhFactor.push_back(1);
    transition.push_back(0);

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

/*! \brief getter function for the layer
 *
 *  coarse grids contain the interface at the fine<->coarse and coarse<->fine transition 
 */
template <typename ValueType>
IndexType KITGPI::Acquisition::Coordinates<ValueType>::getLayer(coordinate3D coordinate) const
{
    IndexType layer = 0;

    for (layer = 0; layer < numLayers; layer++) {
        if ((coordinate.y < interface[layer + 1]) && (coordinate.y > interface[layer])) {
            break;
        } else if (coordinate.y == interface[layer + 1]) {
            layer += transition.at(layer);
            break;
        }
    }
    return (layer);
}
/*! \brief getter function for number of layers
 *

 */
template <typename ValueType>
IndexType KITGPI::Acquisition::Coordinates<ValueType>::getNumLayers() const
{
    return (numLayers);
}
/*! \brief getter function for DH
 \param layer layer of the variable grid
 *
 */
template <typename ValueType>
ValueType KITGPI::Acquisition::Coordinates<ValueType>::getDH(IndexType layer) const
{
    return (varDH[layer]);
}

/*! \brief getter function for DH
 \param coordinate 3D coordinate struct
 *
 */
template <typename ValueType>
ValueType KITGPI::Acquisition::Coordinates<ValueType>::getDH(coordinate3D coordinate) const
{
    return (varDH[getLayer(coordinate)]);
}

/*! \brief getter function for dhFactor
 *
 \param coordinate 3D coordinate struct
 */
template <typename ValueType>
IndexType KITGPI::Acquisition::Coordinates<ValueType>::getDHFactor(coordinate3D coordinate) const
{
    return (dhFactor[getLayer(coordinate)]);
}

/*! \brief getter function for dhFactor
 *
  \param layer layer of the variable grid
 */
template <typename ValueType>
IndexType KITGPI::Acquisition::Coordinates<ValueType>::getDHFactor(IndexType layer) const
{
    return (dhFactor[layer]);
}

/*! \brief check if coordinate is on an variable grid interface
 *
  \param coordinate 3D coordinate struct
 */
template <typename ValueType>
bool KITGPI::Acquisition::Coordinates<ValueType>::locatedOnInterface(coordinate3D coordinate) const
{
    bool isOnInterface = false;
    for (int layer = 0; layer < numLayers; layer++) {
        if (coordinate.y == interface[layer]) {
            isOnInterface = true;
        }
    }
    return (isOnInterface);
}

/*! \brief calculate distance in gridpoints of the finest grid to the variable grid interfaces
 *
 \param Y y coordinate
 */
template <typename ValueType>
IndexType KITGPI::Acquisition::Coordinates<ValueType>::distToInterface(IndexType Y) const
{
    IndexType dist = NY;
    for (auto i = interface.begin() + 1; i != interface.end() - 1; ++i) {
        IndexType temp = std::abs(Y - *i);
        if (temp < dist)
            dist = temp;
    }
    return (dist);
}

/*! \brief getter function for NZ
  \param coordinate 3D coordinate struct
 *
 */
template <typename ValueType>
bool KITGPI::Acquisition::Coordinates<ValueType>::getTransition(coordinate3D coordinate) const
{
    bool fineToCoarse = false;
    for (int layer = 0; layer < numLayers; layer++) {
        if (coordinate.y == interface[layer + 1]) {
            fineToCoarse = transition[layer];
            // std::cout << layer << " " << fineToCoarse << "  "<< interface[layer+1] << " " << coordinate.y << std::endl;
        }
    }

    return (fineToCoarse);
}

/*! \brief getter function for coordinates
 *
 \param dist Distribution
 \param ctx lama context 
 */
template <typename ValueType>
std::vector<scai::lama::DenseVector<ValueType>> KITGPI::Acquisition::Coordinates<ValueType>::getCoordinates(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx) const
{
    SCAI_ASSERT_ERROR(dist, "no distribution given");

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

/*! \brief getter function for coordinates
 *
 \param dist Distribution
 \param ctx lama context 
 \param filename string with the filename to write the coordinates
 */
template <typename ValueType>
void KITGPI::Acquisition::Coordinates<ValueType>::writeCoordinates(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, std::string filename) const
{
    auto coords = getCoordinates(dist, ctx);
    coords[0].writeToFile(filename + "X.mtx");
    coords[1].writeToFile(filename + "Y.mtx");
    coords[2].writeToFile(filename + "Z.mtx");
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

    SCAI_ASSERT(X < NX, "X=" << X << " Y=" << Y << " Z=" << Z << " NX=" << NX << " Could not map from coordinate to index!");
    SCAI_ASSERT(Y < NY, "X=" << X << " Y=" << Y << " Z=" << Z << " NY=" << NY << " Could not map from coordinate to index!");
    SCAI_ASSERT(Z < NZ, "X=" << X << " Y=" << Y << " Z=" << Z << " NZ=" << NZ << " Could not map from coordinate to index!");
    SCAI_ASSERT(X >= 0, "X=" << X << " Y=" << Y << " Z=" << Z << " Could not map from coordinate to index!");
    SCAI_ASSERT(Y >= 0, "X=" << X << " Y=" << Y << " Z=" << Z << " Could not map from coordinate to index!");
    SCAI_ASSERT(Z >= 0, "X=" << X << " Y=" << Y << " Z=" << Z << " Could not map from coordinate to index!");

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
    SCAI_ASSERT(X < NX, "X=" << X << " Y=" << Y << " Z=" << Z << " NX=" << NX << " Could not map from coordinate to index!");
    SCAI_ASSERT(Y < NY, "X=" << X << " Y=" << Y << " Z=" << Z << " NY=" << NY << " Could not map from coordinate to index!");
    SCAI_ASSERT(Z < NZ, "X=" << X << " Y=" << Y << " Z=" << Z << " NZ=" << NZ << " Could not map from coordinate to index!");
    SCAI_ASSERT(X >= 0, "X=" << X << " Y=" << Y << " Z=" << Z << " Could not map from coordinate to index!");
    SCAI_ASSERT(Y >= 0, "X=" << X << " Y=" << Y << " Z=" << Z << " Could not map from coordinate to index!");
    SCAI_ASSERT(Z >= 0, "X=" << X << " Y=" << Y << " Z=" << Z << " Could not map from coordinate to index!");

    coordinate3D distance;

    distance.x = !((NX - 1 - X) < (X)) ? (X) : (NX - 1 - X);
    distance.y = !((NY - 1 - Y) < (Y)) ? (Y) : (NY - 1 - Y);
    distance.z = !((NZ - 1 - Z) < (Z)) ? (Z) : (NZ - 1 - Z);

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
