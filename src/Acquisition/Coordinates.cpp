#include "Coordinates.hpp"
using namespace scai;

/* ------- */
/* Mapping */
/* ------- */

/*! \brief constructor for regular grid
 *
 \param NX Number of grid points in X
 \param NY Number of grid points in Y
 \param NZ Number of grid points in Z
 *
 */
template <typename ValueType>
KITGPI::Acquisition::Coordinates<ValueType>::Coordinates(scai::IndexType NX, scai::IndexType NY, scai::IndexType NZ, ValueType DH) : NX(NX), NY(NY), NZ(NZ), DH(DH)
{
    SCAI_ASSERT_ERROR(NX > 0, "NX<=0");
    SCAI_ASSERT_ERROR(NY > 0, "NY<=0");
    SCAI_ASSERT_ERROR(NZ > 0, "NZ<=0");
}

/*! \brief constructor for regular grid
 *
 \param NX Number of grid points in X
 \param NY Number of grid points in Y
 \param NZ Number of grid points in Z
 *
 */
template <typename ValueType>
KITGPI::Acquisition::Coordinates<ValueType>::Coordinates(scai::IndexType NX, scai::IndexType NY, scai::IndexType NZ,  ValueType DH, scai::dmemo::DistributionPtr dist,scai::hmemo::ContextPtr ctx) : NX(NX), NY(NY), NZ(NZ), DH(DH)
{
    SCAI_ASSERT_ERROR(NX > 0, "NX<=0");
    SCAI_ASSERT_ERROR(NY > 0, "NY<=0");
    SCAI_ASSERT_ERROR(NZ > 0, "NZ<=0");
    
    hmemo::HArray<IndexType> ownedIndexes; // all (global) points owned by this process
    dist->getOwnedIndexes(ownedIndexes);

    lama::VectorAssembly<IndexType> xAssembly;
    xAssembly.reserve(ownedIndexes.size());
    lama::VectorAssembly<IndexType> yAssembly;
    yAssembly.reserve(ownedIndexes.size());
    lama::VectorAssembly<IndexType> zAssembly;
    zAssembly.reserve(ownedIndexes.size()); 
    

    for (IndexType ownedIndex : hmemo::hostReadAccess(ownedIndexes)) {
    coordinate3D coordinate = index2coordinate(ownedIndex);
    xAssembly.push(ownedIndex,coordinate.x);
    yAssembly.push(ownedIndex,coordinate.y);
    zAssembly.push(ownedIndex,coordinate.z);
    }
    
    xCoordinates.allocate(dist);
    xCoordinates.fillFromAssembly(xAssembly);
    xCoordinates.setContextPtr(ctx);

    
    xCoordinates.writeToFile("xcoord.mtx");
    
    yCoordinates.allocate(dist);
    yCoordinates.fillFromAssembly(yAssembly);
    yCoordinates.setContextPtr(ctx);

    yCoordinates.writeToFile("ycoord.mtx");
    
    zCoordinates.allocate(dist);
    zCoordinates.fillFromAssembly(zAssembly);
    zCoordinates.setContextPtr(ctx);

    
    zCoordinates.writeToFile("zcoord.mtx");
    
    useCoordinateVectors=true;
    
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

/*! \brief getter function for x coordinates
 *
 *
 */
template <typename ValueType>
scai::lama::DenseVector<IndexType> KITGPI::Acquisition::Coordinates<ValueType>::getXCoordinates() const
{
    return (xCoordinates);
}

/*! \brief getter function for y coordinates
 *
 *
 */
template <typename ValueType>
scai::lama::DenseVector<IndexType> KITGPI::Acquisition::Coordinates<ValueType>::getYCoordinates() const
{
    return (yCoordinates);
}

/*! \brief getter function for z coordinates
 *
 *
 */
template <typename ValueType>
scai::lama::DenseVector<IndexType> KITGPI::Acquisition::Coordinates<ValueType>::getZCoordinates() const
{
    return (zCoordinates);
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
    coordinate3D result;

    
//     if (useCoordinateVectors) {
//         if (index==10){
//                std::cout <<"message map true: xcoord " << index << " " <<  xCoordinates.getValue(index) << " " << std::endl;
//                std::cout <<"message map true: xcoord " << index << " " <<  xCoordinates.getValue(index) << " " << std::endl;
//                std::cout <<"message map true: xcoord " << index << " " <<  xCoordinates.getValue(index) << " " << std::endl;
//                std::cout <<"message map true: xcoord " << index << " " <<  xCoordinates.getValue(index) << " " << std::endl;
//                std::cout <<"message map true: xcoord " << index << " " <<  xCoordinates.getValue(index) << " " << std::endl;
//                std::cout <<"message map true: xcoord " << index << " " <<  xCoordinates.getValue(index) << " " << std::endl;
//                std::cout <<"message map true: xcoord " << index << " " <<  xCoordinates.getValue(index) << " " << std::endl;
//                std::cout <<"message map true: xcoord " << index << " " <<  xCoordinates.getValue(index) << " " << std::endl;     
//         }
//         result.x = xCoordinates.getValue(index);
//         result.y = yCoordinates.getValue(index);
//         result.z = zCoordinates.getValue(index);
//     }else {
    
    result.y = IndexType(index / (NX * NZ));
    index -= result.y * (NX * NZ);

    result.z = IndexType(index / (NX));
    index -= result.z * (NX);

    result.x = index;
//     }
    
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

    SCAI_ASSERT(Z < NZ, "Could not map from coordinate to index!");
    SCAI_ASSERT(X < NX, "Could not map from coordinate to index!");
    SCAI_ASSERT(Y < NY, "Could not map from coordinate to index!");
    SCAI_ASSERT(Z >= 0, "Could not map from coordinate to index!");
    SCAI_ASSERT(Y >= 0, "Could not map from coordinate to index!");
    SCAI_ASSERT(X >= 0, "Could not map from coordinate to index!");

    return ((X) + (Z)*NX + (Y)*NX * NZ);
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

    SCAI_ASSERT(Z < NZ, "No valid argument!");
    SCAI_ASSERT(X < NX, "No valid argument!");
    SCAI_ASSERT(Y < NY, "No valid argument!");
    SCAI_ASSERT(Z >= 0, "No valid argument!");
    SCAI_ASSERT(Y >= 0, "No valid argument!");
    SCAI_ASSERT(X >= 0, "No valid argument!");

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
