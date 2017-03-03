#include "Coordinates.hpp"
using namespace scai;

/* ------- */
/* Mapping */
/* ------- */

/*! \brief Returns bool if given coordinate is located on the surface
 *
 * This method determines if a given coordinate is located on the surface of the modelling domain.
 \param coordinate 1-D coordinate
 \param NX Number of grid points in X
 \param NY Number of grid points in Y
 \param NZ Number of grid points in Z
 *
 */
template <typename ValueType>
bool KITGPI::Acquisition::Coordinates<ValueType>::locatedOnSurface(IndexType coordinate, IndexType NX, IndexType NY, IndexType /*NZ*/)
{
    coordinate3D result;
    result = map3Dindex2coordinate(coordinate, NX, NY);
    if (result.y == 0) {
        return (true);
    } else {
        return (false);
    }
}

/*! \brief General mapping from 1-D coordinate to 3-D coordinate
 *
 \param coordinate 1-D coordinate
 \param NX Number of grid points in X
 \param NY Number of grid points in Y
 *
 */
template <typename ValueType>
KITGPI::Acquisition::coordinate3D KITGPI::Acquisition::Coordinates<ValueType>::map3Dindex2coordinate(IndexType coordinate, IndexType NX, IndexType NY)
{
    coordinate3D result;

    result.z = IndexType(coordinate) / (NX * NY);
    coordinate -= result.z * (NX * NY);

    result.y = IndexType(coordinate) / (NX);
    coordinate -= result.y * (NX);

    result.x = coordinate;

    return (result);
}

/*! \brief General mapping from 1-D coordinate to 3-D coordinate
 *
 * Maps a 1-D coordinate back into 3-D coordinates.
 * The 3-D grid starts at 0 and runs to (NX-1), (NY-1) or (NZ-1).
 *
 \param coordinate 1-D coordinate
 \param NX Number of grid points in X
 \param NY Number of grid points in Y
 \param NZ Number of grid points in Z
 *
 */
template <typename ValueType>
KITGPI::Acquisition::coordinate3D KITGPI::Acquisition::Coordinates<ValueType>::index2coordinate(IndexType coordinate, IndexType NX, IndexType NY, IndexType /*NZ*/)
{
    return (map3Dindex2coordinate(coordinate, NX, NY));
}

/*! \brief General mapping from 3-D coordinates to 1-D coordinate
 *
 \param X 3-D coordinate in X (Horizontal 1)
 \param Y 3-D coordinate in Y (Depth)
 \param Z 3-D coordinate in Z (Horizontal 2)
 \param NX Total number of grid points in X
 \param NY Total number of grid points in Y
 \param NZ Total number of grid points in Z
 */
template <typename ValueType>
IndexType KITGPI::Acquisition::Coordinates<ValueType>::map3Dcoordinate2index(IndexType X, IndexType Y, IndexType Z, IndexType NX, IndexType NY, IndexType NZ)
{

    SCAI_ASSERT(Z < NZ, "Could not map from coordinate to index!");
    SCAI_ASSERT(X < NX, "Could not map from coordinate to index!");
    SCAI_ASSERT(Y < NY, "Could not map from coordinate to index!");
    SCAI_ASSERT(Z >= 0, "Could not map from coordinate to index!");
    SCAI_ASSERT(Y >= 0, "Could not map from coordinate to index!");
    SCAI_ASSERT(X >= 0, "Could not map from coordinate to index!");

    return ((X) + (Y)*NX + (Z)*NX * NY);
}

/* ---------- */
/* Interfaces */
/* ---------- */

/*! \brief Convert 3-D coordinates to 1-D coordinates
 *
 * This method returns the 1-D coordinate of 3-D coordinates.
 \param X 3-D coordinate in X (Horizontal 1)
 \param Y 3-D coordinate in Y (Depth)
 \param Z 3-D coordinate in Z (Horizontal 2)
 \param NX Total number of grid points in X
 \param NY Total number of grid points in Y
 \param NZ Total number of grid points in Z
 \return 1-D coordinate
 */
template <typename ValueType>
IndexType KITGPI::Acquisition::Coordinates<ValueType>::coordinate2index(IndexType X, IndexType Y, IndexType Z, IndexType NX, IndexType NY, IndexType NZ)
{
    return (map3Dcoordinate2index(X, Y, Z, NX, NY, NZ));
}

/*! \brief Convert 3-D coordinates to 1-D coordinates
 *
 * This method returns the 1-D coordinate of 3-D coordinates.
 \param coordinate as a coordinate3D struct
 \param NX Total number of grid points in X (Horizontal 1)
 \param NY Total number of grid points in Y (Depth)
 \param NZ Total number of grid points in Z (Horizontal 2)
 \return 1-D coordinate
 */
template <typename ValueType>
IndexType KITGPI::Acquisition::Coordinates<ValueType>::coordinate2index(coordinate3D coordinate, IndexType NX, IndexType NY, IndexType NZ)
{
    return (map3Dcoordinate2index(coordinate.x, coordinate.y, coordinate.z, NX, NY, NZ));
}

/*! \brief Determination of local coordinates based on given global coordinates
 *
 * Calculate the number of coordinates within the local processing unit as well as
 * the coordinates of the local coordinates.
 *
 \param coordinatesglobal DenseVector with global coordinates
 \param localIndices DenseVector with local coordinates
 \param dist Distribution of global grid
 */
template <typename ValueType>
void KITGPI::Acquisition::Coordinates<ValueType>::Global2Local(scai::lama::Vector const &coordinatesglobal, scai::hmemo::HArray<IndexType> &localIndices, scai::dmemo::DistributionPtr dist) const
{

    IndexType n_global = coordinatesglobal.size(); // Number of global entries

    IndexType coordinatetemp_int;
    scai::lama::Scalar coordinatetemp_scalar = 0;

    IndexType i = 0;
    for (IndexType n = 0; n < n_global; n++) {

        coordinatetemp_scalar = coordinatesglobal.getValue(n);
        coordinatetemp_int = coordinatetemp_scalar.getValue<IndexType>();

        if (dist->isLocal(coordinatetemp_int)) {
            i++;
        }
    }

    /* Determine coordinates of local receivers in the global coordinate vector */
    localIndices.resize(i);
    hmemo::WriteAccess<IndexType> write_localIndices(localIndices);
    i = 0;
    for (IndexType n = 0; n < n_global; n++) {

        coordinatetemp_scalar = coordinatesglobal.getValue(n);
        coordinatetemp_int = coordinatetemp_scalar.getValue<IndexType>();
        if (dist->isLocal(coordinatetemp_int)) {
            write_localIndices[i] = n;
            i++;
        }
    }
}

/*! \brief Calculation of distance to boundaries of the modelling domain
 *
 * This method calculates the distance of a given coordinate to the boundaries of the modelling domain.
 \param X Coordinate X
 \param Y Coordinate Y
 \param Z Coordinate Z
 \param NX Total number of grid points in X
 \param NY Total number of grid points in Y
 \param NZ Total number of grid points in Z
 */
template <typename ValueType>
KITGPI::Acquisition::coordinate3D KITGPI::Acquisition::Coordinates<ValueType>::estimateDistanceToEdges3D(IndexType X, IndexType Y, IndexType Z, IndexType NX, IndexType NY, IndexType NZ)
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
    distance.z = !((NX - 1 - Z) < (Z)) ? (Z) : (NZ - 1 - Z);

    return (distance);
}

/*! \brief Determination of distance to boundaries
 *
 * This method calculates the distance of a given coordinate3D to the boundaries of the modelling domain.
 \param coordinate 3D-coordinate structs
 \param NX Total number of grid points in X
 \param NY Total number of grid points in Y
 \param NZ Total number of grid points in Z
 */
template <typename ValueType>
KITGPI::Acquisition::coordinate3D KITGPI::Acquisition::Coordinates<ValueType>::edgeDistance(coordinate3D coordinate, IndexType NX, IndexType NY, IndexType NZ)
{
    return (estimateDistanceToEdges3D(coordinate.x, coordinate.y, coordinate.z, NX, NY, NZ));
}

template class KITGPI::Acquisition::Coordinates<double>;
template class KITGPI::Acquisition::Coordinates<float>;
