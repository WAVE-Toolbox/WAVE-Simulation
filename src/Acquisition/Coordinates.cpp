#include "Coordinates.hpp"
using namespace scai;

/* ------- */
/* Mapping */
/* ------- */

/*! \brief constructor for regular grid
 *
 \param NumX Number of grid points in X
 \param NumY Number of grid points in Y
 \param NumZ Number of grid points in Z
 *
 */
KITGPI::Acquisition::Coordinates::Coordinates(scai::IndexType NumX, scai::IndexType NumY, scai::IndexType NumZ)
{
    NX = NumX;
    NY = NumY;
    NZ = NumZ;
}

/*! \brief Returns bool if given coordinate is located on the surface
 *
 * This method determines if a given coordinate is located on the surface of the modelling domain.
 \param coordinate 1-D coordinate
 \param NX Number of grid points in X
 \param NY Number of grid points in Y
 \param NZ Number of grid points in Z
 *
 */
bool KITGPI::Acquisition::Coordinates::locatedOnSurface(IndexType coordinate)
{
    coordinate3D result;
    result = map3Dindex2coordinate(coordinate);
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
KITGPI::Acquisition::coordinate3D KITGPI::Acquisition::Coordinates::map3Dindex2coordinate(IndexType coordinate)
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
KITGPI::Acquisition::coordinate3D KITGPI::Acquisition::Coordinates::index2coordinate(IndexType coordinate)
{
    return (map3Dindex2coordinate(coordinate));
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
IndexType KITGPI::Acquisition::Coordinates::map3Dcoordinate2index(IndexType X, IndexType Y, IndexType Z)
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
IndexType KITGPI::Acquisition::Coordinates::coordinate2index(IndexType X, IndexType Y, IndexType Z)
{
    return (map3Dcoordinate2index(X, Y, Z));
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
IndexType KITGPI::Acquisition::Coordinates::coordinate2index(coordinate3D coordinate)
{
    return (map3Dcoordinate2index(coordinate.x, coordinate.y, coordinate.z));
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
KITGPI::Acquisition::coordinate3D KITGPI::Acquisition::Coordinates::estimateDistanceToEdges3D(IndexType X, IndexType Y, IndexType Z)
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
 \param NX Total number of grid points in X
 \param NY Total number of grid points in Y
 \param NZ Total number of grid points in Z
 */
KITGPI::Acquisition::coordinate3D KITGPI::Acquisition::Coordinates::edgeDistance(coordinate3D coordinate)
{
    return (estimateDistanceToEdges3D(coordinate.x, coordinate.y, coordinate.z));
}
