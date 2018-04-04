#pragma once

#include <scai/lama.hpp>

namespace KITGPI
{

    namespace Acquisition
    {

        /*! \brief Struct to save 3-D coordinates
         *
         * This struct saves coordinates as 3-D components.
         */
        struct coordinate3D {
            scai::IndexType x; //!< x Position in X-direction in grid points (Horizontal 1)
            scai::IndexType y; //!< y Position in Y-direction in grid points (Depth)
            scai::IndexType z; //!< z Position in Z-direction in grid points (Horizontal 2)

            /*! \brief Return the minimum of all three values */
            scai::IndexType min()
            {
                scai::IndexType temp = 0;
                if (x < y) {
                    temp = x;
                } else {
                    temp = y;
                }
                if (z < temp) {
                    temp = z;
                }
                return (temp);
            }
        };

        /*! \brief This class manages the transformation of Coordinates from 3-D to 1-D and vice-versa
         */
        class Coordinates
        {

          public:
            // Coordinate --> Index:
            // Interfaces 3-D
            scai::IndexType coordinate2index(coordinate3D coordinate, scai::IndexType NX, scai::IndexType NY, scai::IndexType NZ);
            scai::IndexType coordinate2index(scai::IndexType X, scai::IndexType Y, scai::IndexType Z, scai::IndexType NX, scai::IndexType NY, scai::IndexType NZ);

            // Index --> Coordinate:
            coordinate3D index2coordinate(scai::IndexType coordinate, scai::IndexType NX, scai::IndexType NY, scai::IndexType NZ);

            coordinate3D edgeDistance(coordinate3D coordinate, scai::IndexType NX, scai::IndexType NY, scai::IndexType NZ);

            bool locatedOnSurface(scai::IndexType coordinate, scai::IndexType NX, scai::IndexType NY, scai::IndexType NZ);

            void Global2Local(scai::lama::Vector<scai::IndexType> const &coordinatesglobal, scai::hmemo::HArray<scai::IndexType> &coordinateslocal, scai::dmemo::DistributionPtr dist) const;

          private:
            // Coordinate --> Index:
            scai::IndexType map3Dcoordinate2index(scai::IndexType X, scai::IndexType Y, scai::IndexType Z, scai::IndexType NX, scai::IndexType NY, scai::IndexType NZ);

            // Index --> Coordinate:
            coordinate3D map3Dindex2coordinate(scai::IndexType coordinate, scai::IndexType NX, scai::IndexType NY);

            coordinate3D estimateDistanceToEdges3D(scai::IndexType X, scai::IndexType Y, scai::IndexType Z, scai::IndexType NX, scai::IndexType NY, scai::IndexType NZ);
        };
    }
}
