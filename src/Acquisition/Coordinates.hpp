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
            IndexType x; //!< x Position in X-direction in grid points (Horizontal 1)
            IndexType y; //!< y Position in Y-direction in grid points (Depth)
            IndexType z; //!< z Position in Z-direction in grid points (Horizontal 2)

            /*! \brief Return the minimum of all three values */
            IndexType min()
            {
                IndexType temp = 0;
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
        template <typename ValueType>
        class Coordinates
        {

          public:
            // Coordinate --> Index:
            // Interfaces 3-D
            IndexType coordinate2index(coordinate3D coordinate, IndexType NX, IndexType NY, IndexType NZ);
            IndexType coordinate2index(IndexType X, IndexType Y, IndexType Z, IndexType NX, IndexType NY, IndexType NZ);

            // Index --> Coordinate:
            coordinate3D index2coordinate(IndexType coordinate, IndexType NX, IndexType NY, IndexType NZ);

            coordinate3D edgeDistance(coordinate3D coordinate, IndexType NX, IndexType NY, IndexType NZ);

            bool locatedOnSurface(IndexType coordinate, IndexType NX, IndexType NY, IndexType NZ);

            void Global2Local(scai::lama::Vector const &coordinatesglobal, scai::hmemo::HArray<IndexType> &coordinateslocal, scai::dmemo::DistributionPtr dist) const;

          private:
            // Coordinate --> Index:
            IndexType map3Dcoordinate2index(IndexType X, IndexType Y, IndexType Z, IndexType NX, IndexType NY, IndexType NZ);

            // Index --> Coordinate:
            coordinate3D map3Dindex2coordinate(IndexType coordinate, IndexType NX, IndexType NY);

            coordinate3D estimateDistanceToEdges3D(IndexType X, IndexType Y, IndexType Z, IndexType NX, IndexType NY, IndexType NZ);
        };
    }
}
