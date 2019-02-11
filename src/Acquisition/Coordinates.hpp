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
        template <typename ValueType>
        class Coordinates
        {

          public:
            //! \brief Default constructor
            Coordinates(){};

            //! Destructor, releases all allocated resources.
            ~Coordinates(){};
            // constructor for regular grid
            Coordinates(scai::IndexType NX, scai::IndexType NY, scai::IndexType NZ, ValueType DH);
            Coordinates(scai::IndexType NX, scai::IndexType NY, scai::IndexType NZ, ValueType DH,scai::dmemo::DistributionPtr dist,scai::hmemo::ContextPtr ctx);
            
            ValueType getDH() const;
            scai::IndexType getNX() const;
            scai::IndexType getNY() const;
            scai::IndexType getNZ() const;
            std::vector<scai::lama::DenseVector<ValueType>> getCoordinates() const;

            
            // Coordinate --> Index:
            // Interfaces 3-D
            scai::IndexType coordinate2index(coordinate3D coordinate) const;
            scai::IndexType coordinate2index(scai::IndexType X, scai::IndexType Y, scai::IndexType Z) const;

            // Index --> Coordinate:
            coordinate3D index2coordinate(scai::IndexType index) const;

            coordinate3D edgeDistance(coordinate3D coordinate) const;

            bool locatedOnSurface(scai::IndexType index) const;

          private:
            scai::IndexType NX;
            scai::IndexType NY;
            scai::IndexType NZ;
            ValueType DH;
             
            std::vector<scai::lama::DenseVector<ValueType>> coordinateVector;
            
            // Coordinate --> Index:
            scai::IndexType map3Dcoordinate2index(scai::IndexType X, scai::IndexType Y, scai::IndexType Z) const;

            // Index --> Coordinate:
            coordinate3D mapIndex2coordinate(scai::IndexType index) const;

            coordinate3D estimateDistanceToEdges3D(scai::IndexType X, scai::IndexType Y, scai::IndexType Z) const;
        };
    }
}
