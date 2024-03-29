#pragma once

#include "../Configuration/Configuration.hpp"
#include <scai/lama.hpp>

using namespace scai;

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

            friend std::ostream& operator<< (std::ostream &out, const  coordinate3D &coordinate3D) {
                out << "(x,y,z) = (" << coordinate3D.x << ", " << coordinate3D.y << ", " << coordinate3D.z << ")";
                return out;
            }

            bool operator==(const coordinate3D& other) const
            {
                if((other.x==x) && (other.y==y) && (other.z==z))
                    return (true);
                else
                    return (false);
            }

            bool operator!=(const coordinate3D& other) const
            {
                    return (!(*this==other));
            }

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
            //! \brief Default constructor
            Coordinates(): NX(0), NY(0), NZ(0), DH(0.0){};

            //! Destructor, releases all allocated resources.
            ~Coordinates(){};

            Coordinates(Configuration::Configuration const &config, IndexType DHInversion = 1, ValueType NXPerShot = 0);
            
            void init(Configuration::Configuration const &config, IndexType DHInversion = 1, ValueType NXPerShot = 0);
            // constructor for regular grid
            Coordinates(IndexType NX, IndexType NY, IndexType NZ, ValueType DH);

            Coordinates(IndexType nx, IndexType ny, IndexType nz, ValueType dh, std::vector<IndexType> &dhFactor, std::vector<int> &interface);
            
            ValueType getDH() const;
            ValueType getX0() const;
            ValueType getY0() const;
            ValueType getZ0() const;
            IndexType getNX() const;
            IndexType getNY() const;
            IndexType getNZ() const;
            IndexType getNGridpoints() const;
            IndexType getNGridpoints(IndexType layer) const;
            ValueType getDH(coordinate3D coordinate) const;
            ValueType getDH(IndexType layer) const;
            IndexType getLayer(coordinate3D coordinate) const;
            IndexType getNumLayers() const;
            IndexType getDHFactor(coordinate3D coordinate) const;
            IndexType getDHFactor(IndexType layer) const;
            std::vector<int> getInterfaceVec() const;

            std::vector<scai::lama::DenseVector<ValueType>> getCoordinates(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx) const;
            void writeCoordinates(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, std::string filename, IndexType fileFormat) const;

            // Coordinate --> Index:
            // Interfaces 3-D
            IndexType coordinate2index(coordinate3D coordinate) const;
            IndexType coordinate2index(IndexType X, IndexType Y, IndexType Z) const;

            // Index --> Coordinate:
            coordinate3D index2coordinate(IndexType index) const;

            coordinate3D edgeDistance(coordinate3D coordinate) const;

            bool locatedOnSurface(IndexType index) const;

            IndexType distToInterface(IndexType Y) const;

            bool locatedOnInterface(IndexType yCoordinate) const;
            bool locatedOnInterface(coordinate3D coordinate) const;

            int getTransition(coordinate3D coordinate) const;
            int getTransition(IndexType yCoordinate) const;
            bool isVariable() const { return (VariableGrid); };

          private:
            IndexType NX; //!< Number of gridpoints in x direction
            IndexType NY; //!< Number of gridpoints in y direction
            IndexType NZ; //!< Number of gridpoints in z direction
            ValueType DH;       //!< Gridspacing in m
            ValueType x0 = 0.0;       //!< the starting x coordinate in m
            ValueType y0 = 0.0;       //!< the starting y coordinate in m
            ValueType z0 = 0.0;       //!< the starting z coordinate in m

            std::vector<IndexType> dhFactor; //!< factors of the gridspacing relativ to the finest gridspacing in the model (must be 3^n)
            std::vector<int> interface;            //!< interfaces of the variable grid relative to the aquidistant grid with the gridspacing 1*DH

            std::vector<ValueType> varDH;       //!< DH values per layer
            std::vector<IndexType> varNX; //!< Number of gridpoints in x direction per layer
            std::vector<IndexType> varNY; //!< Number of gridpoints in y direction per layer
            std::vector<IndexType> varNZ; //!< Number of gridpoints in z direction per layer

            IndexType numLayers;               //!< Number of layers of the variable grid
            std::vector<IndexType> layerStart; //!< start position of each layer
            std::vector<IndexType> layerEnd;   //!< end position of each layer
            std::vector<int> transition;             //!< transition of the interfaces 1= fine to coarse, 0 = coarse to fine

            IndexType nGridpoints;                      //!< total number of gridpoints
            std::vector<IndexType> nGridpointsPerLayer; //!< number of gridpoints per layer

            void init();
            void init(std::vector<IndexType> &dhFactor, std::vector<int> &interface);

            // Coordinate --> Index:
            IndexType map3Dcoordinate2index(IndexType X, IndexType Y, IndexType Z) const;

            // Index --> Coordinate:
            coordinate3D mapIndex2coordinate(IndexType index) const;

            coordinate3D estimateDistanceToEdges3D(IndexType X, IndexType Y, IndexType Z) const;

            bool VariableGrid;
        };
    }
}
