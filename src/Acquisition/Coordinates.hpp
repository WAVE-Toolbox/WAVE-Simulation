#pragma once

#include "../Configuration/Configuration.hpp"
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

            Coordinates(Configuration::Configuration const &config);
            // constructor for regular grid
            Coordinates(scai::IndexType NX, scai::IndexType NY, scai::IndexType NZ, ValueType DH);

            Coordinates(scai::IndexType nx, scai::IndexType ny, scai::IndexType nz, ValueType dh, std::vector<scai::IndexType> &dhFactor, std::vector<scai::IndexType> &interface);

            ValueType getDH() const;
            scai::IndexType getNX() const;
            scai::IndexType getNY() const;
            scai::IndexType getNZ() const;
            scai::IndexType getNGridpoints() const;
            ValueType getDH(coordinate3D coordinate) const;
            ValueType getDH(scai::IndexType layer) const;
            scai::IndexType getLayer(coordinate3D coordinate) const;
            scai::IndexType getNumLayers() const;
            scai::IndexType getDHFactor(coordinate3D coordinate) const;
            scai::IndexType getDHFactor(scai::IndexType layer) const;

            std::vector<scai::lama::DenseVector<ValueType>> getCoordinates(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx) const;
            void writeCoordinates(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, std::string filename) const;

            // Coordinate --> Index:
            // Interfaces 3-D
            scai::IndexType coordinate2index(coordinate3D coordinate) const;
            scai::IndexType coordinate2index(scai::IndexType X, scai::IndexType Y, scai::IndexType Z) const;

            // Index --> Coordinate:
            coordinate3D index2coordinate(scai::IndexType index) const;

            coordinate3D edgeDistance(coordinate3D coordinate) const;

            bool locatedOnSurface(scai::IndexType index) const;

            scai::IndexType distToInterface(scai::IndexType Y) const;
            bool locatedOnInterface(coordinate3D coordinate) const;
            bool getTransition(coordinate3D coordinate) const;
            bool isVariable() const { return (VariableGrid); };

          private:
            scai::IndexType NX; //!< Number of gridpoints in x direction
            scai::IndexType NY; //!< Number of gridpoints in y direction
            scai::IndexType NZ; //!< Number of gridpoints in z direction
            ValueType DH;       //!< Gridspacing in m

            std::vector<scai::IndexType> dhFactor;  //!< factors of the gridspacing relativ to the finest gridspacing in the model (must be 3^n)
            std::vector<scai::IndexType> interface; //!< interfaces of the variable grid relative to the aquidistant grid with the gridspacing 1*DH

            std::vector<ValueType> varDH;       //!< DH values per layer
            std::vector<scai::IndexType> varNX; //!< Number of gridpoints in x direction per layer
            std::vector<scai::IndexType> varNY; //!< Number of gridpoints in y direction per layer
            std::vector<scai::IndexType> varNZ; //!< Number of gridpoints in z direction per layer

            scai::IndexType numLayers;               //!< Number of layers of the variable grid
            std::vector<scai::IndexType> layerStart; //!< start position of each layer
            std::vector<scai::IndexType> layerEnd;   //!< end position of each layer
            std::vector<scai::IndexType> transition; //!< transition of the interfaces 1= fine to coarse, 0 = coarse to fine

            scai::IndexType nGridpoints;                      //!< total number of gridpoints
            std::vector<scai::IndexType> nGridpointsPerLayer; //!< number of gridpoints per layer

            void init();
            void init(std::vector<scai::IndexType> &dhFactor, std::vector<scai::IndexType> &interface);

            // Coordinate --> Index:
            scai::IndexType map3Dcoordinate2index(scai::IndexType X, scai::IndexType Y, scai::IndexType Z) const;

            // Index --> Coordinate:
            coordinate3D mapIndex2coordinate(scai::IndexType index) const;

            coordinate3D estimateDistanceToEdges3D(scai::IndexType X, scai::IndexType Y, scai::IndexType Z) const;

            bool VariableGrid;
        };
    }
}
