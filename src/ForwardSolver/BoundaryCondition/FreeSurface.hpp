#pragma once

#include "../../Acquisition/Coordinates.hpp"
#include "../../Common/HostPrint.hpp"
#include "../Derivatives/Derivatives.hpp"
#include "FreeSurface.hpp"

namespace KITGPI
{

    namespace ForwardSolver
    {

        //! \brief BoundaryCondition namespace
        namespace BoundaryCondition
        {

            //! \brief Abstract free surface class
            template <typename ValueType>
            class FreeSurface
            {
              public:
                //! Default constructor
                FreeSurface() : active(false){};

                //! Default destructor
                ~FreeSurface(){};

                /*! \brief Initialitation of the free surface
                 *
                 *
                 \param dist Distribution of wavefields
                 \param derivatives Derivative class
                 \param NX Number of grid points in X-Direction
                 \param NY Number of grid points in Y-Direction (Depth)
                 \param NZ Number of grid points in Z-Direction
                 \param modelCoordinates Coordinate class, which eg. maps 3D coordinates to 1D model indices
                 \param DT Temporal Sampling
                 \param DH Distance between grid points
                 */
                virtual void init(scai::dmemo::DistributionPtr dist, Derivatives::Derivatives<ValueType> &derivatives, scai::IndexType NX, scai::IndexType NY, scai::IndexType NZ, Acquisition::Coordinates<ValueType> const &modelCoordinates, ValueType DT, ValueType DH) = 0;

                /*! \brief Getter method for active bool
                 *
                 *
                 */
                virtual bool getActive() const;

              protected:
                bool active; //!< Bool if this free surface is active and initialized (==ready-to use)
            };
        } /* end namespace BoundaryCondition */
    }     /* end namespace ForwardSolver */
} /* end namespace KITGPI */
