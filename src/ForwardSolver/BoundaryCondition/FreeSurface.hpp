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
                 \param modelCoordinates Coordinate class, which eg. maps 3D coordinates to 1D model indices
                 \param DT Temporal Sampling
                 */
                virtual void init(scai::dmemo::DistributionPtr dist, Derivatives::Derivatives<ValueType> &derivatives, Acquisition::Coordinates<ValueType> const &modelCoordinates, ValueType DT) = 0;

                void setSurfaceZero(scai::lama::Vector<ValueType> &vector);

                /*! \brief Getter method for active bool
                 *
                 *
                 */
                virtual bool getActive() const;

              protected:
                scai::lama::SparseVector<ValueType> setZeroFreeSurface; //!< Vector, which sets everything besides the free surface to zero
                bool active;                                            //!< Bool if this free surface is active and initialized (==ready-to use)
            };
        } /* end namespace BoundaryCondition */
    }     /* end namespace ForwardSolver */
} /* end namespace KITGPI */
