#pragma once

#include "../../Acquisition/Coordinates.hpp"
#include "../../Common/HostPrint.hpp"
#include <scai/dmemo.hpp>
#include <scai/hmemo.hpp>
#include <scai/lama.hpp>

namespace KITGPI
{

    namespace ForwardSolver
    {

        namespace BoundaryCondition
        {

            //! \brief Abstract class for the calculation and application of cpml boundaries
            template <typename ValueType>
            class CPMLEM
            {
              public:
                //! \brief Default constructor
                CPMLEM(){};

                //! \brief Default destructor
                ~CPMLEM(){};

                virtual ValueType estimateMemory(scai::IndexType BoundaryWidth, scai::IndexType useFreeSurface, scai::dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> const &modelCoordinates) = 0;

                //! init CPMLEM coefficient vectors and CPMLEM memory variables
                virtual void init(scai::dmemo::DistributionPtr const dist, scai::hmemo::ContextPtr const ctx, Acquisition::Coordinates<ValueType> const &modelCoordinates, ValueType const DT, scai::IndexType const BoundaryWidth, ValueType const NPower, ValueType const CenterFrequencyCPML, ValueType const VMaxCPML, scai::IndexType const useFreeSurface) = 0;

              protected:
                typedef typename scai::lama::SparseVector<ValueType> VectorType;

                void resetVector(scai::lama::Vector<ValueType> &vector);

                void initVector(scai::lama::Vector<ValueType> &vector, scai::hmemo::ContextPtr const ctx, scai::dmemo::DistributionPtr const dist);

                void calcCoeffCPMLEM(std::vector<ValueType> &a, std::vector<ValueType> &b, ValueType const NPower, ValueType const CenterFrequencyCPML, ValueType const VMaxCPML, ValueType const DT, ValueType const DH, bool const shiftGrid = false);

                /*inline*/ void applyCPMLEM(scai::lama::DenseVector<ValueType> &Vec, VectorType &Psi, VectorType const &a, VectorType const &b);

                VectorType temp; //!< temporary vector for pml application

                bool active; //!< Bool if CPMLEM is active
            };
        } /* end namespace BoundaryCondition  */
    }     /* end namespace ForwardSolver */
} /* end namespace KITGPI */
