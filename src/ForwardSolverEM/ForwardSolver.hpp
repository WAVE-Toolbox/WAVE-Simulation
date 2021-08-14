
#pragma once
#include <scai/dmemo.hpp>
#include <scai/hmemo.hpp>

#include "../Acquisition/AcquisitionGeometry.hpp"

#include "../Common/HostPrint.hpp"
#include "../Modelparameter/Modelparameter.hpp"
#include "../Wavefields/Wavefields.hpp"
#include "../ForwardSolver/Derivatives/Derivatives.hpp"

#include "BoundaryCondition/ABS.hpp"
#include "BoundaryCondition/CPML.hpp"
#include "BoundaryCondition/FreeSurface.hpp"

namespace KITGPI
{

    //! \brief ForwardSolver namespace
    namespace ForwardSolver
    {

        //! \brief Abstract class for forward solver
        template <typename ValueType>
        class ForwardSolverEM
        {
          public:
            //! \brief Declare ForwardSolverEM pointer
            typedef std::shared_ptr<ForwardSolverEM<ValueType>> ForwardSolverPtr;

            //! Default constructor
            ForwardSolverEM() : useFreeSurface(false), useDampingBoundary(false), useConvPML(false){};

            //! Default destructor
            ~ForwardSolverEM(){};

            /*! \brief Running the foward solver
             *
             * Start the forward solver as defined by the given parameters
             *
             \param receiver Configuration of the receivers
             \param sources Configuration of the sources
             \param model Configuration of the modelparameter
             \param wavefield Wavefields for the modelling
             \param derivatives Derivations matrices to calculate the spatial derivatives
             \param NT Total number of time steps
             \param DT Temporal Sampling intervall in seconds
             */
            virtual void run(Acquisition::AcquisitionGeometry<ValueType> &receiver, Acquisition::AcquisitionGeometry<ValueType> const &sources, Modelparameter::Modelparameter<ValueType> const &model, Wavefields::Wavefields<ValueType> &wavefield, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives, scai::IndexType t) = 0;
            virtual void runAdjoint(Acquisition::AcquisitionGeometry<ValueType> &receiver, Acquisition::AcquisitionGeometry<ValueType> const &sources, Modelparameter::Modelparameter<ValueType> const &model, Wavefields::Wavefields<ValueType> &wavefield, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives, scai::IndexType t) = 0;

            /*! \brief Reset PML memory variables if PML is used 
            *
            */
            virtual void resetCPML() = 0;

            /*! \brief Sets Scaling facrors for horizontal updates at the Free Surface
            *
            *
            \param model model parameter object
            */
            virtual void prepareForModelling(Modelparameter::Modelparameter<ValueType> const &model, ValueType DT) = 0;

            ValueType estimateBoundaryMemory(Configuration::Configuration const &config, scai::dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> const &modelCoordinates, BoundaryCondition::ABSEM<ValueType> &DampingBoundary, BoundaryCondition::CPMLEM<ValueType> &ConvPML);

            virtual ValueType estimateMemory(Configuration::Configuration const &config, scai::dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> const &modelCoordinates) = 0;

            /*! \brief Initialitation of the boundary conditions
             *
             *
             \param config Configuration
             \param modelCoordinates Coordinate class, which eg. maps 3D coordinates to 1D model indices
             \param derivatives Derivatives matrices
             \param dist Distribution of the wave fields
             \param ctx Context
             */
            virtual void prepareBoundaryConditions(Configuration::Configuration const &config, Acquisition::Coordinates<ValueType> const &modelCoordinates, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> &derivatives, scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx) = 0;

            void prepareBoundaries(Configuration::Configuration const &config, Acquisition::Coordinates<ValueType> const &modelCoordinates, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> &derivatives, scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, KITGPI::ForwardSolver::BoundaryCondition::FreeSurfaceEM<ValueType> &FreeSurface, BoundaryCondition::ABSEM<ValueType> &DampingBoundary, BoundaryCondition::CPMLEM<ValueType> &ConvPML);

            /*! \brief Initialitation of the ForwardSolver
             *
             *
             \param config Configuration
             \param derivatives Derivatives matrices
             \param wavefield Wavefields for the modelling
             \param ctx Context
             */
            virtual void initForwardSolver(Configuration::Configuration const &config, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> &derivatives, Wavefields::Wavefields<ValueType> &wavefield, Modelparameter::Modelparameter<ValueType> const &model, Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::hmemo::ContextPtr ctx, ValueType DT) = 0;

          protected:
            scai::IndexType useFreeSurface; //!< Indicator which free surface is in use
            bool useDampingBoundary;        //!< Bool if damping boundary is in use
            bool useConvPML;                //!< Bool if CPMLEM is in use            

            scai::lama::DenseVector<ValueType> getAveragedCinv(scai::lama::Vector<ValueType> const &vecAvDielectricPermittivity, scai::lama::Vector<ValueType> const &vecAvElectricConductivity, ValueType DT);
            scai::lama::DenseVector<ValueType> getAveragedCa(scai::lama::Vector<ValueType> const &vecAvDielectricPermittivity, scai::lama::Vector<ValueType> const &vecAvElectricConductivity, ValueType DT);
            scai::lama::DenseVector<ValueType> getAveragedCb(scai::lama::Vector<ValueType> const &vecAvDielectricPermittivity, scai::lama::Vector<ValueType> const &vecAvElectricConductivity, ValueType DT);
            ValueType getCc(ValueType tauElectricDisplacement, ValueType DT);
            scai::lama::DenseVector<ValueType> getAveragedCd(scai::lama::Vector<ValueType> const &vecAvDielectricPermittivitystatic, scai::lama::Vector<ValueType> const &vecAvTauDielectricPermittivity, scai::IndexType numRelaxationMechanisms, ValueType tauElectricDisplacement, ValueType DT);
                  
            scai::lama::DenseVector<ValueType> update;
            scai::lama::DenseVector<ValueType> update_temp; 
            ValueType DT_temp; //!< Vector storing DT.
            scai::lama::DenseVector<ValueType> CaAverageX; //!< Vector storing averaged coefficient Ca in x-direction.
            scai::lama::DenseVector<ValueType> CaAverageY; //!< Vector storing averaged coefficient Ca in y-direction.
            scai::lama::DenseVector<ValueType> CaAverageZ; //!< Vector storing averaged coefficient Ca in z-direction.
            scai::lama::DenseVector<ValueType> CbAverageX; //!< Vector storing averaged coefficient Cb in x-direction.
            scai::lama::DenseVector<ValueType> CbAverageY; //!< Vector storing averaged coefficient Cb in y-direction.
            scai::lama::DenseVector<ValueType> CbAverageZ; //!< Vector storing averaged coefficient Cb in z-direction.            
            ValueType Cc; //!< Vector storing averaged coefficient Cc.
            scai::lama::DenseVector<ValueType> CdAverageX; //!< Vector storing averaged coefficient Cd in x-direction.
            scai::lama::DenseVector<ValueType> CdAverageY; //!< Vector storing averaged coefficient Cd in y-direction.
            scai::lama::DenseVector<ValueType> CdAverageZ; //!< Vector storing averaged coefficient Cd in z-direction.
            
        };
    } /* end namespace ForwardSolver */
} /* end namespace KITGPI */
