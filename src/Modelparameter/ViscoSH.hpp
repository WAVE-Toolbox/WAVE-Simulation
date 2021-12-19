
#pragma once

#include <scai/lama.hpp>

#include <scai/lama/DenseVector.hpp>
#include <scai/lama/matrix/all.hpp>
#include <scai/lama/matutils/MatrixCreator.hpp>
#include <scai/lama/norm/L2Norm.hpp>

#include <scai/dmemo/BlockDistribution.hpp>

#include <scai/hmemo/HArray.hpp>
#include <scai/hmemo/ReadAccess.hpp>
#include <scai/hmemo/WriteAccess.hpp>

#include <scai/tracing.hpp>

#include <scai/common/Walltime.hpp>
#include <scai/logging.hpp>

#include <iostream>

#include "ModelparameterSeismic.hpp"

namespace KITGPI
{

    //! \brief Modelparameter namespace
    namespace Modelparameter
    {

        //! Class for Modelparameter for viscosh simulations (Subsurface properties)
        /*!
         This class handels the modelparameter for the viscosh finite-difference simulation.
         */
        template <typename ValueType>
        class ViscoSH : public ModelparameterSeismic<ValueType>
        {
          public:
            //! Default constructor.
            ViscoSH() { equationType = "viscosh"; };

            //! Destructor, releases all allocated resources.
            ~ViscoSH(){};

            explicit ViscoSH(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> const &modelCoordinates);
            explicit ViscoSH(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, ValueType velocityS_const, ValueType rho_const, ValueType tauS_const, IndexType numRelaxationMechanisms_in, std::vector<ValueType> relaxationFrequency_in, ValueType centerFrequencyCPML_in);
            explicit ViscoSH(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, std::string filename, scai::IndexType fileFormat);

            //! Copy Constructor.
            ViscoSH(const ViscoSH &rhs);

            ValueType estimateMemory(scai::dmemo::DistributionPtr dist) override;

            void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, ValueType velocityS_const, ValueType rho_const, ValueType tauS_const, IndexType numRelaxationMechanisms_in, std::vector<ValueType> relaxationFrequency_in, ValueType centerFrequencyCPML_in);
            void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, std::string filename, scai::IndexType fileFormat) override;
            void init(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> const &modelCoordinates) override;

            void initRelaxationMechanisms(scai::IndexType numRelaxationMechanisms_in, std::vector<ValueType> relaxationFrequency_in, ValueType centerFrequencyCPML_in);
            
            void write(std::string filename, scai::IndexType fileFormat) const override;
            void writeRockMatrixParameter(std::string filename, scai::IndexType fileFormat) override;

            std::string getEquationType() const;

            /* Getter methods for not required parameters */
            scai::lama::Vector<ValueType> const &getVelocityP() const override;
            scai::lama::Vector<ValueType> const &getPWaveModulus() override;
            scai::lama::Vector<ValueType> const &getPWaveModulus() const override;
            scai::lama::Vector<ValueType> const &getBulkModulusRockMatrix() const override;
            scai::lama::Vector<ValueType> const &getTauP() const override;
            scai::lama::Vector<ValueType> const &getTauSAverageXY() override;
            scai::lama::Vector<ValueType> const &getTauSAverageXY() const override;
            scai::lama::Vector<ValueType> const &getSWaveModulus() override;
            scai::lama::Vector<ValueType> const &getSWaveModulus() const override;
            
            void calcRockMatrixParameter(Configuration::Configuration const &config) override;
            void calcWaveModulusFromPetrophysics() override;
            void calcPetrophysicsFromWaveModulus() override;
            void calcReflectivity(Acquisition::Coordinates<ValueType> const &modelCoordinates, KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType> const &derivatives, ValueType DT) override;

            void prepareForModelling(Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, scai::dmemo::CommunicatorPtr comm) override;

            void applyThresholds(Configuration::Configuration const &config) override;

            void getModelPerShot(KITGPI::Modelparameter::Modelparameter<ValueType> &modelPerShot, Acquisition::Coordinates<ValueType> const &modelCoordinates, Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, Acquisition::coordinate3D const cutCoordinate) override;

            void minusAssign(KITGPI::Modelparameter::Modelparameter<ValueType> const &rhs) override;
            void plusAssign(KITGPI::Modelparameter::Modelparameter<ValueType> const &rhs) override;
            void assign(KITGPI::Modelparameter::Modelparameter<ValueType> const &rhs) override;

            /* Overloading Operators */
            KITGPI::Modelparameter::ViscoSH<ValueType> operator*(ValueType rhs);
            KITGPI::Modelparameter::ViscoSH<ValueType> &operator*=(ValueType const &rhs);
            KITGPI::Modelparameter::ViscoSH<ValueType> operator+(KITGPI::Modelparameter::ViscoSH<ValueType> const &rhs);
            KITGPI::Modelparameter::ViscoSH<ValueType> &operator+=(KITGPI::Modelparameter::ViscoSH<ValueType> const &rhs);
            KITGPI::Modelparameter::ViscoSH<ValueType> operator-(KITGPI::Modelparameter::ViscoSH<ValueType> const &rhs);
            KITGPI::Modelparameter::ViscoSH<ValueType> &operator-=(KITGPI::Modelparameter::ViscoSH<ValueType> const &rhs);
            KITGPI::Modelparameter::ViscoSH<ValueType> &operator=(KITGPI::Modelparameter::ViscoSH<ValueType> const &rhs);

          private:
            void init(scai::dmemo::DistributionPtr variableDist, Acquisition::Coordinates<ValueType> const &variableCoordinates, Acquisition::Coordinates<ValueType> const &regularCoordinates){COMMON_THROWEXCEPTION("variable grid is not implemented in the viscosh case")};
            void calculateAveraging() override;

            using Modelparameter<ValueType>::equationType;

            using Modelparameter<ValueType>::dirtyFlagInverseDensity;
            using Modelparameter<ValueType>::dirtyFlagSWaveModulus;
            using Modelparameter<ValueType>::dirtyFlagAveraging;
            using Modelparameter<ValueType>::sWaveModulus;
            using Modelparameter<ValueType>::density;
            using Modelparameter<ValueType>::inverseDensity;
            using Modelparameter<ValueType>::velocityS;
            using Modelparameter<ValueType>::tauS;
            using Modelparameter<ValueType>::shearModulusRockMatrix;   //!< Vector storing S-wave modulus.
            using Modelparameter<ValueType>::densityRockMatrix;        //!< Vector storing Density.
            using Modelparameter<ValueType>::porosity; 
            using Modelparameter<ValueType>::saturation;
            using Modelparameter<ValueType>::reflectivity; //!< Vector storing reflectivity. 
            using Modelparameter<ValueType>::DensityWater;
            using Modelparameter<ValueType>::DensityAir;
            using Modelparameter<ValueType>::CriticalPorosity;

            void initializeMatrices(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::CommunicatorPtr comm) override;

            void purgeMatrices() override;

            using Modelparameter<ValueType>::averageMatrixX;
            using Modelparameter<ValueType>::averageMatrixY;
            using Modelparameter<ValueType>::averageMatrixZ;
            using Modelparameter<ValueType>::averageMatrixXY;
            using Modelparameter<ValueType>::averageMatrixXZ;
            using Modelparameter<ValueType>::averageMatrixYZ;

            using Modelparameter<ValueType>::inverseDensityAverageX;
            using Modelparameter<ValueType>::inverseDensityAverageY;
            using Modelparameter<ValueType>::inverseDensityAverageZ;
            using Modelparameter<ValueType>::sWaveModulusAverageXY;
            using Modelparameter<ValueType>::sWaveModulusAverageXZ;
            using Modelparameter<ValueType>::sWaveModulusAverageYZ;
            using Modelparameter<ValueType>::relaxationFrequency;
            using Modelparameter<ValueType>::centerFrequencyCPML;
            using Modelparameter<ValueType>::numRelaxationMechanisms;
            using Modelparameter<ValueType>::tauSAverageXY;
            using Modelparameter<ValueType>::tauSAverageXZ;
            using Modelparameter<ValueType>::tauSAverageYZ;

            /* Not required parameters */
            using Modelparameter<ValueType>::pWaveModulus;
            using Modelparameter<ValueType>::bulkModulusRockMatrix;
            using Modelparameter<ValueType>::velocityP;
            using Modelparameter<ValueType>::tauP;
        };
    }
}
