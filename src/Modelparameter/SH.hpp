
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

#include "../PartitionedInOut/PartitionedInOut.hpp"
#include "Modelparameter.hpp"

namespace KITGPI
{

    //! \brief Modelparameter namespace
    namespace Modelparameter
    {

        //! Class for Modelparameter for sh simulations (Subsurface properties)
        /*!
         This class handels the modelparameter for the sh finite-difference simulation.
         */
        template <typename ValueType>
        class SH : public Modelparameter<ValueType>
        {
          public:
            //! Default constructor.
            SH() { equationType = "sh"; };

            //! Destructor, releases all allocated resources.
            ~SH(){};

            explicit SH(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist);
            explicit SH(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, ValueType velocityS_const, ValueType rho);
            explicit SH(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, std::string filename, scai::IndexType partitionedIn);

            //! Copy Constructor.
            SH(const SH &rhs);

            void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, ValueType velocityS_const, ValueType rho_const);
            void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, std::string filename, scai::IndexType partitionedIn) override;
            void init(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist) override;

            void write(std::string filename, scai::IndexType partitionedOut) const override;

            std::string getEquationType() const;

            /* Getter methods for not requiered parameters */

            scai::lama::Vector<ValueType> const &getVelocityP() const override;
            scai::lama::Vector<ValueType> const &getPWaveModulus() override;
            scai::lama::Vector<ValueType> const &getPWaveModulus() const override;
            scai::lama::Vector<ValueType> const &getTauP() const override;
            scai::lama::Vector<ValueType> const &getTauS() const override;
            scai::lama::Vector<ValueType> const &getTauSAverageXY() override;
            scai::lama::Vector<ValueType> const &getTauSAverageXY() const override;
            scai::lama::Vector<ValueType> const &getTauSAverageXZ() override;
            scai::lama::Vector<ValueType> const &getTauSAverageXZ() const override;
            scai::lama::Vector<ValueType> const &getTauSAverageYZ() override;
            scai::lama::Vector<ValueType> const &getTauSAverageYZ() const override;
            
            scai::IndexType getNumRelaxationMechanisms() const override;
            ValueType getRelaxationFrequency() const override;

            void prepareForModelling(Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, scai::dmemo::CommunicatorPtr comm) override;

            void applyThresholds(Configuration::Configuration const &config) override;

            void minusAssign(KITGPI::Modelparameter::Modelparameter<ValueType> const &rhs) override;
            void plusAssign(KITGPI::Modelparameter::Modelparameter<ValueType> const &rhs) override;
            void assign(KITGPI::Modelparameter::Modelparameter<ValueType> const &rhs) override;

            /* Overloading Operators */
            KITGPI::Modelparameter::SH<ValueType> operator*(ValueType rhs);
            KITGPI::Modelparameter::SH<ValueType> &operator*=(ValueType const &rhs);
            KITGPI::Modelparameter::SH<ValueType> operator+(KITGPI::Modelparameter::SH<ValueType> const &rhs);
            KITGPI::Modelparameter::SH<ValueType> &operator+=(KITGPI::Modelparameter::SH<ValueType> const &rhs);
            KITGPI::Modelparameter::SH<ValueType> operator-(KITGPI::Modelparameter::SH<ValueType> const &rhs);
            KITGPI::Modelparameter::SH<ValueType> &operator-=(KITGPI::Modelparameter::SH<ValueType> const &rhs);
            KITGPI::Modelparameter::SH<ValueType> &operator=(KITGPI::Modelparameter::SH<ValueType> const &rhs);

          private:
            void calculateAveraging() override;

            using Modelparameter<ValueType>::equationType;

            using Modelparameter<ValueType>::dirtyFlagInverseDensity;
            using Modelparameter<ValueType>::dirtyFlagSWaveModulus;
            using Modelparameter<ValueType>::dirtyFlagAveraging;
            using Modelparameter<ValueType>::sWaveModulus;
            using Modelparameter<ValueType>::density;
            using Modelparameter<ValueType>::inverseDensity;
            using Modelparameter<ValueType>::velocityS;

            void initializeMatrices(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::CommunicatorPtr comm) override;

            using Modelparameter<ValueType>::DensityAverageMatrixX;
            using Modelparameter<ValueType>::DensityAverageMatrixY;
            using Modelparameter<ValueType>::DensityAverageMatrixZ;
            using Modelparameter<ValueType>::sWaveModulusAverageMatrixXY;
            using Modelparameter<ValueType>::sWaveModulusAverageMatrixXZ;
            using Modelparameter<ValueType>::sWaveModulusAverageMatrixYZ;

            using Modelparameter<ValueType>::inverseDensityAverageX;
            using Modelparameter<ValueType>::inverseDensityAverageY;
            using Modelparameter<ValueType>::inverseDensityAverageZ;
            using Modelparameter<ValueType>::sWaveModulusAverageXY;
            using Modelparameter<ValueType>::sWaveModulusAverageXZ;
            using Modelparameter<ValueType>::sWaveModulusAverageYZ;

            /* Not requiered parameters */
            using Modelparameter<ValueType>::pWaveModulus;
            using Modelparameter<ValueType>::velocityP;
            using Modelparameter<ValueType>::tauP;
            using Modelparameter<ValueType>::tauS;
            using Modelparameter<ValueType>::relaxationFrequency;
            using Modelparameter<ValueType>::numRelaxationMechanisms;
            using Modelparameter<ValueType>::tauSAverageXY;
            using Modelparameter<ValueType>::tauSAverageXZ;
            using Modelparameter<ValueType>::tauSAverageYZ;
        };
    }
}
