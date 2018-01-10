
#pragma once

#include <scai/lama.hpp>

#include <scai/lama/DenseVector.hpp>
#include <scai/lama/expression/all.hpp>
#include <scai/lama/matrix/all.hpp>
#include <scai/lama/matutils/MatrixCreator.hpp>
#include <scai/lama/norm/L2Norm.hpp>

#include <scai/dmemo/BlockDistribution.hpp>

#include <scai/hmemo/HArray.hpp>
#include <scai/hmemo/ReadAccess.hpp>
#include <scai/hmemo/WriteAccess.hpp>

#include <scai/tracing.hpp>

#include <scai/common/Walltime.hpp>
#include <scai/common/unique_ptr.hpp>
#include <scai/logging.hpp>

#include <iostream>

#include "../PartitionedInOut/PartitionedInOut.hpp"
#include "Modelparameter.hpp"

namespace KITGPI
{

    //! \brief Modelparameter namespace
    namespace Modelparameter
    {

        //! Class for Modelparameter for elastic simulations (Subsurface properties)
        /*!
         This class handels the modelparameter for the elastic finite-difference simulation.
         */
        template <typename ValueType>
        class Elastic : public Modelparameter<ValueType>
        {
          public:
            //! Default constructor.
            Elastic(){};

            //! Destructor, releases all allocated resources.
            ~Elastic(){};

            explicit Elastic(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist);
            explicit Elastic(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, scai::lama::Scalar pWaveModulus_const, scai::lama::Scalar sWaveModulus_const, scai::lama::Scalar rho);
            explicit Elastic(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, std::string filename, IndexType partitionedIn);

            //! Copy Constructor.
            Elastic(const Elastic &rhs);

            void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, scai::lama::Scalar pWaveModulus, scai::lama::Scalar sWaveModulus, scai::lama::Scalar rho);
            void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, std::string filename, IndexType partitionedIn) override;
            void init(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist) override;


            void write(std::string filename, IndexType partitionedOut) const override;

            /* Getter methods for not requiered parameters */
            scai::lama::Vector const &getTauP() const override;
            scai::lama::Vector const &getTauS() const override;
            scai::lama::Vector const &getTauSAverageXY() override;
	    scai::lama::Vector const &getTauSAverageXY() const override;
            scai::lama::Vector const &getTauSAverageXZ() override;
	    scai::lama::Vector const &getTauSAverageXZ() const override;
            scai::lama::Vector const &getTauSAverageYZ() override;
	    scai::lama::Vector const &getTauSAverageYZ() const override;
            IndexType getNumRelaxationMechanisms() const override;
            ValueType getRelaxationFrequency() const override;

            void prepareForModelling(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, scai::dmemo::CommunicatorPtr comm) override;

            /* Overloading Operators */
            KITGPI::Modelparameter::Elastic<ValueType> operator*(scai::lama::Scalar rhs);
            KITGPI::Modelparameter::Elastic<ValueType> &operator*=(scai::lama::Scalar const &rhs);
            KITGPI::Modelparameter::Elastic<ValueType> operator+(KITGPI::Modelparameter::Elastic<ValueType> const &rhs);
            KITGPI::Modelparameter::Elastic<ValueType> &operator+=(KITGPI::Modelparameter::Elastic<ValueType> const &rhs);
            KITGPI::Modelparameter::Elastic<ValueType> operator-(KITGPI::Modelparameter::Elastic<ValueType> const &rhs);
            KITGPI::Modelparameter::Elastic<ValueType> &operator-=(KITGPI::Modelparameter::Elastic<ValueType> const &rhs);
	    KITGPI::Modelparameter::Elastic<ValueType> &operator=(KITGPI::Modelparameter::Elastic<ValueType> const &rhs);

          private:
            void calculateAveraging() override;

            using Modelparameter<ValueType>::dirtyFlagInverseDensity;
            using Modelparameter<ValueType>::dirtyFlagPWaveModulus;
	    using Modelparameter<ValueType>::dirtyFlagSWaveModulus;
            using Modelparameter<ValueType>::dirtyFlagAveraging;
            using Modelparameter<ValueType>::pWaveModulus;
            using Modelparameter<ValueType>::sWaveModulus;
            using Modelparameter<ValueType>::density;
            using Modelparameter<ValueType>::inverseDensity;
            using Modelparameter<ValueType>::velocityP;
            using Modelparameter<ValueType>::velocityS;

            void initializeMatrices(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, IndexType NX, IndexType NY, IndexType NZ, ValueType DH, ValueType DT, scai::dmemo::CommunicatorPtr comm) override;

            void initializeMatrices(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, Configuration::Configuration config, scai::dmemo::CommunicatorPtr comm);

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
