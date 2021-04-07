
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

#include "Modelparameter.hpp"

namespace KITGPI
{

    //! \brief Modelparameter namespace
    namespace Modelparameter
    {

        //! Class for Modelparameter for acoustic simulations (Subsurface properties)
        /*!
         This class handels the modelparameter for the acoustic finite-difference simulation.
         */
        template <typename ValueType>
        class Acoustic : public Modelparameter<ValueType>
        {
          public:
            //! Default constructor.
            Acoustic() { equationType = "acoustic"; };

            //! Destructor, releases all allocated resources.
            ~Acoustic(){};

            explicit Acoustic(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> const &modelCoordinates);
            explicit Acoustic(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, ValueType pWaveModulus_const, ValueType rho_const);
            explicit Acoustic(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, std::string filename, scai::IndexType fileFormat);

            //! Copy Constructor.
            Acoustic(const Acoustic &rhs);

            ValueType estimateMemory(scai::dmemo::DistributionPtr dist) override;

            void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, ValueType pWaveModulus_const, ValueType rho_const);
            void init(Configuration::Configuration const &config, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, Acquisition::Coordinates<ValueType> const &modelCoordinates) override;
            void init(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, std::string filename, scai::IndexType fileFormat) override;

            void write(std::string filename, scai::IndexType fileFormat) const override;

            std::string getEquationType() const;

            /* Getter methods for not requiered parameters */
            scai::lama::Vector<ValueType> const &getInverseDensity() override;
            scai::lama::Vector<ValueType> const &getSWaveModulus() override;
            scai::lama::Vector<ValueType> const &getSWaveModulus() const override;
            scai::lama::DenseVector<ValueType> const &getVelocityS() const override;
            scai::lama::Vector<ValueType> const &getTauP() const override;
            scai::lama::Vector<ValueType> const &getTauS() const override;
            scai::lama::Vector<ValueType> const &getSWaveModulusAverageXY() override;
            scai::lama::Vector<ValueType> const &getSWaveModulusAverageXY() const override;
            scai::lama::Vector<ValueType> const &getSWaveModulusAverageXZ() override;
            scai::lama::Vector<ValueType> const &getSWaveModulusAverageXZ() const override;
            scai::lama::Vector<ValueType> const &getSWaveModulusAverageYZ() override;
            scai::lama::Vector<ValueType> const &getSWaveModulusAverageYZ() const override;
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
            
            void getModelSubset(KITGPI::Modelparameter::Modelparameter<ValueType> &modelSubset, Acquisition::Coordinates<ValueType> const &modelCoordinates, Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, std::vector<Acquisition::coordinate3D> cutCoordinates, scai::IndexType cutCoordInd) override;
            void setModelSubset(KITGPI::Modelparameter::Modelparameter<ValueType> &invertedModelSubset, Acquisition::Coordinates<ValueType> const &modelCoordinates, Acquisition::Coordinates<ValueType> const &modelCoordinatesBig, std::vector<Acquisition::coordinate3D> cutCoordinates, scai::IndexType cutCoordInd, scai::IndexType smoothRange, scai::IndexType NX, scai::IndexType NY, scai::IndexType NXBig, scai::IndexType NYBig, scai::IndexType boundaryWidth) override;
            
            void minusAssign(KITGPI::Modelparameter::Modelparameter<ValueType> const &rhs) override;
            void plusAssign(KITGPI::Modelparameter::Modelparameter<ValueType> const &rhs) override;
            void assign(KITGPI::Modelparameter::Modelparameter<ValueType> const &rhs) override;

            /* Overloading Operators */
            KITGPI::Modelparameter::Acoustic<ValueType> operator*(ValueType rhs);
            KITGPI::Modelparameter::Acoustic<ValueType> &operator*=(ValueType const &rhs);
            KITGPI::Modelparameter::Acoustic<ValueType> operator+(KITGPI::Modelparameter::Acoustic<ValueType> const &rhs);
            KITGPI::Modelparameter::Acoustic<ValueType> &operator+=(KITGPI::Modelparameter::Acoustic<ValueType> const &rhs);
            KITGPI::Modelparameter::Acoustic<ValueType> operator-(KITGPI::Modelparameter::Acoustic<ValueType> const &rhs);
            KITGPI::Modelparameter::Acoustic<ValueType> &operator-=(KITGPI::Modelparameter::Acoustic<ValueType> const &rhs);
            KITGPI::Modelparameter::Acoustic<ValueType> &operator=(KITGPI::Modelparameter::Acoustic<ValueType> const &rhs);

          private:
            void init(scai::dmemo::DistributionPtr variableDist, Acquisition::Coordinates<ValueType> const &variableCoordinates, Acquisition::Coordinates<ValueType> const &regularCoordinates);
            void calculateAveraging() override;

            using Modelparameter<ValueType>::equationType;

            using Modelparameter<ValueType>::dirtyFlagInverseDensity;
            using Modelparameter<ValueType>::dirtyFlagPWaveModulus;
            using Modelparameter<ValueType>::dirtyFlagAveraging;
            using Modelparameter<ValueType>::pWaveModulus;
            using Modelparameter<ValueType>::density;
            using Modelparameter<ValueType>::inverseDensity;
            using Modelparameter<ValueType>::velocityP;

            void initializeMatrices(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::CommunicatorPtr comm) override;
            void purgeMatrices() override;

            using Modelparameter<ValueType>::averageMatrixX;
            using Modelparameter<ValueType>::averageMatrixY;
            using Modelparameter<ValueType>::averageMatrixZ;

            using Modelparameter<ValueType>::inverseDensityAverageX;
            using Modelparameter<ValueType>::inverseDensityAverageY;
            using Modelparameter<ValueType>::inverseDensityAverageZ;

            /* Not requiered parameters */
            using Modelparameter<ValueType>::velocityS;
            using Modelparameter<ValueType>::sWaveModulus;
            using Modelparameter<ValueType>::tauP;
            using Modelparameter<ValueType>::tauS;
            using Modelparameter<ValueType>::relaxationFrequency;
            using Modelparameter<ValueType>::numRelaxationMechanisms;
            using Modelparameter<ValueType>::sWaveModulusAverageXY;
            using Modelparameter<ValueType>::sWaveModulusAverageXZ;
            using Modelparameter<ValueType>::sWaveModulusAverageYZ;
            using Modelparameter<ValueType>::tauSAverageXY;
            using Modelparameter<ValueType>::tauSAverageXZ;
            using Modelparameter<ValueType>::tauSAverageYZ;
        };
    }
}
