#pragma once

#include "../Acquisition/Receivers.hpp"
#include "../Acquisition/Sources.hpp"
#include "../Acquisition/suHandler.hpp"
#include "../Common/Common.hpp"
#include "../Configuration/Configuration.hpp"
#include "../Modelparameter/Modelparameter.hpp"
#include <scai/lama.hpp>

namespace KITGPI
{
    //! \brief CheckParameter namespace
    namespace CheckParameter
    {

        /*! \brief check number of processes
        *
        \param config configuration class
        \param commAll lama communicator
        */
        void checkNumberOfProcesses(const KITGPI::Configuration::Configuration &config, scai::dmemo::CommunicatorPtr commAll)
        {
            IndexType npS = config.get<IndexType>("ProcNS");
            IndexType npM = commAll->getSize() / npS;

            SCAI_ASSERT_ERROR(commAll->getSize() == npS * npM, "\n Error: Number of MPI processes (" << commAll->getSize() << ") is not multiple of shot domains in configuration"
                                                                                                                 << ": ProcNS = " << npS << "\n")
        }

        /*! \brief check Courant-Friedrichs-Lewy-Criterion
        *
        \param dt Temporal sampling interval in seconds. 
        \param DH Spatial sampling interval in meters.
        \param vpMax Maximum P wave velocity.
        \param dimension Dimension specified in config.txt 
        \param spFDo Spatial FD order
        */
        template <typename ValueType>
        void checkStabilityCriterion(ValueType dt, ValueType DH, ValueType vpMax, std::string dimension, scai::IndexType spFDo, scai::dmemo::CommunicatorPtr comm)
        {
            if (comm->getRank() == MASTERGPI) {
                scai::IndexType D;
                if (dimension.compare("2D") == 0) {
                    D = 2;
                } else if (dimension.compare("3D") == 0) {
                    D = 3;
                } else {
                    SCAI_ASSERT_ERROR(false, "Unknown dimension");
                }

                ValueType h;
                switch (spFDo) {
                case 2:
                    h = 1;
                    break;
                case 4:
                    h = 7 / 6;
                    break;
                case 6:
                    h = 149 / 120;
                    break;
                case 8:
                    h = 2161 / 1680;
                    break;
                case 10:
                    h = 53089 / 40320;
                    break;
                case 12:
                    h = 1187803 / 887040;
                    break;
                default:
                    SCAI_ASSERT_ERROR(false, "Unknown spatial FD order")
                }

                //Assess stability criterion
                SCAI_ASSERT_ERROR(dt <= DH / (h * sqrt(D) * vpMax), "\nCourant-Friedrichs-Lewy-Criterion is not met! \ndt is " << dt << " but should be less than DH/(h*sqrt(D)*vpMax=" << DH / (h * sqrt(D) * vpMax) << "\n\n");
            }
        }

        /*! \brief check criterion to avoid numerical dispersion
        \param DH Spatial sampling interval in meters.
        \param vMin Minimum wave velocity. (in the acoustic case the minimal P wave velocity is used; in the elastic and viscoelastic cases the minimal S wave velocity is used)
        \param fcMax Maximum center frequency of the sources in hertz.
        \param spFDo Spatial FD order.
        */
        template <typename ValueType>
        void checkNumericalDispersion(ValueType DH, ValueType vMin, ValueType fcMax, scai::IndexType spFDo, scai::dmemo::CommunicatorPtr comm)
        {
            if (comm->getRank() == MASTERGPI) {
                scai::IndexType N;
                switch (spFDo) {
                case 2:
                    N = 12;
                    break;
                case 4:
                    N = 8;
                    break;
                case 6:
                    N = 7;
                    break;
                case 8:
                    N = 6;
                    break;
                case 10:
                    N = 5;
                    break;
                case 12:
                    N = 4;
                    break;
                default:
                    SCAI_ASSERT_ERROR(false, "Unknown spatial FD order")
                }

                if (DH > vMin / (2 * fcMax * N)) {
                    std::cerr << "\nCriterion to avoid numerical dispersion is not met! \nDH is " << DH << " but should be less than vMin/(2*fcMax*N)=" << vMin / (2 * fcMax * N) << "\n\n";
                }
            }
        }

        //! \brief Wrapper Function who calls checkStabilityCriterion and checkNumericalDispersion
        template <typename ValueType>
        void checkNumericalArtefeactsAndInstabilities(const KITGPI::Configuration::Configuration &config, Modelparameter::Modelparameter<ValueType> &model, scai::dmemo::CommunicatorPtr comm)
        {
            if (!config.get<bool>("initSourcesFromSU")) {
                ValueType vMaxTmp;
                vMaxTmp = (config.get<std::string>("equationType").compare("sh") == 0) ? model.getVelocityS().max() : model.getVelocityP().max();

                ValueType vMinTmp;
                scai::lama::DenseVector<ValueType> velocityTmp;
                velocityTmp = (config.get<std::string>("equationType").compare("acoustic") == 0) ? model.getVelocityP() : model.getVelocityS();
                KITGPI::Common::searchAndReplace<ValueType>(velocityTmp, 0.0, vMaxTmp, 5);
                vMinTmp = velocityTmp.min();

                scai::lama::DenseMatrix<ValueType> acquisition_temp;
                scai::lama::DenseVector<ValueType> wavelet_fc;
                acquisition_temp.readFromFile(config.get<std::string>("SourceFilename") + ".mtx");
                acquisition_temp.getColumn(wavelet_fc, 8);
                ValueType fcMax = wavelet_fc.max();

                checkStabilityCriterion<ValueType>(config.get<ValueType>("DT"), config.get<ValueType>("DH"), vMaxTmp, config.get<std::string>("dimension"), config.get<scai::IndexType>("spatialFDorder"), comm);
                checkNumericalDispersion<ValueType>(config.get<ValueType>("DH"), vMinTmp, fcMax, config.get<scai::IndexType>("spatialFDorder"), comm);
            }
        }

        /*! \brief check if sources are located within the grid
        \param NX Number of gridpoints in x-direction.
        \param NY Number of gridpoints in y-direction.
        \param NZ Number of gridpoints in z-direction.
        \param sourcefile Name of the source file 
        */

        template <typename ValueType>
        void checkSources(Configuration::Configuration const &config, Acquisition::Sources<ValueType> const &sources, scai::dmemo::CommunicatorPtr comm)
        {
            scai::IndexType NX = config.get<scai::IndexType>("NX");
            scai::IndexType NY = config.get<scai::IndexType>("NY");
            scai::IndexType NZ = config.get<scai::IndexType>("NZ");

            scai::lama::DenseMatrix<ValueType> acquisition_temp;
            sources.getAcquisitionMat(config, acquisition_temp);
            scai::IndexType numRows = acquisition_temp.getNumRows();
            scai::lama::DenseVector<ValueType> row_temp;
            scai::IndexType X_temp;
            scai::IndexType Y_temp;
            scai::IndexType Z_temp;

            for (scai::IndexType row_ind = 0; row_ind < numRows; row_ind++) {
                acquisition_temp.getRow(row_temp, row_ind);
                X_temp = row_temp.getValue(0);
                Y_temp = row_temp.getValue(1);
                Z_temp = row_temp.getValue(2);
                if (comm->getRank() == MASTERGPI) {
                    SCAI_ASSERT_ERROR(X_temp >= 0, "X coordinate (defined as gridpoint) of source #" << row_ind + 1 << " must be positive!")
                    SCAI_ASSERT_ERROR(Y_temp >= 0, "Y coordinate (defined as gridpoint) of source #" << row_ind + 1 << " must be positive!")
                    SCAI_ASSERT_ERROR(Z_temp >= 0, "Z coordinate (defined as gridpoint) of source #" << row_ind + 1 << " must be positive!")
                    SCAI_ASSERT_ERROR(X_temp < NX, "X coordinate (defined as gridpoint) of source #" << row_ind + 1 << " must be smaller than NX!")
                    SCAI_ASSERT_ERROR(Y_temp < NY, "Y coordinate (defined as gridpoint) of source #" << row_ind + 1 << " must be smaller than NY!")
                    SCAI_ASSERT_ERROR(Z_temp < NZ, "Z coordinate (defined as gridpoint) of source #" << row_ind + 1 << " must be smaller than NZ!")
                }
            }
        }

        /*! \brief check if receivers are located within the grid
        \param NX Number of gridpoints in x-direction.
        \param NY Number of gridpoints in y-direction.
        \param NZ Number of gridpoints in z-direction.
        \param receiverfile Name of the receiver file 
        */

        template <typename ValueType>
        void checkReceivers(Configuration::Configuration const &config, Acquisition::Receivers<ValueType> const &receiver, scai::dmemo::CommunicatorPtr comm)
        {
            scai::IndexType NX = config.get<scai::IndexType>("NX");
            scai::IndexType NY = config.get<scai::IndexType>("NY");
            scai::IndexType NZ = config.get<scai::IndexType>("NZ");

            scai::lama::DenseMatrix<ValueType> acquisition_temp;
            receiver.getAcquisitionMat(config, acquisition_temp);
            scai::IndexType numRows = acquisition_temp.getNumRows();
            scai::lama::DenseVector<ValueType> row_temp;
            scai::IndexType X_temp;
            scai::IndexType Y_temp;
            scai::IndexType Z_temp;

            for (scai::IndexType row_ind = 0; row_ind < numRows; row_ind++) {
                acquisition_temp.getRow(row_temp, row_ind);
                X_temp = row_temp.getValue(0);
                Y_temp = row_temp.getValue(1);
                Z_temp = row_temp.getValue(2);
                if (comm->getRank() == MASTERGPI) {
                    SCAI_ASSERT_ERROR(X_temp >= 0, "X coordinate (defined as gridpoint) of receiver #" << row_ind + 1 << " must be positive!")
                    SCAI_ASSERT_ERROR(Y_temp >= 0, "Y coordinate (defined as gridpoint) of receiver #" << row_ind + 1 << " must be positive!")
                    SCAI_ASSERT_ERROR(Z_temp >= 0, "Z coordinate (defined as gridpoint) of receiver #" << row_ind + 1 << " must be positive!")
                    SCAI_ASSERT_ERROR(X_temp < NX, "X coordinate (defined as gridpoint) of receiver #" << row_ind + 1 << " must be smaller than NX!")
                    SCAI_ASSERT_ERROR(Y_temp < NY, "Y coordinate (defined as gridpoint) of receiver #" << row_ind + 1 << " must be smaller than NY!")
                    SCAI_ASSERT_ERROR(Z_temp < NZ, "Z coordinate (defined as gridpoint) of receiver #" << row_ind + 1 << " must be smaller than NZ!")
                }
            }
        }
    }
}
