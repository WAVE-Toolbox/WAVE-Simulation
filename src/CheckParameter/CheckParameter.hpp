#pragma once

#include "../Acquisition/Receivers.hpp"
#include "../Acquisition/Sources.hpp"
#include "../Common/Common.hpp"
#include "../Configuration/Configuration.hpp"
#include "../Modelparameter/Modelparameter.hpp"
#include <scai/lama.hpp>

namespace KITGPI
{
    //! \brief CheckParameter namespace
    namespace CheckParameter
    {

        /*! \brief check variable grid
        *
        \param config configuration class
        \param commAll lama communicator
        \param modelCoordinates modelCoordinates object
        */
        template <typename ValueType>
        void checkVariableGrid(const KITGPI::Configuration::Configuration &config, scai::dmemo::CommunicatorPtr commAll, Acquisition::Coordinates<ValueType> const &modelCoordinates)
        {

            IndexType NX = config.get<IndexType>("NX");
            IndexType NY = config.get<IndexType>("NY");
            IndexType NZ = config.get<IndexType>("NZ");
            IndexType newNX = modelCoordinates.getNX();
            IndexType newNY = modelCoordinates.getNY();
            IndexType newNZ = modelCoordinates.getNZ();
            if ((NX != newNX) || (NY != newNY) || (NZ != newNZ)) {
                HOST_PRINT(commAll, "\nIn order to fit the variable grid, the model dimension had been altered from:" << NX << " X " << NY << " X " << NZ << " to " << newNX << " X " << newNY << " X " << newNZ << " \n");
            }
            auto newInterfaces = modelCoordinates.getInterfaceVec();
            std::vector<int> interface;
            std::ifstream is(config.get<std::string>("interfaceFilename"));
            std::istream_iterator<IndexType> start(is), end;
            interface.assign(start, end);

            for (unsigned int i = 0; i < interface.size(); i++) {
                if (interface[i] != newInterfaces[i + 1])
                    HOST_PRINT(commAll, "In order to fit the variable grid, the interface Nr." << i + 1 << " has benn moved from Y=" << interface[i] << " to Y=" << newInterfaces[i + 1] << "\n\n")
            }
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
        void checkStabilityCriterion(ValueType dt, ValueType DH, ValueType vpMax, std::string dimension, scai::IndexType spFDo, scai::dmemo::CommunicatorPtr comm, IndexType shotNumber = -1, IndexType layer = -1)
        {
            if (comm->getRank() == 0) {
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

                if (dt > DH / (h * sqrt(D) * vpMax)) {
                    std::stringstream message;
                    message << "! \ndt is " << dt << " but should be less than DH/(h*sqrt(D)*vpMax"
                            << "=" << DH << "/ (" << h << "* sqrt(" << D << ") * " << vpMax << ") =" << DH / (h * sqrt(D) * vpMax);

                    if ((shotNumber >= 0) && (layer >= 0))
                        HOST_PRINT(comm, "\nCourant-Friedrichs-Lewy-Criterion is not met for shot number: " << shotNumber << " in layer: " << layer << message.str() << "\n\n");
                    if ((shotNumber >= 0) && (layer < 0))
                        HOST_PRINT(comm, "\nCourant-Friedrichs-Lewy-Criterion is not met for shot number: " << shotNumber << message.str() << "\n\n");
                    if ((shotNumber < 0) && (layer < 0))
                        HOST_PRINT(comm, "\nCourant-Friedrichs-Lewy-Criterion is not met" << message.str() << "\n\n");
                }

                SCAI_ASSERT_ERROR(dt <= DH / (h * sqrt(D) * vpMax), "\n\nCourant-Friedrichs-Lewy-Criterion is not met! \n\n");
            }
        }

        /*! \brief check criterion to avoid numerical dispersion
        \param DH Spatial sampling interval in meters.
        \param vMin Minimum wave velocity. (in the acoustic case the minimal P wave velocity is used; in the elastic and viscoelastic cases the minimal S wave velocity is used)
        \param fcMax Maximum center frequency of the sources in hertz.
        \param spFDo Spatial FD order.
        */
        template <typename ValueType>
        void checkNumericalDispersion(ValueType DH, ValueType vMin, ValueType fcMax, scai::IndexType spFDo, scai::dmemo::CommunicatorPtr comm, IndexType shotNumber = -1, IndexType layer = -1)
        {
            if (comm->getRank() == 0) {
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
                    std::stringstream message;
                    message << "! \nDH is " << DH << " but should be less than vMin/(2*fcMax*N)=" << vMin << " / (2 * " << fcMax << " * " << N << ") = " << vMin / (2 * fcMax * N);

                    if ((shotNumber >= 0) && (layer >= 0))
                        HOST_PRINT(comm, "\nCriterion to avoid numerical dispersion is not met for shot number: " << shotNumber << " in layer: " << layer << message.str() << "\n\n");
                    if ((shotNumber >= 0) && (layer < 0))
                        HOST_PRINT(comm, "\nCriterion to avoid numerical dispersion is not met for shot number: " << shotNumber << message.str() << "\n\n");
                    if ((shotNumber < 0) && (layer < 0))
                        HOST_PRINT(comm, "\nCriterion to avoid numerical dispersion is not met" << message.str() << "\n\n");
                }
            }
        }

        //! \brief Wrapper Function who calls checkStabilityCriterion and checkNumericalDispersion
        template <typename ValueType>
        void checkNumericalArtefeactsAndInstabilities(const KITGPI::Configuration::Configuration &config, std::vector<Acquisition::sourceSettings<ValueType>> sourceSettings, Modelparameter::Modelparameter<ValueType> &model, Acquisition::Coordinates<ValueType> const &modelCoordinates, IndexType shotNumber = -1)
        {
            if (!config.get<bool>("initSourcesFromSU")) {
                auto dist = model.getDensity().getDistributionPtr();
                auto commShot = dist->getCommunicatorPtr();
                hmemo::HArray<IndexType> ownedIndexes; // all (global) points owned by this process
                dist->getOwnedIndexes(ownedIndexes);

                IndexType numlayer = modelCoordinates.getNumLayers();

                ValueType vMax[numlayer] = {0};
                ValueType vMin[numlayer] = {0};
                std::fill_n(vMin, numlayer, 3e8);

                IndexType localIndex = 0;

                scai::lama::DenseVector<ValueType> vMaxTmp;
                vMaxTmp = (config.get<std::string>("equationType").compare("sh") == 0) ? model.getVelocityS() : model.getVelocityP();
                auto read_vMaxTmp = hmemo::hostReadAccess(vMaxTmp.getLocalValues());

                scai::lama::DenseVector<ValueType> vMinTmp;
                vMinTmp = (config.get<std::string>("equationType").compare("acoustic") == 0) ? model.getVelocityP() : model.getVelocityS();
                auto read_vMinTmp = hmemo::hostReadAccess(vMinTmp.getLocalValues());

                for (IndexType ownedIndex : hmemo::hostReadAccess(ownedIndexes)) {

                    Acquisition::coordinate3D coordinate = modelCoordinates.index2coordinate(ownedIndex);
                    auto layer = modelCoordinates.getLayer(coordinate);

                    if ((config.get<std::string>("equationType").compare("elastic") == 0) || (config.get<std::string>("equationType").compare("visco") == 0)) {
                        SCAI_ASSERT_ERROR(read_vMaxTmp[localIndex] / read_vMinTmp[localIndex] >= 1, "\n vp/vs (" << read_vMaxTmp[localIndex] << "/" << read_vMinTmp[localIndex] << ") < 1 at X,Y,Z =" << coordinate.x << "," << coordinate.y << "," << coordinate.z << "\n\n");
                    }
                    if (read_vMaxTmp[localIndex] > vMax[layer]) {
                        vMax[layer] = read_vMaxTmp[localIndex];
                    }

                    if ((read_vMinTmp[localIndex] < vMin[layer]) && (read_vMinTmp[localIndex] > 0.0)) {
                        vMin[layer] = read_vMinTmp[localIndex];
                    }

                    localIndex++;
                }

                // communicate vMin and vMax and find the minimu and maximum of all processes
                commShot->minImpl(vMin, vMin, numlayer, common::TypeTraits<ValueType>::stype);
                commShot->maxImpl(vMax, vMax, numlayer, common::TypeTraits<ValueType>::stype);

                ValueType fcMax = 0;
                for (unsigned i = 0; i < sourceSettings.size(); i++) {
                    if (sourceSettings[i].fc > fcMax)
                        fcMax = sourceSettings[i].fc;
                }

                std::vector<IndexType> spatialFDorderVec;

                if (config.get<bool>("useVariableFDoperators")) {

                    std::ifstream is(config.get<std::string>("spatialFDorderFilename"));
                    if (!is)
                        COMMON_THROWEXCEPTION(" could not open " << config.get<std::string>("spatialFDorderFilename"));
                    std::istream_iterator<IndexType> start(is), end;
                    spatialFDorderVec.assign(start, end);
                    if (spatialFDorderVec.empty()) {
                        COMMON_THROWEXCEPTION("FDorder file is empty");
                    }
                } else {
                    spatialFDorderVec.assign(numlayer,config.get<scai::IndexType>("spatialFDorder"));
                }

                for (IndexType layer = 0; layer < numlayer; layer++) {
                    if (config.get<bool>("useVariableGrid")) {
                        checkStabilityCriterion<ValueType>(config.get<ValueType>("DT"), modelCoordinates.getDH(layer), vMax[layer], config.get<std::string>("dimension"), spatialFDorderVec[layer], commShot, shotNumber, layer);
                        checkNumericalDispersion<ValueType>(modelCoordinates.getDH(layer), vMin[layer], fcMax, spatialFDorderVec[layer], commShot, shotNumber, layer);
                    } else {
                        checkStabilityCriterion<ValueType>(config.get<ValueType>("DT"), modelCoordinates.getDH(layer), vMax[layer], config.get<std::string>("dimension"), spatialFDorderVec[layer], commShot, shotNumber);
                        checkNumericalDispersion<ValueType>(modelCoordinates.getDH(layer), vMin[layer], fcMax, spatialFDorderVec[layer], commShot, shotNumber);
                    }
                }
            }
        }

        /*! \brief check if sources are located within the grid
        \param NX Number of gridpoints in x-direction.
        \param NY Number of gridpoints in y-direction.
        \param NZ Number of gridpoints in z-direction.
        \param sourcefile Name of the source file 
        */

        template <typename ValueType>
        void checkSources(Configuration::Configuration const &config, std::vector<Acquisition::sourceSettings<ValueType>> sourceSettings_temp, scai::dmemo::CommunicatorPtr comm)
        {
            scai::IndexType NX = config.get<scai::IndexType>("NX");
            scai::IndexType NY = config.get<scai::IndexType>("NY");
            scai::IndexType NZ = config.get<scai::IndexType>("NZ");

            scai::IndexType numRows = sourceSettings_temp.size();
            scai::IndexType X_temp;
            scai::IndexType Y_temp;
            scai::IndexType Z_temp;

            for (scai::IndexType row_ind = 0; row_ind < numRows; row_ind++) {
                X_temp = sourceSettings_temp[row_ind].sourceCoords.x;
                Y_temp = sourceSettings_temp[row_ind].sourceCoords.y;
                Z_temp = sourceSettings_temp[row_ind].sourceCoords.z;

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
        void checkReceivers(Configuration::Configuration const &config, Acquisition::Receivers<ValueType> &receiver, scai::dmemo::CommunicatorPtr comm)
        {
            scai::IndexType NX = config.get<scai::IndexType>("NX");
            scai::IndexType NY = config.get<scai::IndexType>("NY");
            scai::IndexType NZ = config.get<scai::IndexType>("NZ");

            //             scai::lama::DenseMatrix<ValueType> acquisition_temp;
            std::vector<Acquisition::receiverSettings> acquisition_temp;
            receiver.getAcquisitionMat(config, acquisition_temp);
            scai::IndexType numRows = acquisition_temp.size();
            scai::IndexType X_temp;
            scai::IndexType Y_temp;
            scai::IndexType Z_temp;

            for (scai::IndexType row_ind = 0; row_ind < numRows; row_ind++) {
                X_temp = acquisition_temp[row_ind].receiverCoords.x;
                Y_temp = acquisition_temp[row_ind].receiverCoords.y;
                Z_temp = acquisition_temp[row_ind].receiverCoords.z;
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
