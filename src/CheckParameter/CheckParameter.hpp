#pragma once

#include "../Acquisition/Receivers.hpp"
#include "../Acquisition/Sources.hpp"
#include "../Common/Common.hpp"
#include "../Configuration/Configuration.hpp"
#include "../Modelparameter/Modelparameter.hpp"
#include "../ModelparameterEM/Modelparameter.hpp"
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
        void checkVariableGrid(Configuration::Configuration const &config, scai::dmemo::CommunicatorPtr commAll, Acquisition::Coordinates<ValueType> const &modelCoordinates)
        {
            scai::IndexType NX = config.get<scai::IndexType>("NX");
            scai::IndexType NY = config.get<scai::IndexType>("NY");
            scai::IndexType NZ = config.get<scai::IndexType>("NZ");
            scai::IndexType newNX = modelCoordinates.getNX();
            scai::IndexType newNY = modelCoordinates.getNY();
            scai::IndexType newNZ = modelCoordinates.getNZ();
            if ((NX != newNX) || (NY != newNY) || (NZ != newNZ)) {
                HOST_PRINT(commAll, "\nIn order to fit the variable grid, the model dimension had been altered from:" << NX << " X " << NY << " X " << NZ << " to " << newNX << " X " << newNY << " X " << newNZ << " \n");
            }
            auto newInterfaces = modelCoordinates.getInterfaceVec();
            std::vector<int> interface;
            std::string gridConfigFileName = config.get<std::string>("gridConfigurationFilename");

            unsigned int column = 0;
            Common::readColumnFromFile(gridConfigFileName, interface, column);
            if (interface.at(0) == 0) {
                interface.erase(interface.begin());
            } else {
                COMMON_THROWEXCEPTION("First interface must by at y=0 ");
            }

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
        void checkStabilityCriterion(ValueType dt, ValueType DH, ValueType vpMax, std::string dimension, scai::IndexType spFDo, scai::dmemo::CommunicatorPtr comm, scai::IndexType shotNumber = -1, scai::IndexType layer = -1)
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
        void checkNumericalDispersion(ValueType DH, ValueType vMin, ValueType fcMax, scai::IndexType spFDo, scai::dmemo::CommunicatorPtr comm, scai::IndexType shotNumber = -1, scai::IndexType layer = -1)
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
                        HOST_PRINT(comm, "", "\nCriterion to avoid numerical dispersion is not met for shot number: " << shotNumber << " in layer: " << layer << message.str() << "\n\n");
                    if ((shotNumber >= 0) && (layer < 0))
                        HOST_PRINT(comm, "", "\nCriterion to avoid numerical dispersion is not met for shot number: " << shotNumber << message.str() << "\n\n");
                    if ((shotNumber < 0) && (layer < 0))
                        HOST_PRINT(comm, "", "\nCriterion to avoid numerical dispersion is not met" << message.str() << "\n\n");
                }
            }
        }

        //! \brief Wrapper Function who calls checkStabilityCriterion and checkNumericalDispersion
        template <typename ValueType>
        void checkNumericalArtefactsAndInstabilities(Configuration::Configuration const &config, std::vector<Acquisition::sourceSettings<ValueType>> sourceSettings, Modelparameter::Modelparameter<ValueType> &model, Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::IndexType shotNumber = -1)
        {
            if (!config.get<bool>("initSourcesFromSU")) {
                auto dist = model.getDensity().getDistributionPtr();
                auto commShot = dist->getCommunicatorPtr();
                scai::hmemo::HArray<scai::IndexType> ownedIndexes; // all (global) points owned by this process
                dist->getOwnedIndexes(ownedIndexes);

                scai::IndexType numlayer = modelCoordinates.getNumLayers();

                ValueType vMax[numlayer] = {0};
                ValueType vMin[numlayer] = {0};
                std::fill_n(vMin, numlayer, 3e8);

                scai::IndexType localIndex = 0;

                scai::lama::DenseVector<ValueType> vMaxTmp;
                vMaxTmp = (config.get<std::string>("equationType").compare("sh") == 0) ? model.getVelocityS() : model.getVelocityP();
                auto read_vMaxTmp = scai::hmemo::hostReadAccess(vMaxTmp.getLocalValues());

                scai::lama::DenseVector<ValueType> vMinTmp;
                vMinTmp = (config.get<std::string>("equationType").compare("acoustic") == 0) ? model.getVelocityP() : model.getVelocityS();
                auto read_vMinTmp = scai::hmemo::hostReadAccess(vMinTmp.getLocalValues());

                for (scai::IndexType ownedIndex : scai::hmemo::hostReadAccess(ownedIndexes)) {

                    Acquisition::coordinate3D coordinate = modelCoordinates.index2coordinate(ownedIndex);
                    auto layer = modelCoordinates.getLayer(coordinate);

                    if ((config.get<std::string>("equationType").compare("elastic") == 0) || (config.get<std::string>("equationType").compare("viscoelastic") == 0)) {
                        SCAI_ASSERT_ERROR(read_vMaxTmp[localIndex] / read_vMinTmp[localIndex] >= sqrt(2.0), "\n vp/vs (" << read_vMaxTmp[localIndex] << "/" << read_vMinTmp[localIndex] << ") < sqrt(2.0) at X,Y,Z =" << coordinate.x << "," << coordinate.y << "," << coordinate.z << "\n\n");
                        // vp/vs = sqrt((lambda+2*mu)/mu) >= sqrt(2.0)
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
                commShot->minImpl(vMin, vMin, numlayer, scai::common::TypeTraits<ValueType>::stype);
                commShot->maxImpl(vMax, vMax, numlayer, scai::common::TypeTraits<ValueType>::stype);

                ValueType fcMax = 0;
                for (unsigned i = 0; i < sourceSettings.size(); i++) {
                    if (sourceSettings[i].fc > fcMax)
                        fcMax = sourceSettings[i].fc;
                }

                std::vector<scai::IndexType> spatialFDorderVec;

                if (config.get<bool>("useVariableFDoperators")) {

                    unsigned int column = 2;
                    Common::readColumnFromFile(config.get<std::string>("gridConfigurationFilename"), spatialFDorderVec, column);
                } else {
                    spatialFDorderVec.assign(numlayer, config.get<scai::IndexType>("spatialFDorder"));
                }

                for (scai::IndexType layer = 0; layer < numlayer; layer++) {
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

        //! \brief Wrapper Function who calls checkStabilityCriterion and checkNumericalDispersion
        template <typename ValueType>
        void checkNumericalArtefactsAndInstabilities(Configuration::Configuration const &config, std::vector<Acquisition::sourceSettings<ValueType>> sourceSettings, Modelparameter::ModelparameterEM<ValueType> &model, Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::IndexType shotNumber = -1)
        {
            if (!config.get<bool>("initSourcesFromSU")) {
                auto dist = model.getMagneticPermeabilityEM().getDistributionPtr();
                auto commShot = dist->getCommunicatorPtr();
                scai::hmemo::HArray<scai::IndexType> ownedIndexes; // all (global) points owned by this process
                dist->getOwnedIndexes(ownedIndexes);

                scai::IndexType numlayer = modelCoordinates.getNumLayers();

                ValueType vMax[numlayer] = {0};
                ValueType vMin[numlayer] = {0};
                std::fill_n(vMin, numlayer, 3e8);

                scai::IndexType localIndex = 0;

                scai::lama::DenseVector<ValueType> vMaxTmp;
                vMaxTmp = model.getVelocityEM();
                auto read_vMaxTmp = scai::hmemo::hostReadAccess(vMaxTmp.getLocalValues());

                scai::lama::DenseVector<ValueType> vMinTmp;
                vMinTmp = model.getVelocityEM();
                auto read_vMinTmp = scai::hmemo::hostReadAccess(vMinTmp.getLocalValues());

                for (scai::IndexType ownedIndex : scai::hmemo::hostReadAccess(ownedIndexes)) {

                    Acquisition::coordinate3D coordinate = modelCoordinates.index2coordinate(ownedIndex);
                    auto layer = modelCoordinates.getLayer(coordinate);

                    if (read_vMaxTmp[localIndex] > vMax[layer]) {
                        vMax[layer] = read_vMaxTmp[localIndex];
                    }

                    if ((read_vMinTmp[localIndex] < vMin[layer]) && (read_vMinTmp[localIndex] > 0.0)) {
                        vMin[layer] = read_vMinTmp[localIndex];
                    }

                    localIndex++;
                }

                // communicate vMin and vMax and find the minimu and maximum of all processes
                commShot->minImpl(vMin, vMin, numlayer, scai::common::TypeTraits<ValueType>::stype);
                commShot->maxImpl(vMax, vMax, numlayer, scai::common::TypeTraits<ValueType>::stype);

                ValueType fcMax = 0;
                for (unsigned i = 0; i < sourceSettings.size(); i++) {
                    if (sourceSettings[i].fc > fcMax)
                        fcMax = sourceSettings[i].fc;
                }

                std::vector<scai::IndexType> spatialFDorderVec;

                if (config.get<bool>("useVariableFDoperators")) {

                    unsigned int column = 2;
                    Common::readColumnFromFile(config.get<std::string>("gridConfigurationFilename"), spatialFDorderVec, column);
                } else {
                    spatialFDorderVec.assign(numlayer, config.get<scai::IndexType>("spatialFDorder"));
                }

                for (scai::IndexType layer = 0; layer < numlayer; layer++) {
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
        void checkSources(std::vector<Acquisition::sourceSettings<ValueType>> sourceSettings_temp, Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::CommunicatorPtr comm)
        {
            for (auto const &source : sourceSettings_temp) {

                if ((source.sourceCoords.x >= modelCoordinates.getNX() || source.sourceCoords.x < 0 || source.sourceCoords.y >= modelCoordinates.getNY() || source.sourceCoords.y < 0 || source.sourceCoords.z >= modelCoordinates.getNZ() || source.sourceCoords.z < 0)) {
                    HOST_PRINT(comm, "\nsource Nr. " << source.sourceNo << " with coordinate: " << source.sourceCoords << " is not located on the grid \n\n");
                    COMMON_THROWEXCEPTION("source Nr. " << source.sourceNo << " with coordinate: " << source.sourceCoords << " is not located on the grid");
                }

                int index = modelCoordinates.coordinate2index(source.sourceCoords);
                Acquisition::coordinate3D coordinate = modelCoordinates.index2coordinate(index);

                //check for the variable grid
                if (coordinate != source.sourceCoords) {
                    HOST_PRINT(comm, "\nsource Nr. " << source.sourceNo << " with coordinate: " << source.sourceCoords << " is not located on the grid \n\n");
                    COMMON_THROWEXCEPTION("source Nr. " << source.sourceNo << " with coordinate: " << source.sourceCoords << " is not located on the grid");
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
        void checkReceivers(std::vector<Acquisition::receiverSettings> receiverSettings_temp, Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::CommunicatorPtr comm)
        {
            int i = 0;
            for (auto const &receiver : receiverSettings_temp) {
                i++;

                if ((receiver.receiverCoords.x >= modelCoordinates.getNX() || receiver.receiverCoords.x < 0 || receiver.receiverCoords.y >= modelCoordinates.getNY() || receiver.receiverCoords.y < 0 || receiver.receiverCoords.z >= modelCoordinates.getNZ() || receiver.receiverCoords.z < 0)) {
                    HOST_PRINT(comm, "\nreceiver Nr. " << i << " with coordinate: " << receiver.receiverCoords << " is not located on the model grid \n\n");
                    COMMON_THROWEXCEPTION("\nreceiver Nr. " << i << " with coordinate: " << receiver.receiverCoords << " is not located on the model grid \n\n");
                }

                int index = modelCoordinates.coordinate2index(receiver.receiverCoords);
                Acquisition::coordinate3D coordinate = modelCoordinates.index2coordinate(index);

                //check for the variable grid
                if (coordinate != receiver.receiverCoords) {
                    HOST_PRINT(comm, "\nreceiver Nr. " << i << " with coordinate: " << receiver.receiverCoords << " is not located on the model grid \n\n");
                    COMMON_THROWEXCEPTION("\nreceiver Nr. " << i << " with coordinate: " << receiver.receiverCoords << " is not located on the model grid \n\n");
                }
            }
        }
    }
}
