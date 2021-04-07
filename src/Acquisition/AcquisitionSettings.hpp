#pragma once

#include "Acquisition.hpp"
#include "Coordinates.hpp"
#include <fstream>
#include <scai/dmemo/GeneralDistribution.hpp>
#include <scai/lama.hpp>
#include <algorithm>
#include "../Common/HostPrint.hpp"

namespace KITGPI
{

    namespace Acquisition
    {
        /*! \brief Struct to save settings of a source
         *
         * This struct saves the coordinates and parameters of a source
         */
        template <typename ValueType>
        struct sourceSettings {
            scai::IndexType sourceNo;
            coordinate3D sourceCoords;
            scai::IndexType sourceType;
            scai::IndexType waveletType;
            scai::IndexType waveletShape;
            ValueType fc;
            ValueType amp;
            ValueType tShift;
            coordinate3D getCoords() { return sourceCoords; }
            scai::IndexType getType() { return sourceType; }
            scai::IndexType row;
        };

        /*! \brief Struct to save 3-D coordinates
         *
         * This struct saves the coordinates and parameters of a receiver
         */
        struct receiverSettings {
            coordinate3D receiverCoords;
            scai::IndexType receiverType;
            coordinate3D getCoords() { return receiverCoords; }
            scai::IndexType getType() { return receiverType; }
        };

        /*! \brief Read all source settings into vector of sourceSettings
        *
         \param allSettings vector of sourceSettings structs
        \param fileName Name of source file
        */
        template <typename ValueType>
        inline void readAllSettings(std::vector<sourceSettings<ValueType>> &allSettings, std::string fileName)
        {
            std::ifstream istream(fileName, std::ios::binary);

            allSettings.clear();
            sourceSettings<ValueType> thisSettings;
            scai::IndexType row_index = -1;
            if (istream.is_open()) {
                std::string line;
                while (getline(istream, line)) {
                    std::stringstream strings(line);
                    std::vector<std::string> vecStrings;

                    char firstChar = strings.peek();
                    if ((firstChar == '#') || (line.empty() || (std::all_of(line.begin(), line.end(), isspace)))) {
                        continue;
                    } else {
                        row_index++;
                        std::string tempStr;
                        while (strings >> tempStr) {
                            vecStrings.push_back(tempStr);
                        }

                        if (vecStrings.size() != 10) {
                            COMMON_THROWEXCEPTION("Wrong number of parameters in line of source acquisition file (" << fileName << ")")
                        }

                        try {
                            thisSettings.sourceNo = std::stoi(vecStrings[0]);
                            thisSettings.sourceCoords.x = std::stoi(vecStrings[1]);
                            thisSettings.sourceCoords.y = std::stoi(vecStrings[2]);
                            thisSettings.sourceCoords.z = std::stoi(vecStrings[3]);
                            thisSettings.sourceType = std::stoi(vecStrings[4]);
                            thisSettings.waveletType = std::stoi(vecStrings[5]);
                            thisSettings.waveletShape = std::stoi(vecStrings[6]);
                            thisSettings.fc = std::stof(vecStrings[7]);
                            thisSettings.amp = std::stof(vecStrings[8]);
                            thisSettings.tShift = std::stof(vecStrings[9]);
                        }

                        catch (const std::invalid_argument &ia) {
                            COMMON_THROWEXCEPTION("Invalid argument while reading file " << fileName << " Bad line: " << line << " Message: " << ia.what());
                        }

                        catch (const std::out_of_range &oor) {
                            COMMON_THROWEXCEPTION("Argument out of range while reading file " << fileName << " Bad line: " << line << oor.what());
                        }
                        thisSettings.row = row_index;
                        allSettings.push_back(thisSettings);
                    }
                }
            } else {
                COMMON_THROWEXCEPTION("Could not open source acquisition file " << fileName)
            }
        }

        /*! \brief Cut settings for all source settings into settings for one shot number
        *
        \param settings vector of sourceSettings structs corresponding to shotNumber
        \param allSettings vector of sourceSettings structs with settings for all shots
        \param shotNumber shotNumber to extract corresponding source settings
        
        */
        template <typename ValueType>
        inline void createSettingsForShot(std::vector<KITGPI::Acquisition::sourceSettings<ValueType>> &settings, std::vector<KITGPI::Acquisition::sourceSettings<ValueType>> allSettings, scai::IndexType shotNumber)
        {
            settings.clear();
            for (unsigned i = 0; i < allSettings.size(); i++) {
                if (allSettings[i].sourceNo == shotNumber) {
                    settings.push_back(allSettings[i]);
                }
            }
        }

        /*! \brief adjust sourceSettings for big model to settings for one model pershot
        *
        \param settings vector of sourceSettings structs corresponding to model pershot 
        \param allSettings vector of sourceSettings structs with settings for big model
        \param cutCoordinates coordinates to extract corresponding model pershot        
        */
        template <typename ValueType>
        inline void getSettingsPerShot(std::vector<KITGPI::Acquisition::sourceSettings<ValueType>> &settings, std::vector<KITGPI::Acquisition::sourceSettings<ValueType>> const allSettings, std::vector<KITGPI::Acquisition::coordinate3D> const cutCoordinates)
        {
            settings.clear();
            settings = allSettings;
            for (unsigned i = 0; i < allSettings.size(); i++) {
                settings[i].sourceCoords.x -= cutCoordinates[i].x;
                settings[i].sourceCoords.y -= cutCoordinates[i].y;
                settings[i].sourceCoords.z -= cutCoordinates[i].z;
            }
        }
        
        /*! \brief adjust receiverSettings for big model to settings for one model pershot
        *
        \param settings vector of receiverSettings structs corresponding to model pershot 
        \param allSettings vector of receiverSettings structs with settings for big model
        \param cutCoordinates coordinates to extract corresponding model pershot        
        */
        template <typename ValueType>
        inline void getSettingsPerShot(std::vector<KITGPI::Acquisition::receiverSettings> &settings, std::vector<KITGPI::Acquisition::receiverSettings> const allSettings, std::vector<KITGPI::Acquisition::coordinate3D> const cutCoordinates, scai::IndexType shotIndTrue)
        {
            settings.clear();
            settings = allSettings;
            for (unsigned i = 0; i < allSettings.size(); i++) {
                settings[i].receiverCoords.x -= cutCoordinates[shotIndTrue].x;
                settings[i].receiverCoords.y -= cutCoordinates[shotIndTrue].y;
                settings[i].receiverCoords.z -= cutCoordinates[shotIndTrue].z;
            }
        }

        /*! \brief compute vector of unique shot numbers
        *
        \param uniqueShotNo vector with all shot numbers, each included only once
        \param sourceSettings vector of sourceSettings structs
        */
        template <typename ValueType>
        inline void calcuniqueShotNo(std::vector<scai::IndexType> &uniqueShotNo, std::vector<sourceSettings<ValueType>> sourceSettings)
        {
            uniqueShotNo.clear();
            uniqueShotNo.push_back(sourceSettings[0].sourceNo);
            for (unsigned i = 0; i < sourceSettings.size(); i++) {
                if (std::find(uniqueShotNo.begin(), uniqueShotNo.end(), sourceSettings[i].sourceNo) != uniqueShotNo.end()) {
                    // shotNo already included
                } else {
                    uniqueShotNo.push_back(sourceSettings[i].sourceNo);
                }
            }
        }

        /*! \brief compute vector of unique shot numbers
        *
        \param uniqueShotNo vector with all shot numbers, each included only once
        \param sourceSettings vector of sourceSettings structs
        */
        inline void getuniqueShotInd(scai::IndexType &shotInd, std::vector<scai::IndexType> uniqueShotNo, scai::IndexType shotNumber)
        {
            for (unsigned i = 0; i < uniqueShotNo.size(); i++) {
                if (uniqueShotNo[i] == shotNumber) {
                    shotInd = i;
                    break;
                }
            }
        }
        
        /*! \brief Read all receiver settings into vector of receiverSettings
         *
         \param allSettings vector of receiverSettings structs
         \param fileName Name of receiver file
         */
        inline void readAllSettings(std::vector<receiverSettings> &allSettings, std::string fileName)
        {
            std::ifstream istream(fileName, std::ios::binary);

            allSettings.clear();
            receiverSettings thisSettings;
            if (istream.is_open()) {
                std::string line;
                while (getline(istream, line)) {
                    std::stringstream strings(line);
                    std::vector<std::string> vecStrings;

                    char firstChar = strings.peek();
                    if ((firstChar == '#') || (line.empty() || (std::all_of(line.begin(), line.end(), isspace)))) {
                        continue;
                    } else {
                        std::string tempStr;
                        while (strings >> tempStr) {
                            vecStrings.push_back(tempStr);
                        }

                        if (vecStrings.size() != 4) {
                            COMMON_THROWEXCEPTION("Wrong number of parameters in line of receiver acquisition file (" << fileName << ")")
                        }

                        try {
                            thisSettings.receiverCoords.x = std::stoi(vecStrings[0]);
                            thisSettings.receiverCoords.y = std::stoi(vecStrings[1]);
                            thisSettings.receiverCoords.z = std::stoi(vecStrings[2]);
                            thisSettings.receiverType = std::stoi(vecStrings[3]);
                        }

                        catch (const std::invalid_argument &ia) {
                            COMMON_THROWEXCEPTION("Invalid argument while reading file " << fileName << " Bad line: " << line << " Message: " << ia.what());
                        }

                        catch (const std::out_of_range &oor) {
                            COMMON_THROWEXCEPTION("Argument out of range while reading file " << fileName << " Bad line: " << line << oor.what());
                        }
                        allSettings.push_back(thisSettings);
                    }
                }
            } else {
                COMMON_THROWEXCEPTION("Could not open receiver acquisition file " << fileName)
            }
        }

        /*! \brief Determination of local indices based on given global indeces
        *
        * Calculate the number of indeces within the local processing unit as well as
        * the indeces of the local index.
        *
        \param coordinatesglobal DenseVector with global coordinates
        \param localIndices DenseVector with local coordinates
        \param dist Distribution of global grid
        */
        inline void Global2Local(scai::lama::Vector<scai::IndexType> const &coordinatesglobal, scai::hmemo::HArray<scai::IndexType> &localIndices, scai::dmemo::DistributionPtr dist)
        {

            scai::IndexType n_global = coordinatesglobal.size(); // Number of global entries

            scai::IndexType coordinatetemp_int;

            scai::IndexType i = 0;
            for (scai::IndexType n = 0; n < n_global; n++) {

                coordinatetemp_int = coordinatesglobal.getValue(n);

                SCAI_ASSERT(coordinatetemp_int >= 0 && coordinatetemp_int < dist->getGlobalSize(), "Message from AcquisitionSettings.hpp function Global2Local : Index " << coordinatetemp_int << " is not inside the model grid");
                if (dist->isLocal(coordinatetemp_int)) {
                    i++;
                }
            }

            /* Determine coordinates of local receivers in the global coordinate vector */
            localIndices.resize(i);
            scai::hmemo::WriteAccess<scai::IndexType> write_localIndices(localIndices);
            i = 0;
            for (scai::IndexType n = 0; n < n_global; n++) {

                coordinatetemp_int = coordinatesglobal.getValue(n);
                if (dist->isLocal(coordinatetemp_int)) {
                    write_localIndices[i] = n;
                    i++;
                }
            }
        }

        /*! \brief Getter method for distribution of local traces 
        *
        \param coordinates coordinates
        \param dist_wavefield Distribution of the wavefields
        */
        scai::dmemo::DistributionPtr inline calcDistribution(scai::lama::DenseVector<scai::IndexType> const &coordinates, scai::dmemo::DistributionPtr const dist_wavefield)
        {
            SCAI_ASSERT_DEBUG(coordinates.size() > 0, "The vector coordinates does not contain any elements! ");

            scai::hmemo::HArray<scai::IndexType> localIndices;

            Global2Local(coordinates, localIndices, dist_wavefield);

            scai::dmemo::DistributionPtr dist_temp(new scai::dmemo::GeneralDistribution(coordinates.size(), localIndices, true, dist_wavefield->getCommunicatorPtr()));

            return (dist_temp);
        }

        /*! \brief Getter method for source coordinates from sourceSettings 
        *
        \param sourceSettings sourceSettings
        */
        template <typename ValueType>
        scai::lama::DenseVector<scai::IndexType> getsourcecoordinates(std::vector<sourceSettings<ValueType>> &sourceSettings, Coordinates<ValueType> const &modelCoordinates)
        {
            scai::lama::DenseVector<scai::IndexType> sourcecoords1D;
            sourcecoords1D.allocate(sourceSettings.size());
            for (unsigned i = 0; i < sourceSettings.size(); i++) {
                sourcecoords1D.setValue(i, modelCoordinates.Coordinates<ValueType>::coordinate2index(sourceSettings[i].getCoords()));
            }
            return (sourcecoords1D);
        }
        
        /*! \brief coordinates of cutting model per shot
        \param cutCoordinates coordinates of cutting model per shot
        \param sourceSettingsBig sourceSettings for the big model
        */
        template <typename ValueType>
        inline void getCutCoord(std::vector<KITGPI::Acquisition::coordinate3D> &cutCoordinates, std::vector<Acquisition::sourceSettings<ValueType>> sourceSettingsBig)
        {
            cutCoordinates.clear();
            std::vector<scai::IndexType> uniqueShotNos;
            Acquisition::calcuniqueShotNo(uniqueShotNos, sourceSettingsBig);
            Acquisition::coordinate3D coordinate;
            SCAI_ASSERT(sourceSettingsBig.size() == uniqueShotNos.size(), "sourceSettingsBig.size() != uniqueShotNos.size()");
            std::vector<scai::IndexType> uniqueShotNosX;
            for (unsigned i = 0; i < sourceSettingsBig.size(); i++) {               
                uniqueShotNosX.push_back(sourceSettingsBig[i].sourceCoords.x);
            }
            auto minX = min(uniqueShotNosX.begin(), uniqueShotNosX.end());
            for (unsigned i = 0; i < sourceSettingsBig.size(); i++) {               
                coordinate.x = sourceSettingsBig[i].sourceCoords.x - *minX;
                coordinate.y = 0;
                coordinate.z = 0;
                cutCoordinates.push_back(coordinate);
            }
        }
        
        /*! \brief Write to cutCoordinates-file
        *
        \param cutCoordinatesFilename Name of cutCoordinates-file
        \param cutCoordinates coordinates of cutting model per shot
        \param uniqueShotNos unique shot numbers
        */
        inline void writeCutCoordToFile(std::string cutCoordinatesFilename, std::vector<KITGPI::Acquisition::coordinate3D> cutCoordinates, std::vector<scai::IndexType> uniqueShotNos)
        {
            std::ofstream outputFile(cutCoordinatesFilename);
            outputFile << "# Coordinate for cutting model per shot \n";
            outputFile << "# ShotNumber | index_x | index_y | index_z\n";                
            for (unsigned i = 0; i < uniqueShotNos.size(); i++) {  
                outputFile << std::setw(12) << uniqueShotNos[i] << std::setw(10) << cutCoordinates[i].x << std::setw(10) << cutCoordinates[i].y << std::setw(10) << cutCoordinates[i].z << "\n";
            }
            outputFile.close();
        }
        
        /*! \brief Generate a random shot sequence without any repeated elements and with the max appearances of one shot
        *
        \param uniqueShotNosRand random Shot numbers
        \param filterHistoryCount a vector to count how much times the shot appears
        \param uniqueShotNos unique Shot numbers
        \param maxcount max times for one shot
        \param useRandSource useRandSource
        */
        template <typename ValueType>
        void getRandShotNos(std::vector<scai::IndexType> &uniqueShotNosRand, std::vector<scai::IndexType> &filterHistoryCount, std::vector<scai::IndexType> uniqueShotNos, scai::IndexType maxcount, scai::IndexType useRandSource)
        {  
            scai::IndexType numShotDomains = uniqueShotNosRand.size();
            scai::IndexType numshots = uniqueShotNos.size(); 
            scai::IndexType randShotInd;
              
            if (useRandSource == 1) {
                std::vector<scai::IndexType> randShotIndHistory(numShotDomains, 0);
                std::srand((int)time(0));
                for (scai::IndexType shotDomainInd = 0; shotDomainInd < numShotDomains; shotDomainInd++) {                                     
                    bool repeat = false;
                    randShotInd = std::rand() % numshots;
                    randShotIndHistory[shotDomainInd] = randShotInd;
                    for (scai::IndexType i = 0; i < shotDomainInd; i++) {
                        if (randShotIndHistory[i] == randShotInd) {
                            shotDomainInd--;
                            repeat = true;
                            break;
                        }
                    }
                    if (filterHistoryCount[randShotInd] >= maxcount) {
                        if (!repeat) 
                            shotDomainInd--;
                    } else if (!repeat) {
                        uniqueShotNosRand[shotDomainInd] = uniqueShotNos[randShotInd];
                        filterHistoryCount[randShotInd]++;
                    }
                }
            } else if (useRandSource == 2) {
                scai::IndexType startShotInd = 0;
                if (filterHistoryCount[0] != filterHistoryCount[numshots-1]) {
                    for (scai::IndexType shotInd = 0; shotInd < numshots-1; shotInd++) { 
                        if (filterHistoryCount[shotInd] != filterHistoryCount[shotInd+1]) {
                            startShotInd = shotInd + 1;
                        }
                    }
                }
                for (scai::IndexType shotDomainInd = 0; shotDomainInd < numShotDomains; shotDomainInd++) {                                     
                    randShotInd = startShotInd + shotDomainInd;
                    uniqueShotNosRand[shotDomainInd] = uniqueShotNos[randShotInd];
                    filterHistoryCount[randShotInd]++;
                }
            }
        }
                
        /*! \brief Write to randSource-file
        *
        \param comm Communicator
        \param randSourceFilename Name of randSource-file
        \param uniqueShotNosRand unique Shot numbers randomly
        \param stage inversion stage
        \param iteration inversion iteration
        \param useRandSource useRandSource
        */
        inline void writeRandShotNosToFile(scai::dmemo::CommunicatorPtr comm, std::string randSourceFilename, std::vector<scai::IndexType> uniqueShotNosRand, scai::IndexType stage, scai::IndexType iteration, scai::IndexType useRandSource)
        {      
            int myRank = comm->getRank();  
            if (useRandSource != 0 && myRank == MASTERGPI) {
                std::ofstream outputFile; 
                if (stage == 1 && iteration == 1) {
                    outputFile.open(randSourceFilename);
                    outputFile << "# ShotNumber records during inversion\n"; 
                    outputFile << "# random source type = " << useRandSource << " (0=all sequential shot, 1=numShotDomains random shot, 2=numShotDomains sequential shot)\n"; 
                    outputFile << "# Stage | Iteration | shotNumbers\n"; 
                } else {                    
                    outputFile.open(randSourceFilename, std::ios_base::app);
                    outputFile << std::scientific;
                }
                outputFile << std::setw(5) << stage << std::setw(10) << iteration;
                for (unsigned i = 0; i < uniqueShotNosRand.size(); i++) {  
                    if (i == 0) {
                        outputFile << std::setw(9) << uniqueShotNosRand[i];
                    } else {
                        outputFile << std::setw(4) << uniqueShotNosRand[i];
                    }
                }
                outputFile << "\n";
                outputFile.close();
            }
        }
    }
}
