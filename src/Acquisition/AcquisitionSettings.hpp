#pragma once

#include "Acquisition.hpp"
#include "Coordinates.hpp"
#include <fstream>
#include <scai/dmemo/GeneralDistribution.hpp>
#include <scai/lama.hpp>
#include <algorithm>
#include "../Common/HostPrint.hpp"
#include <numeric>

using namespace scai;

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
            IndexType sourceNo;
            coordinate3D sourceCoords;
            IndexType sourceType;
            IndexType waveletType;
            IndexType waveletShape;
            ValueType fc;
            ValueType amp;
            ValueType tShift;
            coordinate3D getCoords() { return sourceCoords; }
            IndexType getType() { return sourceType; }
            IndexType row;
        };

        /*! \brief Struct to save 3-D coordinates
         *
         * This struct saves the coordinates and parameters of a receiver
         */
        struct receiverSettings {
            coordinate3D receiverCoords;
            IndexType receiverType;
            coordinate3D getCoords() { return receiverCoords; }
            IndexType getType() { return receiverType; }
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
            IndexType row_index = -1;
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
        inline void createSettingsForShot(std::vector<KITGPI::Acquisition::sourceSettings<ValueType>> &settings, std::vector<KITGPI::Acquisition::sourceSettings<ValueType>> allSettings, IndexType shotNumber)
        {
            settings.clear();
            for (unsigned i = 0; i < allSettings.size(); i++) {
                if (std::abs(allSettings[i].sourceNo) == shotNumber) {
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
        inline void getSettingsPerShot(std::vector<KITGPI::Acquisition::receiverSettings> &settings, std::vector<KITGPI::Acquisition::receiverSettings> const allSettings, KITGPI::Acquisition::coordinate3D const cutCoordinate)
        {
            settings.clear();
            settings = allSettings;
            for (unsigned i = 0; i < allSettings.size(); i++) {
                settings[i].receiverCoords.x -= cutCoordinate.x;
                settings[i].receiverCoords.y -= cutCoordinate.y;
                settings[i].receiverCoords.z -= cutCoordinate.z;
            }
        }

        /*! \brief compute vector of unique shot numbers
        *
        \param uniqueShotNos vector with all shot numbers, each included only once
        \param sourceSettings vector of sourceSettings structs
        */
        template <typename ValueType>
        inline void calcuniqueShotNo(std::vector<IndexType> &uniqueShotNos, std::vector<sourceSettings<ValueType>> sourceSettings)
        {
            uniqueShotNos.clear();
            uniqueShotNos.push_back(std::abs(sourceSettings[0].sourceNo));
            for (unsigned i = 0; i < sourceSettings.size(); i++) {
                if (std::find(uniqueShotNos.begin(), uniqueShotNos.end(), std::abs(sourceSettings[i].sourceNo)) != uniqueShotNos.end()) {
                    // shotNo already included
                } else { // only the positive shotnr can be used
                    uniqueShotNos.push_back(std::abs(sourceSettings[i].sourceNo));
                }
            }
        }
        
        /*! \brief compute vector of unique shot numbers
        *
        \param uniqueShotNos vector with all shot numbers, each included only once
        */
        inline void getuniqueShotInd(IndexType &shotInd, std::vector<IndexType> uniqueShotNos, IndexType shotNumber)
        {
            for (unsigned i = 0; i < uniqueShotNos.size(); i++) {
                if (uniqueShotNos[i] == shotNumber) {
                    shotInd = i;
                    break;
                }
            }
        }
        
        /*! \brief compute vector of unique shot numbers
        *
        \param sourceSettingsEncode vector of sourceSettings structs
        */
        template <typename ValueType>
        inline void getuniqueShotInd(IndexType &shotInd, std::vector<sourceSettings<ValueType>> sourceSettingsEncode, IndexType shotNumber)
        {
            for (unsigned i = 0; i < sourceSettingsEncode.size(); i++) {
                if (std::abs(sourceSettingsEncode[i].sourceNo) == shotNumber) {
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
        inline void Global2Local(scai::lama::Vector<IndexType> const &coordinatesglobal, scai::hmemo::HArray<IndexType> &localIndices, scai::dmemo::DistributionPtr dist)
        {

            IndexType n_global = coordinatesglobal.size(); // Number of global entries

            IndexType coordinatetemp_int;

            IndexType i = 0;
            for (IndexType n = 0; n < n_global; n++) {

                coordinatetemp_int = coordinatesglobal.getValue(n);

                SCAI_ASSERT(coordinatetemp_int >= 0 && coordinatetemp_int < dist->getGlobalSize(), "Message from AcquisitionSettings.hpp function Global2Local : Index " << coordinatetemp_int << " is not inside the model grid");
                if (dist->isLocal(coordinatetemp_int)) {
                    i++;
                }
            }

            /* Determine coordinates of local receivers in the global coordinate vector */
            localIndices.resize(i);
            scai::hmemo::WriteAccess<IndexType> write_localIndices(localIndices);
            i = 0;
            for (IndexType n = 0; n < n_global; n++) {

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
        scai::dmemo::DistributionPtr inline calcDistribution(scai::lama::DenseVector<IndexType> const &coordinates, scai::dmemo::DistributionPtr const dist_wavefield)
        {
            SCAI_ASSERT_DEBUG(coordinates.size() > 0, "The vector coordinates does not contain any elements! ");

            scai::hmemo::HArray<IndexType> localIndices;

            Global2Local(coordinates, localIndices, dist_wavefield);

            scai::dmemo::DistributionPtr dist_temp(new scai::dmemo::GeneralDistribution(coordinates.size(), localIndices, true, dist_wavefield->getCommunicatorPtr()));

            return (dist_temp);
        }

        /*! \brief Getter method for source coordinates from sourceSettings 
        *
        \param sourceSettings sourceSettings
        */
        template <typename ValueType>
        scai::lama::DenseVector<IndexType> getsourcecoordinates(std::vector<sourceSettings<ValueType>> &sourceSettings, Coordinates<ValueType> const &modelCoordinates)
        {
            scai::lama::DenseVector<IndexType> sourcecoords1D;
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
        inline void getCutCoord(Configuration::Configuration const &config, std::vector<KITGPI::Acquisition::coordinate3D> &cutCoordinates, std::vector<Acquisition::sourceSettings<ValueType>> sourceSettingsBig, Coordinates<ValueType> const &modelCoordinates, Coordinates<ValueType> const &modelCoordinatesBig)
        {
            cutCoordinates.clear();
            std::vector<IndexType> uniqueShotNos;
            Acquisition::calcuniqueShotNo(uniqueShotNos, sourceSettingsBig);
            Acquisition::coordinate3D coordinate;
            SCAI_ASSERT(sourceSettingsBig.size() == uniqueShotNos.size(), "sourceSettingsBig.size() != uniqueShotNos.size()");
            IndexType numshotsIncr = sourceSettingsBig.size();
            ValueType x0 = modelCoordinates.getX0();
            ValueType x0Big = modelCoordinatesBig.getX0();
            ValueType DH = modelCoordinates.getDH();
            ValueType DHBig = modelCoordinatesBig.getDH();
            SCAI_ASSERT(x0 == x0Big, "x0 != x0Big");
            SCAI_ASSERT(DH == DHBig, "DH != DHBig");
            IndexType sourceCoordsX = 0;
            IndexType useSourceEncode = config.getAndCatch("useSourceEncode", 0);
            IndexType numShotDomains = config.get<IndexType>("NumShotDomains"); // the number of supershot
            IndexType numShotPerSuperShot = ceil(ValueType(numshotsIncr) / numShotDomains); 
                
            IndexType minX = sourceSettingsBig[0].sourceCoords.x;
            for (int i = 0; i < numshotsIncr; i++) {               
                if (minX > sourceSettingsBig[i].sourceCoords.x)
                    minX = sourceSettingsBig[i].sourceCoords.x;
            }
            for (int i = 0; i < numshotsIncr; i++) {    
                if (sourceSettingsBig[i].sourceNo >= 0) {
                    if (useSourceEncode == 0) {
                        sourceCoordsX = sourceSettingsBig[i].sourceCoords.x;
                    } else if (useSourceEncode == 3 && i % numShotPerSuperShot == 0) {
                        sourceCoordsX = sourceSettingsBig[i].sourceCoords.x;
                    }
                } // if sourceSettingsBig[i].sourceNo < 0, the previous sourceCoordsX will be used.
                SCAI_ASSERT(sourceCoordsX != 0, "sourceCoordsX cannot be 0 when sourceSettingsBig[i].sourceNo < 0");
                coordinate.x = sourceCoordsX - minX;
                coordinate.y = 0;
                coordinate.z = 0;
                cutCoordinates.push_back(coordinate);
            }
        }
        
        /*! \brief Write to cutCoordinates-file
        *
        \param config config
        \param cutCoordinates coordinates of cutting model per shot
        \param uniqueShotNos unique shot numbers
        */
        template <typename ValueType>
        inline void writeCutCoordToFile(scai::dmemo::CommunicatorPtr comm, Configuration::Configuration const &config, std::vector<KITGPI::Acquisition::coordinate3D> cutCoordinates, std::vector<IndexType> uniqueShotNos, ValueType NXPerShot)
        {
            int myRank = comm->getRank(); 
            bool useStreamConfig = config.getAndCatch<bool>("useStreamConfig", false);
            if (useStreamConfig && myRank == MASTERGPI) {
                std::ofstream outputFile(config.get<std::string>("cutCoordinatesFilename"));
                outputFile << "# Coordinate for cutting model per shot (the first row is the size of modelPerShot)\n";
                outputFile << "# ShotNumber | index_x | index_y | index_z\n";  
                outputFile << std::setw(12) << NXPerShot*config.get<IndexType>("NY")*config.get<IndexType>("NZ") << std::setw(10) << NXPerShot << std::setw(10) << config.get<IndexType>("NY") << std::setw(10) << config.get<IndexType>("NZ") << "\n";
                for (unsigned i = 0; i < uniqueShotNos.size(); i++) {  
                    outputFile << std::setw(12) << uniqueShotNos[i] << std::setw(10) << cutCoordinates[i].x << std::setw(10) << cutCoordinates[i].y << std::setw(10) << cutCoordinates[i].z << "\n";
                }
                outputFile.close();
            }
        }
        
        /*! \brief Generate a random shot sequence without any repeated elements and with the max appearances of one shot
        *
        \param uniqueShotInds unique shot indexes
        \param shotHistory a vector to count how much times the shot appears
        \param numshots numbers of shots
        \param maxcount max times for one shot
        \param useRandomSource useRandomSource
        */
        template <typename ValueType>
        void getRandomShotInds(std::vector<IndexType> &uniqueShotInds, std::vector<IndexType> &shotHistory, IndexType numshots, IndexType maxcount, IndexType useRandomSource, IndexType &seedtime)
        {  
            IndexType numShotDomains = uniqueShotInds.size();
            IndexType randomShotInd = 0;
              
            if (useRandomSource == 1) {
                std::vector<IndexType> randomShotIndHistory(numShotDomains, 0);
                std::srand(seedtime);
                seedtime++;
                for (IndexType shotDomainInd = 0; shotDomainInd < numShotDomains; shotDomainInd++) {                                     
                    bool repeat = false;
                    randomShotInd = std::rand() % numshots;
                    randomShotIndHistory[shotDomainInd] = randomShotInd;
                    for (IndexType i = 0; i < shotDomainInd; i++) {
                        if (randomShotIndHistory[i] == randomShotInd) {
                            repeat = true;
                            break;
                        }
                    }
                    if (shotHistory[randomShotInd] >= maxcount || repeat) {
                        shotDomainInd--;
                    } else {
                        uniqueShotInds[shotDomainInd] = randomShotInd;
                        shotHistory[randomShotInd]++;
                    }
                }
            } else if (useRandomSource == 2) { // sequentially to cover global area
                IndexType sum = accumulate(shotHistory.begin(), shotHistory.end(), 0);
                IndexType step = ceil(ValueType(numshots) / numShotDomains);
                sum /= numShotDomains; // iteration
                for (IndexType shotDomainInd = 0; shotDomainInd < numShotDomains; shotDomainInd++) { 
                    randomShotInd = sum + shotDomainInd * step;
                    randomShotInd = randomShotInd % numshots;
                    uniqueShotInds[shotDomainInd] = randomShotInd;
                    shotHistory[randomShotInd]++;
                }
            } else if (useRandomSource == 3) { // sequentially to cover local area
                IndexType sum = accumulate(shotHistory.begin(), shotHistory.end(), 0);
                sum /= numShotDomains; // iteration
                for (IndexType shotDomainInd = 0; shotDomainInd < numShotDomains; shotDomainInd++) { 
                    randomShotInd = sum * numShotDomains + shotDomainInd;
                    randomShotInd = randomShotInd % numshots;
                    uniqueShotInds[shotDomainInd] = randomShotInd;
                    shotHistory[randomShotInd]++;
                }
            }
        }
                
        /*! \brief Write to randomSource-file
        *
        \param comm Communicator
        \param logFilename Name of log-file
        \param uniqueShotNos unique shot numbers
        \param uniqueShotInds unique shot indexes
        \param stage inversion stage
        \param iteration inversion iteration
        \param useRandomSource useRandomSource
        */
        inline void writeRandomShotNosToFile(scai::dmemo::CommunicatorPtr comm, std::string logFilename, std::vector<IndexType> uniqueShotNos, std::vector<IndexType> uniqueShotInds, IndexType stage, IndexType iteration, IndexType useRandomSource)
        {      
            int myRank = comm->getRank();  
            if (useRandomSource != 0 && myRank == MASTERGPI) {
                std::ofstream outputFile; 
                std::string randomSourceFilename = logFilename.substr(0, logFilename.length()-4) + ".randomSource" + logFilename.substr(logFilename.length()-4, 4);
                if (stage == 1 && iteration == 0) {
                    outputFile.open(randomSourceFilename);
                    outputFile << "# Shot number records during inversion\n"; 
                    outputFile << "# random source type = " << useRandomSource << " (0=all sequential shot, 1=numShotDomains random shot, 2=numShotDomains sequential shot)\n"; 
                    outputFile << "# Stage | Iteration | shot number\n"; 
                } else {                    
                    outputFile.open(randomSourceFilename, std::ios_base::app);
                    outputFile << std::scientific;
                }
                outputFile << std::setw(5) << stage << std::setw(10) << iteration;
                for (unsigned i = 0; i < uniqueShotInds.size(); i++) {  
                    if (i == 0) {
                        outputFile << std::setw(9) << uniqueShotNos[uniqueShotInds[i]];
                    } else {
                        outputFile << std::setw(4) << uniqueShotNos[uniqueShotInds[i]];
                    }
                }
                outputFile << "\n";
                outputFile.close();
            }
        }
    }
}
