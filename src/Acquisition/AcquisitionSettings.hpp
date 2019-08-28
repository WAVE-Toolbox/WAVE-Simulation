#pragma once

#include "Acquisition.hpp"
#include "Coordinates.hpp"
#include <fstream>
#include <scai/lama.hpp>

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
            coordinate3D getCoords() {return sourceCoords;}
            scai::IndexType getType() {return sourceType;}
            scai::IndexType row;
        };
        
        /*! \brief Struct to save 3-D coordinates
         *
         * This struct saves the coordinates and parameters of a receiver
         */
        struct receiverSettings {
            coordinate3D receiverCoords;
            scai::IndexType receiverType;
            coordinate3D getCoords() {return receiverCoords;}
            scai::IndexType getType() {return receiverType;}
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
                    if (firstChar == '#') {
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
                //       std::cout << "Source acquisition file (" << fileName << ") read in." << std::endl;
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
                    if (firstChar == '#') {
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
                //  std::cout << "Receiver acquisition file (" << fileName << ") read in." << std::endl;
            } else {
                COMMON_THROWEXCEPTION("Could not open receiver acquisition file " << fileName)
            }
        }
    }
}
