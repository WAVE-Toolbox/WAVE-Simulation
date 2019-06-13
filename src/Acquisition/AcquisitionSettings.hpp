#pragma once

#include <scai/lama.hpp>
#include "Coordinates.hpp"
#include "Acquisition.hpp"
#include <fstream>

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
            coordinate3D sourceCoords;
            scai::IndexType sourceType;
            scai::IndexType waveletType;
            scai::IndexType waveletShape;
            ValueType fc;
            ValueType amp;
            ValueType tShift;
            coordinate3D getCoords() {return sourceCoords;}
            scai::IndexType getType() {return sourceType;}
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
        
        /*! \brief Read the settings of a source into sourceSettings
        *
        \param settings sourceSettings struct
        \param fileName Name of source file
        \param i number of source to read (which line of source file)
        */
        template <typename ValueType>
        inline void readSettings(sourceSettings<ValueType> &settings, std::string fileName, scai::IndexType i) {
            std::ifstream istream(fileName, std::ios::binary);
            if (istream.is_open()) {
                istream.seekg(i * sizeof(settings));
                istream >> settings.sourceCoords.x;
                istream >> settings.sourceCoords.y;
                istream >> settings.sourceCoords.z;
                istream >> settings.sourceType;
                istream >> settings.waveletType;
                istream >> settings.waveletShape;
                istream >> settings.fc;
                istream >> settings.amp;
                istream >> settings.tShift;
            }
            else {
                COMMON_THROWEXCEPTION("Could not open source acquisition file")
                
            }
        }
        
        inline void readSettings(receiverSettings &settings, std::string fileName, scai::IndexType i) {
            std::ifstream istream(fileName, std::ios::binary);
            if (istream.is_open()) {
                istream.seekg(i * sizeof(settings));
                istream >> settings.receiverCoords.x;
                istream >> settings.receiverCoords.y;
                istream >> settings.receiverCoords.z;
                istream >> settings.receiverType;
            }
            else {
                COMMON_THROWEXCEPTION("Could not open receiver acquisition file")
            }
        }
        
        /*! \brief Read the settings of a source into sourceSettingssd
        *
        \param settings sourceSettings struct
        \param istream Input file stream
        */
        template <typename ValueType>
        inline void readSettings(sourceSettings<ValueType> &settings, std::ifstream &istream) {
            istream >> settings.sourceCoords.x;
            istream >> settings.sourceCoords.y;
            istream >> settings.sourceCoords.z;
            istream >> settings.sourceType;
            istream >> settings.waveletType;
            istream >> settings.waveletShape;
            istream >> settings.fc;
            istream >> settings.amp;
            istream >> settings.tShift;
        }
        
        /*! \brief Read the settings of a receiver into receiverSettings
        *
        \param settings receiverSettings struct
        \param istream Input file stream
        */
        inline void readSettings(receiverSettings &settings, std::ifstream &istream) {
            istream >> settings.receiverCoords.x;
            istream >> settings.receiverCoords.y;
            istream >> settings.receiverCoords.z;
            istream >> settings.receiverType;
        }
        
        template <typename ValueType>
        inline void readAllSettings(std::vector<sourceSettings<ValueType>> &allSettings, std::string fileName) {
            std::ifstream istream(fileName, std::ios::binary);
            allSettings.clear();
            sourceSettings<ValueType> thisSettings;
            if (istream.is_open()) {
                while (!istream.eof()) {
                    readSettings(thisSettings, istream);
                    allSettings.push_back(thisSettings);
                }
                allSettings.pop_back();
                //std::cout << fileName << " " << allSettings.size() << std::endl;
            }
            else {
                COMMON_THROWEXCEPTION("Could not open receiver acquisition file")
            }
        }
        
        inline void readAllSettings(std::vector<receiverSettings> &allSettings, std::string fileName) {
            std::ifstream istream(fileName, std::ios::binary);
            allSettings.clear();
            receiverSettings thisSettings;
            if (istream.is_open()) {
                while (!istream.eof()) {
                    readSettings(thisSettings, istream);
                    allSettings.push_back(thisSettings);
                }
                allSettings.pop_back();
                //std::cout << fileName << " " << allSettings.size() << std::endl;
            }
            else {
                COMMON_THROWEXCEPTION("Could not open receiver acquisition file")
            }
        }
    }
}


