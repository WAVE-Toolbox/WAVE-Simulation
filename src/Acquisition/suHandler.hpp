#pragma once

#include <scai/dmemo/SingleDistribution.hpp>
#include <scai/hmemo.hpp>
#include <scai/lama.hpp>

#include "Acquisition.hpp"
#include "AcquisitionSettings.hpp"
#include "Coordinates.hpp"
#include "segy.hpp"

namespace KITGPI
{

    namespace Acquisition
    {

        //! Handling of I/O in SU format
        /*!
         * This class handels acquisition vectors from su header
         */
        template <typename ValueType>
        class suHandler
        {

          public:
            //! Default constructor
            suHandler(){};

            //! Default destructor
            ~suHandler(){};

            void readAllSettingsFromSU(std::vector<sourceSettings<ValueType>> &allSettings, std::string const &filename, ValueType DH);
            void readAllSettingsFromSU(std::vector<receiverSettings> &allSettings, std::string const &filename, ValueType DH);

          private:
            void readAllSettingsFromSUComp(std::string const &filename, std::vector<sourceSettings<ValueType>> &sourceSettingsVec, ValueType DH);
            void readAllSettingsFromSUComp(std::string const &filename, std::vector<receiverSettings> &receiverSettingsVec, ValueType DH);

            scai::IndexType getComponentFromName(std::string const &filename);
        };
    }
}
