#pragma once

#include <scai/dmemo/SingleDistribution.hpp>
#include <scai/hmemo.hpp>
#include <scai/lama.hpp>

#include "Acquisition.hpp"
#include "Coordinates.hpp"
#include "segy.hpp"
#include "AcquisitionSettings.hpp"

namespace KITGPI
{

    namespace Acquisition
    {

        //! Handling of I/O in SU format
        /*!
         * This class handels reading and writing to SU format.
         */
        template <typename ValueType>
        class suHandler
        {

          public:
            //! Default constructor
            suHandler() : nShots(Acquisition::NUM_ELEMENTS_SEISMOGRAMTYPE){};

            //! Default destructor
            ~suHandler(){};

            // Acquisition Matrix building
            typedef void (suHandler<ValueType>::*buildAcqMatPtr)(std::string const &, scai::lama::DenseMatrix<ValueType> &, ValueType); // function pointer to the functions which are needed to build acquisition matrix for one component
            void buildAcqMatrixSource(std::string const &filename, ValueType DH);
            void buildAcqMatrixReceiver(std::string const &filename, ValueType DH);

            // utility functions
            void locateTrace(std::string &filename, scai::IndexType &traceNumber, scai::IndexType shotNumber);
            std::vector< KITGPI::Acquisition::sourceSettings< ValueType > >& getSourceSettingsVec ();
            std::vector<receiverSettings>  &getReceiverSettingsVec();
            KITGPI::Acquisition::sourceSettings< ValueType >& getSourceSettings (scai::IndexType shotNumber );

            // following member functions are static because they are associated with the class rather than the objects
            static void readDataSU(std::string const &filename, scai::lama::Matrix<ValueType> &data, scai::IndexType ns, scai::IndexType ntr);
            static void readSingleDataSU(std::string const &filename, scai::lama::Vector<ValueType> &data, scai::IndexType traceNumber);
            static void readHeaderSU(std::string const &filename, std::vector<Segy> &header);

            static void writeSU(std::string const &filename, scai::lama::DenseMatrix<ValueType> const &data, scai::lama::DenseVector<scai::IndexType> const &coordinates, ValueType DT, scai::IndexType sourceCoordinate, Coordinates<ValueType> const &modelCoordinates);

          private:
            void buildAcqMatrixSourceComp(std::string const &filename, std::vector<sourceSettings<ValueType>> &sourceSettingsVec, ValueType DH);
            void buildAcqMatrixReceiverComp(std::string const &filename, std::vector<receiverSettings> &receiverSettingsVec, ValueType DH);

            scai::IndexType getComponentFromName(std::string const &filename);
            void indexShots(std::string const &filename);
            static void initSegy(Segy &tr);

            std::vector<sourceSettings<ValueType>> acqSource; // acquisition matrix
            std::vector<receiverSettings> acqReceiver;
            std::vector<scai::IndexType> nShots;       // number of shots per component which is needed to localize a single shot in multiple SU files
        };
    }
}
