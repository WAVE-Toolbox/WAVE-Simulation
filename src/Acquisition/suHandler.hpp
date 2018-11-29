
#pragma once

#include <scai/dmemo/SingleDistribution.hpp>
#include <scai/hmemo.hpp>
#include <scai/lama.hpp>

#include "Acquisition.hpp"
#include "Coordinates.hpp"
#include "segy.hpp"

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
            suHandler(){};

            //! Default destructor
            ~suHandler(){};

            void readDataSU(std::string const &filename, scai::lama::Matrix<ValueType> &data, scai::IndexType ns, scai::IndexType ntr);
            void readSingleDataSU(std::string const &filename, scai::lama::Vector<ValueType> &data, scai::IndexType traceNumber);
            void readHeaderSU(std::string const &filename, std::vector<Segy> &header);
            void writeSU(std::string const &filename, scai::lama::DenseMatrix<ValueType> const &data, scai::lama::DenseVector<scai::IndexType> const &coordinates, ValueType DT, scai::IndexType sourceCoordinate, scai::IndexType NX, scai::IndexType NY, scai::IndexType NZ, ValueType DH) const;

            typedef void (suHandler<ValueType>::*buildAcqMatPtr)(std::string const &, scai::lama::DenseMatrix<ValueType> &, ValueType); // function pointer to the functions which are needed to build acquisition matrix for one component
            void buildAcqMatrix(std::string const &filename, ValueType DH, buildAcqMatPtr buildAcqMat);
            void buildAcqMatrixSource(std::string const &filename, ValueType DH);
            void buildAcqMatrixReceiver(std::string const &filename, ValueType DH);

            void locateTrace(std::string &filename, scai::IndexType &traceNumber, scai::IndexType shotNumber);
            scai::lama::DenseMatrix<ValueType> const &getAcquisition() const;
            void getAcquisitionRow(scai::lama::DenseMatrix<ValueType> &acqRowMat, scai::IndexType shotNumber) const;

          private:
            void buildAcqMatrixSourceComp(std::string const &filename, scai::lama::DenseMatrix<ValueType> &acqMatTmp, ValueType DH);
            void buildAcqMatrixReceiverComp(std::string const &filename, scai::lama::DenseMatrix<ValueType> &acqMatTmp, ValueType DH);

            void initSegy(Segy &tr) const;
            scai::IndexType getComponentFromName(std::string const &filename);
            void indexShots(std::string const &filename);

            scai::lama::DenseMatrix<ValueType> acqMat;
            scai::IndexType nShotsP;
            scai::IndexType nShotsVX;
            scai::IndexType nShotsVY;
            scai::IndexType nShotsVZ;
        };
    }
}