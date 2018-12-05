
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
            suHandler() : nShots(Acquisition::NUM_ELEMENTS_SEISMOGRAMTYPE){};

            //! Default destructor
            ~suHandler(){};

            // Acquisition Matrix building
            typedef void (suHandler<ValueType>::*buildAcqMatPtr)(std::string const &, scai::lama::DenseMatrix<ValueType> &, ValueType); // function pointer to the functions which are needed to build acquisition matrix for one component
            void buildAcqMatrix(std::string const &filename, ValueType DH, buildAcqMatPtr buildAcqMat);
            void buildAcqMatrixSource(std::string const &filename, ValueType DH);
            void buildAcqMatrixReceiver(std::string const &filename, ValueType DH);

            // utility functions
            void locateTrace(std::string &filename, scai::IndexType &traceNumber, scai::IndexType shotNumber);
            scai::lama::DenseMatrix<ValueType> const &getAcquisition() const;
            void getAcquisitionRow(scai::lama::DenseMatrix<ValueType> &acqRowMat, scai::IndexType shotNumber) const;

          private:
            void buildAcqMatrixSourceComp(std::string const &filename, scai::lama::DenseMatrix<ValueType> &acqMatTmp, ValueType DH);
            void buildAcqMatrixReceiverComp(std::string const &filename, scai::lama::DenseMatrix<ValueType> &acqMatTmp, ValueType DH);

            scai::IndexType getComponentFromName(std::string const &filename);
            void indexShots(std::string const &filename);

            scai::lama::DenseMatrix<ValueType> acqMat; // acquisition matrix
            std::vector<scai::IndexType> nShots;       // number of shots per component which is needed to localize a single shot in multiple SU files
        };

        // free functions for SU I/O
        template <typename ValueType>
        void readDataSU(std::string const &filename, scai::lama::Matrix<ValueType> &data, scai::IndexType ns, scai::IndexType ntr);
        template <typename ValueType>
        void readSingleDataSU(std::string const &filename, scai::lama::Vector<ValueType> &data, scai::IndexType traceNumber);
        template <typename ValueType>
        void readHeaderSU(std::string const &filename, std::vector<Segy> &header);
        template <typename ValueType>
        void writeSU(std::string const &filename, scai::lama::DenseMatrix<ValueType> const &data, scai::lama::DenseVector<scai::IndexType> const &coordinates, ValueType DT, scai::IndexType sourceCoordinate, scai::IndexType NX, scai::IndexType NY, scai::IndexType NZ, ValueType DH);
        template <typename ValueType>
        void initSegy(Segy &tr);

        //! \brief Initialize a Segy struct
        /*!
        \param segy Segy struct
        */
        template <typename ValueType>
        void initSegy(Segy &tr)
        {
            /* Define parameters in tr header */
            tr.tracr = 0; /* trace sequence number within reel */
            tr.fldr = 0;  /* field record number */
            tr.tracf = 0; /* trace number within field record */
            tr.ep = 0;    /* energy source point number */
            tr.cdpt = 0;

            /* trace number within CDP ensemble */
            tr.nvs = 0;    /* number of vertically summed traces (see vscode
                                    in bhed structure) */
            tr.nhs = 0;    /* number of horizontally summed traces (see vscode
                                    in bhed structure) */
            tr.duse = 0;   /* data use:
                                    1 = production
                                    2 = test */
            tr.selev = 0;  /* source elevation from sea level
                                (above sea level is positive) */
            tr.gdel = 0;   /* datum elevation at receiver group */
            tr.sdel = 0;   /* datum elevation at source */
            tr.gwdep = 0;  /* water depth at receiver group */
            tr.sy = 0;     /* Y source coordinate */
            tr.gy = 0;     /* Y group coordinate */
            tr.counit = 1; /* coordinate units code:
                                    for previous four entries
                                    1 = length (meters or feet)
                                    2 = seconds of arc (in this case, the
                                    X values are longitude and the Y values
                                    are latitude, a positive value designates
                                    the number of seconds east of Greenwich
                                    or north of the equator */
            tr.wevel = 0;  /* weathering velocity */
            tr.swevel = 0; /* subweathering velocity */
            tr.sut = 0;    /* uphole time at source */
            tr.gut = 0;    /* uphole time at receiver group */
            ;

            tr.sstat = 0;  /* source static correction */
            tr.gstat = 0;  /* group static correction */
            tr.tstat = 0;  /* total static applied */
            tr.laga = 0;   /* lag time A, time in ms between end of 240-
                                byte trace identification header and time
                                break, positive if time break occurs after
                                end of header, time break is defined as
                                the initiation pulse which maybe recorded
                                on an auxiliary trace or as otherwise
                                specified by the recording system */
            tr.lagb = 0;   /* lag time B, time in ms between the time break
                                and the initiation time of the energy source,
                                may be positive or negative */
            tr.delrt = 0;  /* delay recording time, time in ms between
                                initiation time of energy source and time
                                when recording of data samples begins
                                (for deep water work if recording does not
                                start at zero time) */
            tr.muts = 0;   /* mute time--start */
            tr.mute = 0;   /* mute time--end */
            tr.gain = 0;   /* gain type of field instruments code:
                                1 = fixed
                                2 = binary
                                3 = floating point
                                4 ---- N = optional use */
            tr.igc = 0;    /* instrument gain constant */
            tr.igi = 0;    /* instrument early or initial gain */
            tr.corr = 0;   /* correlated:
                                1 = no
                                2 = yes */
            tr.sfs = 0;    /* sweep frequency at start */
            tr.sfe = 0;    /* sweep frequency at end */
            tr.slen = 0;   /* sweep length in ms */
            tr.styp = 0;   /* sweep type code:
                                1 = linear
                                2 = cos-squared
                                3 = other */
            tr.stas = 0;   /* sweep trace length at start in ms */
            tr.stae = 0;   /* sweep trace length at end in ms */
            tr.tatyp = 0;  /* taper type: 1=linear, 2=cos^2, 3=other */
            tr.afilf = 0;  /* alias filter frequency if used */
            tr.afils = 0;  /* alias filter slope */
            tr.nofilf = 0; /* notch filter frequency if used */
            tr.nofils = 0; /* notch filter slope */
            tr.lcf = 0;    /* low cut frequency if used */
            tr.hcf = 0;    /* high cut frequncy if used */
            tr.lcs = 0;    /* low cut slope */
            tr.hcs = 0;    /* high cut slope */
            tr.year = 0;   /* year data recorded */
            tr.day = 0;    /* day of year */
            tr.hour = 0;   /* hour of day (24 hour clock) */
            tr.minute = 0; /* minute of hour */
            tr.sec = 0;    /* second of minute */
            tr.timbas = 0; /* time basis code:
                                1 = local
                                2 = GMT
                                3 = other */
            tr.trwf = 0;   /* trace weighting factor, defined as 1/2^N
                                volts for the least sigificant bit */
            tr.grnors = 0; /* geophone group number of roll switch
                                position one */
            tr.grnofr = 0; /* geophone group number of trace one within
                                original field record */
            tr.grnlof = 0; /* geophone group number of last trace within
                                original field record */
            tr.gaps = 0;   /* gap size (total number of groups dropped) */
            tr.otrav = 0;  /* overtravel taper code:
                                1 = down (or behind)
                                2 = up (or ahead) */

            /* local assignments */

            tr.f1 = 0.0; /* first sample location for non-seismic data */

            tr.d2 = 0.0; /* sample spacing between traces */

            tr.f2 = 0.0; /* first trace location */

            tr.ungpow = 0.0;  /* negative of power used for dynamic
                                range compression */
            tr.unscale = 0.0; /* reciprocal of scaling factor to normalize
                                range */
            tr.mark = 0;
        }

        //! \brief Read a SU file from disk without header
        /*!
        *
        \param filename Filename to read from
        \param data Matrix where the data read is stored in
        \param ns number of samples in trace
        \param ntr number of traces
        */
        template <typename ValueType>
        void readDataSU(std::string const &filename, scai::lama::Matrix<ValueType> &data, scai::IndexType ns, scai::IndexType ntr)
        {
            Segy tr;

            const char *filetemp = filename.c_str();
            FILE *pFile;
            pFile = fopen(filetemp, "rb");

            scai::hmemo::HArray<float> dataTmp;
            scai::lama::DenseVector<float> traceTmp;
            scai::lama::DenseVector<ValueType> trace;

            for (scai::IndexType tracl1 = 0; tracl1 < ntr; tracl1++) {
                fseek(pFile, 240, SEEK_CUR);
                fread(&tr.data[0], 4, ns, pFile);

                dataTmp.setRawData(ns, tr.data);
                traceTmp.assign(dataTmp);
                trace = scai::lama::cast<ValueType, float>(traceTmp);
                data.setRow(trace, tracl1, scai::common::BinaryOp::COPY);
            }
        }

        //! \brief Read a single trace without header form SU
        /*!
        *
        \param filename Filename to read from
        \param data Vector where the data read is stored in
        \param traceNumber number of trace to read
        */
        template <typename ValueType>
        void readSingleDataSU(std::string const &filename, scai::lama::Vector<ValueType> &data, scai::IndexType traceNumber)
        {
            Segy tr;

            const char *filetemp = filename.c_str();
            FILE *pFile;
            pFile = fopen(filetemp, "rb");
            fread(&tr, 1, 240, pFile);

            scai::hmemo::HArray<float> dataTmp;
            scai::lama::DenseVector<float> traceTmp;

            scai::IndexType nSkip = (240 + tr.ns * 4) * traceNumber;
            fseek(pFile, nSkip, SEEK_CUR);
            fread(&tr.data[0], 4, tr.ns, pFile);

            dataTmp.setRawData(tr.ns, tr.data);
            traceTmp.assign(dataTmp);
            data = scai::lama::cast<ValueType, float>(traceTmp);
        }

        //! \brief Write a seismogram to disk in Seismic Unix (SEG-Y) format
        /*!
        *
        * This method writes the seismogram in the Seismic Unix format to disk.
        * Some header information will be calculated based on the input parameters and will be included in the seismic unix file.
        \param filename Filename to write seismogram in Seismic Unix (SEG-Y) format
        \param NX Number of grid points in X direction
        \param NY Number of grid points in Y direction
        \param NZ Number of grid points in Z direction
        \param DH Length of space step in meter
        */
        template <typename ValueType>
        void writeSU(std::string const &filename, scai::lama::DenseMatrix<ValueType> const &data, scai::lama::DenseVector<scai::IndexType> const &coordinates, ValueType DT, scai::IndexType sourceCoordinate, scai::IndexType NX, scai::IndexType NY, scai::IndexType NZ, ValueType DH)
        {
            Segy tr;
            initSegy<ValueType>(tr);

            int tracl1;
            double temp3;
            scai::IndexType temp2;
            float xr, yr, zr, x, y, z;
            float XS = 0.0, YS = 0.0, ZS = 0.0;
            const float xshift = 800.0, yshift = 800.0;
            float dtms = float(DT * 1000000);

            scai::IndexType ns = data.getNumColumns();
            scai::IndexType ntr = data.getNumRows();
            tr.ntr = ntr; /* number of traces */

            const char *filetemp = filename.c_str();
            FILE *pFile;
            pFile = fopen(filetemp, "wb");
            scai::lama::DenseVector<ValueType> tempdata;

            Coordinates coordTransform(NX, NY, NZ);
            coordinate3D coord3Dsrc;
            coordinate3D coord3Drec;
            coord3Dsrc = coordTransform.index2coordinate(sourceCoordinate);

            YS = coord3Dsrc.y;
            XS = coord3Dsrc.x;
            ZS = coord3Dsrc.z;
            YS = YS * DH;
            XS = XS * DH;
            ZS = ZS * DH;

            for (tracl1 = 0; tracl1 < ntr; tracl1++) {
                temp3 = float(coordinates.getValue(tracl1));
                temp2 = floor(temp3);
                coord3Drec = coordTransform.index2coordinate(temp2);
                xr = coord3Drec.x;
                yr = coord3Drec.y;
                zr = coord3Drec.z;
                yr = yr * DH;
                xr = xr * DH;
                zr = zr * DH;
                x = xr - XS; // Taking source position as reference point
                y = yr - YS;
                z = zr - ZS;

                tr.tracl = tracl1 + 1; // trace sequence number within line
                tr.tracr = 1;          // trace sequence number within reel
                tr.ep = 1;
                tr.cdp = ntr;
                tr.trid = (short)1;
                tr.offset = (signed int)round(sqrt((XS - xr) * (XS - xr) + (YS - yr) * (YS - yr) + (ZS - zr) * (ZS - zr)) * 1000.0);
                tr.gelev = (signed int)round(yr * 1000.0);
                tr.sdepth = (signed int)round(YS * 1000.0); /* source depth (positive) */
                /* angle between receiver position and reference point
                    (sperical coordinate system: swdep=theta, gwdep=phi) */
                tr.gdel = (signed int)round(atan2(-y, z) * 180 * 1000.0 / 3.1415926);
                tr.gwdep = (signed int)round(sqrt(z * z + y * y) * 1000.0);
                tr.swdep = round(((360.0 / (2.0 * 3.1415926)) * atan2(x - xshift, y - yshift)) * 1000.0);
                tr.scalel = (signed short)-3;
                tr.scalco = (signed short)-3;
                tr.sx = (signed int)round(XS * 1000.0); /* X source coordinate */
                tr.sy = (signed int)round(ZS * 1000.0); /* Y source coordinate */

                /* group coordinates */
                tr.gx = (signed int)round(xr * 1000.0);
                tr.gy = (signed int)round(zr * 1000.0);
                tr.ns = (unsigned short)ns;          /* number of samples in this trace */
                tr.dt = (unsigned short)round(dtms); /* sample interval in micro-seconds */
                tr.d1 = (float)tr.dt * 1.0e-6;       /* sample spacing for non-seismic data */

                data.getRow(tempdata, tracl1);
                for (scai::IndexType sample = 0; sample < tempdata.size(); sample++) {
                    tr.data[sample] = float(tempdata.getValue(sample));
                }

                fwrite(&tr, 240, 1, pFile);
                fwrite(&tr.data[0], 4, ns, pFile);
            }
            fclose(pFile);
        }

        //! \brief Read all headers of a SU file and store them in a standard vector
        /*!
        \param filename Name of the file
        \param header std::vecor the headers are stored in
        */
        template <typename ValueType>
        void readHeaderSU(std::string const &filename, std::vector<Segy> &header)
        {
            Segy tr;

            header.clear();

            const char *filetemp = filename.c_str();
            FILE *pFile;
            pFile = fopen(filetemp, "rb");

            if (pFile) {
                int c;
                while ((c = fgetc(pFile)) != EOF) {
                    ungetc(c, pFile);
                    fread(&tr, 1, 240, pFile);
                    header.push_back(tr);
                    fseek(pFile, 4 * tr.ns, SEEK_CUR);
                }
            }
        }
    }
}