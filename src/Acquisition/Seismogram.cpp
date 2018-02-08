#include "Seismogram.hpp"
using namespace scai;

//! \brief copy constructor
/*!
 *
 \param rhs seismogram to copy
 */
template <typename ValueType>
KITGPI::Acquisition::Seismogram<ValueType>::Seismogram(const Seismogram &rhs)
{
    numSamples = rhs.numSamples;
    numTracesGlobal = rhs.numTracesGlobal;
    numTracesLocal = rhs.numTracesLocal;
    DT = rhs.DT;
    type = rhs.type;
    coordinates = rhs.coordinates;
    sourceCoordinate = rhs.sourceCoordinate;
    data = rhs.data;
}
//! \brief swap function
/*!
 *
 * swaps all members from rhs and lhs
 \param rhs seismogram to swap with
 */
template <typename ValueType>
void KITGPI::Acquisition::Seismogram<ValueType>::swap(KITGPI::Acquisition::Seismogram<ValueType> &rhs)
{
    std::swap(numSamples, rhs.numSamples);
    std::swap(numTracesGlobal, rhs.numTracesGlobal);
    std::swap(numTracesLocal, rhs.numTracesLocal);
    std::swap(DT, rhs.DT);
    std::swap(type, rhs.type);
    std::swap(coordinates, rhs.coordinates);
    std::swap(sourceCoordinate, rhs.sourceCoordinate);
    data.swap(rhs.data);
}

//! \brief Adding ending to the seismogram-filename-string
/*!
 *
 * This member function adds the #SeismogramType to the filname ending.
 \param filename Filename of output
 */
template <typename ValueType>
std::string KITGPI::Acquisition::Seismogram<ValueType>::addSeismogramTypeToName(std::string const &filename) const
{
    SCAI_ASSERT_DEBUG(type >= 0 && type <= (NUM_ELEMENTS_SEISMOGRAMTYPE - 1), "Wrong Trace Type: " << getTraceType());
    SCAI_ASSERT_DEBUG(strcmp(SeismogramTypeString[SeismogramType::P], "p") == 0, "Error in mapping of SeismogramType to std::string");
    SCAI_ASSERT_DEBUG(strcmp(SeismogramTypeString[SeismogramType::VX], "vx") == 0, "Error in mapping of SeismogramType to std::string");
    SCAI_ASSERT_DEBUG(strcmp(SeismogramTypeString[SeismogramType::VY], "vy") == 0, "Error in mapping of SeismogramType to std::string");
    SCAI_ASSERT_DEBUG(strcmp(SeismogramTypeString[SeismogramType::VZ], "vz") == 0, "Error in mapping of SeismogramType to std::string");

    std::size_t found = filename.find_last_of(".");
    std::string beforeEnding = filename.substr(0, found);
    std::string afterEnding = filename.substr(found);
    std::string traceTypeString = SeismogramTypeString[getTraceType()];
    return (beforeEnding + "." + traceTypeString + afterEnding);
}

//! \brief Setter method for the context ptr
/*!
 *
 * This method sets the Context to the coordinates and to the seismogram data.
 \param ctx Set ContextPtr
 */
template <typename ValueType>
void KITGPI::Acquisition::Seismogram<ValueType>::setContextPtr(scai::hmemo::ContextPtr ctx)
{
    data.setContextPtr(ctx);
    coordinates.setContextPtr(ctx);
}

//! \brief Write the seismogram to disk
/*!
 *
 * This method writes the seismogram data to disk. It uses the Configuration class to determine the filename and some requiered header information.\n
 * The output format is determined by the input parameter `SeismogramFormat`, which will be requested to the Configuration class. \n
 * **Supported Formats:**\n
 * 1. MTX: MatrixMaker format
 * 2. SU: SeismicUnix format
 \param config Configuration class which is used to determine the filename and header information
 */
template <typename ValueType>
void KITGPI::Acquisition::Seismogram<ValueType>::write(Configuration::Configuration const &config, std::string const &filename) const
{
    switch (config.get<IndexType>("SeismogramFormat")) {
    case 1:
        writeToFileRaw(filename + ".mtx");
        break;
    case 2:
        writeToFileSU(filename + ".SU", config.get<IndexType>("NX"), config.get<IndexType>("NY"), config.get<IndexType>("NZ"), config.get<ValueType>("DH"));
        break;
    default:
        COMMON_THROWEXCEPTION(" Unkown SeismogramFormat ")
        break;
    }
}

//! \brief Normalize the seismogram-traces
/*!
 *
 * This methode normalized the traces of the seismogram after the time stepping.
 */
template <typename ValueType>
void KITGPI::Acquisition::Seismogram<ValueType>::normalizeTrace()
{
    if (normalizeTraces == 1) {

        SCAI_ASSERT(data.getNumRows() == numTracesGlobal, " Size of matrix is not matching with number of traces. ");

        scai::lama::DenseVector<ValueType> tempRow;
        scai::lama::Scalar tempMax;
        scai::lama::Scalar tempInverseMax;

        for (IndexType i = 0; i < numTracesGlobal; i++) {
            tempMax = 0.0;
            tempRow.assign(0.0);
            data.getRow(tempRow, i);
            tempMax = tempRow.max();
            tempInverseMax = 1 / tempMax;
            tempRow *= tempInverseMax;
            data.setRow(tempRow, i, scai::common::binary::BinaryOp::COPY);
        }
    }
}

//! \brief Integrate the seismogram-traces
/*!
 *
 * This methode integrate the traces of the seismogram.
 */
template <typename ValueType>
void KITGPI::Acquisition::Seismogram<ValueType>::integrateTraces()
{
    SCAI_ASSERT(data.getNumRows() == numTracesGlobal, " Size of matrix is not matching with number of traces. ");

    scai::lama::DenseVector<ValueType> tempRow;
    for (IndexType i = 0; i < numTracesGlobal; i++) {
	tempRow.assign(0.0);
	data.getRow(tempRow, i);
	for (IndexType j = 0; j <tempRow.size()-1; j++) {
	      tempRow[j+1]=tempRow[j+1]*DT+tempRow[j];
	}
	data.setRow(tempRow, i, scai::common::binary::BinaryOp::COPY);
    }
}

//! \brief Setter method for the temporal sampling DT
/*!
 *
 * This method will set the temporal sampling DT to this class.
 \param newDT Temporal sampling which will be set in seconds
 */
template <typename ValueType>
void KITGPI::Acquisition::Seismogram<ValueType>::setDT(ValueType newDT)
{
    SCAI_ASSERT(newDT >= 0, " DT is smaller zero. ");
    DT = newDT;
}

//! \brief Setter method for the source coordinate
/*!
 *
 * This method sets the source coordinate to this class. The source coordinate will be used for calculation of the header information (e.g. Offset), mainly during seismogram output to disk.
 \param sourceCoord Source coordinate in 1-D format
 */
template <typename ValueType>
void KITGPI::Acquisition::Seismogram<ValueType>::setSourceCoordinate(IndexType sourceCoord)
{
    SCAI_ASSERT_DEBUG(sourceCoord >= 0, "sourceCoord is not valid");
    sourceCoordinate = sourceCoord;
}

//! \brief Setter method for the #SeismogramType
/*!
 *
 * This method will set the #SeismogramType to this class.
 \param trace Trace
 */
template <typename ValueType>
void KITGPI::Acquisition::Seismogram<ValueType>::setTraceType(SeismogramType trace)
{
    type = trace;
};

//! \brief Setter function for the coordinates of the traces
/*!
 *
 * This method will set the coordinates of the traces to this class. 
 * The size of the coordinate vector has to be equal to the number of global traces.
 \param coord DenseVector with coordinates
 */
template <typename ValueType>
void KITGPI::Acquisition::Seismogram<ValueType>::setCoordinates(scai::lama::DenseVector<IndexType> const &coord)
{
    SCAI_ASSERT_ERROR(coord.size() == numTracesGlobal, "Given traceType vector has wrong format");
    coordinates = coord;
};

//! \brief Setter methode to set Index for trace-normalization.
/*!
 *
 * This method sets the index for trace-normalization.
 \param normalizeTrace Index for trace-normalization which will normalize the seismogram traces
 */
template <typename ValueType>
void KITGPI::Acquisition::Seismogram<ValueType>::setNormalizeTraces(IndexType normalize)
{
    SCAI_ASSERT(normalize >= 0 && normalize <= 1, " Index has to be 1 or 0 ");
    normalizeTraces = normalize;
}

//! \brief Getter method for #SeismogramType
/*!
 *
 * This method returns the #SeismogramType of this seismogram.
 *
 * THIS METHOD IS CALLED DURING TIME STEPPING
 * DO NOT WASTE RUNTIME HERE
 *
 */
template <typename ValueType>
KITGPI::Acquisition::SeismogramType KITGPI::Acquisition::Seismogram<ValueType>::getTraceType() const
{
    return (type);
}

//! \brief Getter method for reference to coordinates vector
/*!
 *
 * THIS METHOD IS CALLED DURING TIME STEPPING
 * DO NOT WASTE RUNTIME HERE
 *
 */
template <typename ValueType>
lama::DenseVector<IndexType> const &KITGPI::Acquisition::Seismogram<ValueType>::getCoordinates() const
{
    SCAI_ASSERT_DEBUG(coordinates.size() == numTracesGlobal, "Size mismatch ");
    return (coordinates);
}

//! \brief Getter method for reference to seismogram data
/*!
 *
 * This method returns the DenseMatrix which is used to store the actual seismogram data.
 *
 * THIS METHOD IS CALLED DURING TIME STEPPING
 * DO NOT WASTE RUNTIME HERE
 *
 */
template <typename ValueType>
lama::DenseMatrix<ValueType> &KITGPI::Acquisition::Seismogram<ValueType>::getData()
{
    SCAI_ASSERT_DEBUG(data.getNumRows() * data.getNumColumns() == numTracesGlobal * numSamples, "Size mismatch ");
    return (data);
}

//! \brief Getter method for const reference to seismogram data
/*!
 *
 * This method returns the DenseMatrix which is used to store the actual seismogram data.
 *
 * THIS METHOD IS CALLED DURING TIME STEPPING
 * DO NOT WASTE RUNTIME HERE
 *
 */
template <typename ValueType>
lama::DenseMatrix<ValueType> const &KITGPI::Acquisition::Seismogram<ValueType>::getData() const
{
    SCAI_ASSERT_DEBUG(data.getNumRows() * data.getNumColumns() == numTracesGlobal * numSamples, "Size mismatch ");
    return (data);
}

//! \brief Replicate the seismogram data on all processes
/*!
 * Creates a copy of the seismogram data on all processe
 */
template <typename ValueType>
void KITGPI::Acquisition::Seismogram<ValueType>::replicate()
{
    dmemo::DistributionPtr no_dist_numSamples(new scai::dmemo::NoDistribution(numSamples));
    dmemo::DistributionPtr no_dist_Traces(new scai::dmemo::NoDistribution(numTracesGlobal));

    redistribute(no_dist_Traces, no_dist_numSamples);
}

//! \brief Allocate of the seismogram data
/*!
 * Allocates seismogram based on a given distribution of the traces and the number of samples per trace.
 * The data storage of the seismogram will be distributed according to distTraces. Moreover, the number of
 * local and global traces will be determined based on distTraces.
 *
 \param ctx Context
 \param distTraces Distribution for traces
 \param NT Total number of samples per trace
 */
template <typename ValueType>
void KITGPI::Acquisition::Seismogram<ValueType>::allocate(scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr distTraces, IndexType NT)
{
    SCAI_ASSERT_ERROR(NT > 0, "NT is < 0: No Seismogram allocation ");
    SCAI_ASSERT_ERROR(distTraces != NULL, "No valid distribution");

    numSamples = NT;
    numTracesGlobal = distTraces->getGlobalSize();
    numTracesLocal = distTraces->getLocalSize();

    dmemo::DistributionPtr no_dist_NT(new scai::dmemo::NoDistribution(NT));

    data.setContextPtr(ctx);
    coordinates.setContextPtr(ctx);

    data.allocate(distTraces, no_dist_NT);
    coordinates.allocate(distTraces);
}

//! \brief Reset of the seismogram data
/*!
 * This method sets the seismogra data to zero. 
 * However, the memory will stay allocated, only the content is overwriten by zeros.
 */
template <typename ValueType>
void KITGPI::Acquisition::Seismogram<ValueType>::resetData()
{
    data.scale(0.0);
}

//! \brief Redistribute the seismogram data
/*!
 *
 * Redistribution of the seismogram data according to the given distributions.
 \param distTraces Distribution of traces
 \param distSamples Distribution of temporal samples
 */
template <typename ValueType>
void KITGPI::Acquisition::Seismogram<ValueType>::redistribute(scai::dmemo::DistributionPtr distTraces, scai::dmemo::DistributionPtr distSamples)
{
    if (distSamples == NULL) {
        SCAI_ASSERT_DEBUG(numSamples >= 0, "numSamples not set");
        dmemo::DistributionPtr distSamplestmp(new scai::dmemo::NoDistribution(numSamples));
        distSamples = distSamplestmp;
    }

    data.redistribute(distTraces, distSamples);
    coordinates.redistribute(distTraces);
}

//! \brief Read a seismogram from disk without header
/*!
 *
 \param filename Filename to read seismogram
 \param distTraces Distribution of traces
 \param distSamples Distribution of temporal samples
 */
template <typename ValueType>
void KITGPI::Acquisition::Seismogram<ValueType>::readFromFileRaw(std::string const &filename, scai::dmemo::DistributionPtr distTraces, scai::dmemo::DistributionPtr distSamples)
{
    data.readFromFile(filename);
    IndexType nrow_temp = data.getNumRows();
    IndexType ncolumn_temp = data.getNumColumns();

    numSamples = ncolumn_temp;
    numTracesGlobal = nrow_temp;

    if (distTraces == NULL && distSamples == NULL) {
        replicate();
    } else {
        redistribute(distTraces, distSamples);
    }
}

//! \brief Write a seismogram to disk without header
/*!
 *
 \param filename Filename to write seismogram
 */
template <typename ValueType>
void KITGPI::Acquisition::Seismogram<ValueType>::writeToFileRaw(std::string const &filename) const
{
    if (data.getNumValues() > 0) {
        data.writeToFile(addSeismogramTypeToName(filename));
    }
}

/*! \brief Getter method for reference normalization index
 *
 *
 \return NormalizeTraces Index 
 */
template <typename ValueType>
IndexType KITGPI::Acquisition::Seismogram<ValueType>::getNormalizeTraces() const
{
    return (normalizeTraces);
}

/*! \brief Getter method for the temporal sampling
 *
 * This method returns the temporal sampling DT in seconds.
 \return DT in seconds
 */
template <typename ValueType>
ValueType KITGPI::Acquisition::Seismogram<ValueType>::getDT() const
{
    SCAI_ASSERT_ERROR(DT != 0, "Seismogramm data is not allocated");
    return (DT);
}

/*! \brief Getter method for the number of samples per trace
 *
 *
 \return The number of samples per trace
 */
template <typename ValueType>
IndexType KITGPI::Acquisition::Seismogram<ValueType>::getNumSamples() const
{
    return (numSamples);
}

/*! \brief Getter method for the number of local traces
 *
 *
 \return The number of local traces on this process
 */
template <typename ValueType>
IndexType KITGPI::Acquisition::Seismogram<ValueType>::getNumTracesLocal() const
{
    return (numTracesLocal);
}

/*! \brief Getter method for the number of global traces
*
*
\return The number of global traces of this seismogram
 */
template <typename ValueType>
IndexType KITGPI::Acquisition::Seismogram<ValueType>::getNumTracesGlobal() const
{
    return (numTracesGlobal);
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
void KITGPI::Acquisition::Seismogram<ValueType>::writeToFileSU(std::string const &filename, IndexType NX, IndexType NY, IndexType NZ, ValueType DH) const
{
    if (data.getNumValues() > 0) {
        Segy tr;
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

        int tracl1;
        lama::Scalar temp;
        double temp3;
        IndexType temp2;
        float xr, yr, zr, x, y, z;
        float XS = 0.0, YS = 0.0, ZS = 0.0;
        const float xshift = 800.0, yshift = 800.0;
        float dtms = float(DT * 1000000);

        int ns = int(numSamples);
        int ntr = int(numTracesGlobal);
        tr.ntr = ntr; /* number of traces */
        std::string tempstr = addSeismogramTypeToName(filename);
        const char *filetemp = tempstr.c_str();
        FILE *pFile;
        pFile = fopen(filetemp, "wb");
        lama::DenseVector<ValueType> tempdata;

        Coordinates<ValueType> coordTransform;
        coordinate3D coord3Dsrc;
        coordinate3D coord3Drec;
        coord3Dsrc = coordTransform.index2coordinate(sourceCoordinate, NX, NY, NZ);
        SCAI_ASSERT_DEBUG(coordTransform.index2coordinate(2, 100, 100, 100).x == 2, "")
        SCAI_ASSERT_DEBUG(coordTransform.index2coordinate(102, 100, 100, 100).y == 1, "")
        SCAI_ASSERT_DEBUG(coordTransform.index2coordinate(2, 100, 100, 1).z == 0, "")

        YS = coord3Dsrc.y;
        XS = coord3Dsrc.x;
        ZS = coord3Dsrc.z;
        YS = YS * DH;
        XS = XS * DH;
        ZS = ZS * DH;

        for (tracl1 = 0; tracl1 < ntr; tracl1++) {
            temp = coordinates.getValue(tracl1);
            temp3 = float(temp.getValue<ValueType>());
            temp2 = floor(temp3);
            coord3Drec = coordTransform.index2coordinate(temp2, NX, NY, NZ);
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
            tr.sy = (signed int)round(YS * 1000.0); /* Y source coordinate */

            /* group coordinates */
            tr.gx = (signed int)round(xr * 1000.0);
            tr.gy = (signed int)round(yr * 1000.0);
            tr.ns = (unsigned short)ns;          /* number of samples in this trace */
            tr.dt = (unsigned short)round(dtms); /* sample interval in micro-seconds */
            tr.d1 = (float)tr.dt * 1.0e-6;       /* sample spacing for non-seismic data */

            data.getRow(tempdata, tracl1);
            scai::lama::Scalar tempScalar;
            for (IndexType sample = 0; sample < tempdata.size(); sample++) {
                tempScalar = tempdata.getValue(sample);
                tr.data[sample] = float(tempScalar.getValue<ValueType>());
            }

            fwrite(&tr, 240, 1, pFile);
            fwrite(&tr.data[1], 4, ns, pFile);
        }
        fclose(pFile);
    }
}

/*! \brief Overloading *= Operation
 *
 \param rhs Scalar factor with which the vectors are multiplied.
 */
template <typename ValueType>
KITGPI::Acquisition::Seismogram<ValueType> KITGPI::Acquisition::Seismogram<ValueType>::operator*=(scai::lama::Scalar const &rhs)
{
    data *= rhs;

    return *this;
}

/*! \brief Overloading + Operation
 *
 \param rhs Model which is added.
 */
template <typename ValueType>
KITGPI::Acquisition::Seismogram<ValueType> KITGPI::Acquisition::Seismogram<ValueType>::operator+(KITGPI::Acquisition::Seismogram<ValueType> const &rhs) const
{
    KITGPI::Acquisition::Seismogram<ValueType> result(*this);
    result += rhs;

    return result;
}

/*! \brief Overloading += Operation
 *
 \param rhs Model which is added.
 */
template <typename ValueType>
KITGPI::Acquisition::Seismogram<ValueType> KITGPI::Acquisition::Seismogram<ValueType>::operator+=(KITGPI::Acquisition::Seismogram<ValueType> const &rhs)
{
    data += rhs.data;

    return *this;
}

/*! \brief Overloading - Operation
 *
 \param rhs Model which is subtractet.
 */
template <typename ValueType>
KITGPI::Acquisition::Seismogram<ValueType> KITGPI::Acquisition::Seismogram<ValueType>::operator-(KITGPI::Acquisition::Seismogram<ValueType> const &rhs) const
{
    KITGPI::Acquisition::Seismogram<ValueType> result(*this);
    result -= rhs;

    return result;
}

/*! \brief Overloading -= Operation
 *
 \param rhs Model which is subtractet.
 */
template <typename ValueType>
KITGPI::Acquisition::Seismogram<ValueType> KITGPI::Acquisition::Seismogram<ValueType>::operator-=(KITGPI::Acquisition::Seismogram<ValueType> const &rhs)
{
    data -= rhs.data;

    return *this;
}

/*! \brief Overloading copy assignment operation
 *
 \param rhs Model which is copied.
 */
template <typename ValueType>
KITGPI::Acquisition::Seismogram<ValueType> &KITGPI::Acquisition::Seismogram<ValueType>::operator=(KITGPI::Acquisition::Seismogram<ValueType> const &rhs)
{
    //copy rhs with copy constructor to tmp seismogram
    KITGPI::Acquisition::Seismogram<ValueType> tmp(rhs);
    //swap tmp with *this (lhs)
    swap(tmp);

    return *this;
}

template class KITGPI::Acquisition::Seismogram<double>;
template class KITGPI::Acquisition::Seismogram<float>;
