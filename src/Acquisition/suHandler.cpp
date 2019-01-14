#include "suHandler.hpp"

using namespace scai;

//! \brief Build a Source Acquisition Matrix from all SU files
/*!
 \param filename Name of the file
 \param DH grid spacing
 */
template <typename ValueType>
void KITGPI::Acquisition::suHandler<ValueType>::buildAcqMatrixSource(std::string const &filename, ValueType DH)
{
    buildAcqMatrix(filename, DH, &suHandler<ValueType>::buildAcqMatrixSourceComp);
}

//! \brief Build a Receiver Acquisition Matrix from all SU files
/*!
 \param filename Name of the file
 \param DH grid spacing
 */
template <typename ValueType>
void KITGPI::Acquisition::suHandler<ValueType>::buildAcqMatrixReceiver(std::string const &filename, ValueType DH)
{
    buildAcqMatrix(filename, DH, &suHandler<ValueType>::buildAcqMatrixReceiverComp);
}

//! \brief Build a Acquisition Matrix from all SU files
/*!
 \param filename Name of the file
 \param DH grid spacing
 \param buildAcqMat Function pointer to decide if a source or receiver matrix should be build
 */
template <typename ValueType>
void KITGPI::Acquisition::suHandler<ValueType>::buildAcqMatrix(std::string const &filename, ValueType DH, buildAcqMatPtr buildAcqMat)
{
    acqMat.clear();
    std::vector<lama::DenseMatrix<ValueType>> acqMatVec(NUM_ELEMENTS_SEISMOGRAMTYPE);

    std::string filenameTmp;

    // read all source files
    for (scai::IndexType iComponent = 0; iComponent < NUM_ELEMENTS_SEISMOGRAMTYPE; iComponent++) {
        filenameTmp = filename + "." + std::string(SeismogramTypeString[SeismogramType(iComponent)]) + ".SU";
        (this->*buildAcqMat)(filenameTmp, acqMatVec[iComponent], DH);
        nShots[iComponent] = acqMatVec[iComponent].getNumRows();
    }

    // count number of sources to allocate final acquisition matrix
    IndexType numSources = 0;
    for (auto it = acqMatVec.begin(); it != acqMatVec.end(); it++)
        numSources += (*it).getNumRows();

    auto rowDist = std::make_shared<scai::dmemo::SingleDistribution>(numSources, dmemo::Communicator::getCommunicatorPtr(), 0);
    auto colDist = std::make_shared<scai::dmemo::NoDistribution>(acqMatVec[0].getNumColumns());

    acqMat.allocate(rowDist, colDist);

    // concatenate the acquisition matrices of all components
    auto write_acquisition_HA = hostWriteAccess(acqMat.getLocalStorage().getData());

    IndexType j = 0;
    for (auto it = acqMatVec.begin(); it != acqMatVec.end(); it++) {
        auto read_acquisition_HA = hostReadAccess((*it).getLocalStorage().getData());
        for (IndexType i = 0; i < read_acquisition_HA.size(); i++) {
            write_acquisition_HA[j] = read_acquisition_HA[i];
            j++;
        }
        read_acquisition_HA.release();
    }
    write_acquisition_HA.release();
}

//! \brief Build a source Acquisition Matrix from one SU file
/*!
 \param filename Name of the file
 \param acqMat Acquisition Matrix
 \param DH grid spacing
 */
template <typename ValueType>
void KITGPI::Acquisition::suHandler<ValueType>::buildAcqMatrixSourceComp(std::string const &filename, lama::DenseMatrix<ValueType> &acqMatTmp, ValueType DH)
{
    std::vector<Segy> header;
    readHeaderSU(filename, header);
    Segy thisHeader;

    IndexType component = getComponentFromName(filename);

    auto rowDist = std::make_shared<scai::dmemo::SingleDistribution>(header.size(), dmemo::Communicator::getCommunicatorPtr(), 0);
    auto colDist = std::make_shared<scai::dmemo::NoDistribution>(9);

    acqMatTmp.allocate(rowDist, colDist);

    auto write_acquisition_HA = writeAccess(acqMatTmp.getLocalStorage().getData());

    for (int i = 0; i < write_acquisition_HA.size() / 9; i++) {
        thisHeader = header[i];
        write_acquisition_HA[i * 9 + 0] = common::Math::floor<ValueType>(thisHeader.sx * common::Math::pow<ValueType>(10, thisHeader.scalco) / DH + 0.5);
        write_acquisition_HA[i * 9 + 1] = common::Math::floor<ValueType>(thisHeader.sdepth * common::Math::pow<ValueType>(10, thisHeader.scalel) / DH + 0.5);
        write_acquisition_HA[i * 9 + 2] = common::Math::floor<ValueType>(thisHeader.sy * common::Math::pow<ValueType>(10, thisHeader.scalco) / DH + 0.5);
        write_acquisition_HA[i * 9 + 3] = ValueType(component);
        write_acquisition_HA[i * 9 + 4] = 3; // each source signal should be read from file
    }

    write_acquisition_HA.release();
}

//! \brief Build a source Acquisition Matrix from one SU file
/*!
 \param filename Name of the file
 \param acqMat Acquisition Matrix
 \param DH grid spacing
 */
template <typename ValueType>
void KITGPI::Acquisition::suHandler<ValueType>::buildAcqMatrixReceiverComp(std::string const &filename, lama::DenseMatrix<ValueType> &acqMatTmp, ValueType DH)
{
    std::vector<Segy> header;
    readHeaderSU(filename, header);
    Segy thisHeader;

    IndexType component = getComponentFromName(filename);

    auto rowDist = std::make_shared<scai::dmemo::SingleDistribution>(header.size(), dmemo::Communicator::getCommunicatorPtr(), 0);
    auto colDist = std::make_shared<scai::dmemo::NoDistribution>(4);

    acqMatTmp.allocate(rowDist, colDist);

    auto write_acquisition_HA = writeAccess(acqMatTmp.getLocalStorage().getData());

    for (int i = 0; i < write_acquisition_HA.size() / 4; i++) {
        thisHeader = header[i];
        write_acquisition_HA[i * 4 + 0] = common::Math::floor<ValueType>(thisHeader.gx * common::Math::pow<ValueType>(10, thisHeader.scalco) / DH);
        write_acquisition_HA[i * 4 + 1] = common::Math::floor<ValueType>(thisHeader.gelev * common::Math::pow<ValueType>(10, thisHeader.scalel) / DH);
        write_acquisition_HA[i * 4 + 2] = common::Math::floor<ValueType>(thisHeader.gy * common::Math::pow<ValueType>(10, thisHeader.scalco) / DH);
        write_acquisition_HA[i * 4 + 3] = ValueType(component);
    }

    write_acquisition_HA.release();
}

//! \brief Locate the trace in the SU files based on a shot number
/*!
 \param filename Filename which gets the correct suffix
 \param traceNumber Trace number the shot has in file filename
 \param shotNumber Shot which should be located
 */
template <typename ValueType>
void KITGPI::Acquisition::suHandler<ValueType>::locateTrace(std::string &filename, scai::IndexType &traceNumber, scai::IndexType shotNumber)
{
    IndexType iComponent = 0;
    do {
        shotNumber -= nShots[iComponent];
        iComponent++;
    } while (shotNumber >= 0);

    filename += "." + std::string(SeismogramTypeString[SeismogramType(iComponent - 1)]) + ".SU";
    traceNumber = nShots[iComponent - 1] + shotNumber;
}

//! \brief Derive the component of the SU file based on its filename
/*!
 \param filename Filename where the component is the suffix
 */
template <typename ValueType>
scai::IndexType KITGPI::Acquisition::suHandler<ValueType>::getComponentFromName(std::string const &filename)
{

    IndexType iTmp = filename.find_last_of('.');
    std::string tmpString = filename.substr(0, iTmp);
    iTmp = tmpString.find_last_of('.');
    tmpString = tmpString.substr(iTmp + 1);

    IndexType component = NUM_ELEMENTS_SEISMOGRAMTYPE + 2;

    for (IndexType i = 0; i < NUM_ELEMENTS_SEISMOGRAMTYPE; i++)
        if (tmpString.compare(SeismogramTypeString[SeismogramType(i)]) == 0)
            component = i + 1; // +1 because components in acquisition matrix are defined as p := 1 ...

    SCAI_ASSERT(component < NUM_ELEMENTS_SEISMOGRAMTYPE + 2, "no component found in filename")

    return component;
}

//! \brief Getter member function for the Acquisition Matrix
template <typename ValueType>
lama::DenseMatrix<ValueType> const &KITGPI::Acquisition::suHandler<ValueType>::getAcquisition() const
{
    return (acqMat);
}

//! \brief Getter member function for the Acquisition Matrix
/*!
 \param acqRowMat Acquisition Matrix
 \param shotNumber Shot number to get
 */
template <typename ValueType>
void KITGPI::Acquisition::suHandler<ValueType>::getAcquisitionRow(lama::DenseMatrix<ValueType> &acqRowMat, IndexType shotNumber) const
{
    scai::lama::DenseVector<ValueType> acqRow;
    acqMat.getRow(acqRow, shotNumber);

    auto rowDist = std::make_shared<scai::dmemo::SingleDistribution>(1, dmemo::Communicator::getCommunicatorPtr(), 0);
    auto colDist = std::make_shared<scai::dmemo::NoDistribution>(9);

    acqRowMat.allocate(rowDist, colDist);
    acqRowMat.setRow(acqRow, 0, common::BinaryOp::COPY);
}

//! \brief Initialize a Segy struct
/*!
\param segy Segy struct
*/
template <typename ValueType>
void KITGPI::Acquisition::suHandler<ValueType>::initSegy(Segy &tr)
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
void KITGPI::Acquisition::suHandler<ValueType>::readDataSU(std::string const &filename, scai::lama::Matrix<ValueType> &data, scai::IndexType ns, scai::IndexType ntr)
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
void KITGPI::Acquisition::suHandler<ValueType>::readSingleDataSU(std::string const &filename, scai::lama::Vector<ValueType> &data, scai::IndexType traceNumber)
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
\param data DenseMatrix with traces of one seismogramtype
\param coordinates1D coordinates of the traces
\param sourceCoordinate1D source coordinate (is only !=0 if a single source is used)
\param modelCoordinates Coordinate class, which eg. maps 3D coordinates to 1D model indices
*/
template <typename ValueType>
void KITGPI::Acquisition::suHandler<ValueType>::writeSU(std::string const &filename, scai::lama::DenseMatrix<ValueType> const &data, scai::lama::DenseVector<scai::IndexType> const &coordinates1D, ValueType DT, scai::IndexType sourceCoordinate1D, Coordinates<ValueType> const &modelCoordinates)
{
    Segy tr;
    initSegy(tr);

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

    coordinate3D coord3Dsrc;
    coordinate3D coord3Drec;
    coord3Dsrc = modelCoordinates.index2coordinate(sourceCoordinate1D);
    ValueType DH = modelCoordinates.getDH();

    YS = coord3Dsrc.y;
    XS = coord3Dsrc.x;
    ZS = coord3Dsrc.z;
    YS = YS * DH;
    XS = XS * DH;
    ZS = ZS * DH;

    for (tracl1 = 0; tracl1 < ntr; tracl1++) {
        temp3 = float(coordinates1D.getValue(tracl1));
        temp2 = floor(temp3);
        coord3Drec = modelCoordinates.index2coordinate(temp2);
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
void KITGPI::Acquisition::suHandler<ValueType>::readHeaderSU(std::string const &filename, std::vector<Segy> &header)
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

template class KITGPI::Acquisition::suHandler<double>;
template class KITGPI::Acquisition::suHandler<float>;
