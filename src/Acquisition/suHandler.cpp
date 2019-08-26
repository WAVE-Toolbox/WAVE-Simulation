#include "suHandler.hpp"
#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/dmemo/CollectiveFile.hpp>
using namespace scai;

//! \brief Build a Source Acquisition Matrix from all SU files
/*!
 \param filename Name of the file
 \param DH grid spacing
 */
template <typename ValueType>
void KITGPI::Acquisition::suHandler<ValueType>::buildAcqMatrixSource(std::string const &filename, ValueType DH)
{
    acqSource.clear();
    std::vector<sourceSettings<ValueType>> sourceSettingsVecTmp;

    std::string filenameTmp;

    // read all source files
    for (scai::IndexType iComponent = 0; iComponent < NUM_ELEMENTS_SEISMOGRAMTYPE; iComponent++) {
        filenameTmp = filename + "." + std::string(SeismogramTypeString[SeismogramType(iComponent)]) + ".SU";
        buildAcqMatrixSourceComp(filenameTmp, sourceSettingsVecTmp, DH);
        nShots[iComponent] = sourceSettingsVecTmp.size();
        acqSource.insert(acqSource.end(), sourceSettingsVecTmp.begin(), sourceSettingsVecTmp.end());
    }
}

//! \brief Build a Receiver Acquisition Matrix from all SU files
/*!
 \param filename Name of the file
 \param DH grid spacing
 */
template <typename ValueType>
void KITGPI::Acquisition::suHandler<ValueType>::buildAcqMatrixReceiver(std::string const &filename, ValueType DH)
{
    acqReceiver.clear();
    std::vector<receiverSettings> receiverSettingsVecTmp;

    std::string filenameTmp;

    // read all receiver files
    for (scai::IndexType iComponent = 0; iComponent < NUM_ELEMENTS_SEISMOGRAMTYPE; iComponent++) {
        filenameTmp = filename + "." + std::string(SeismogramTypeString[SeismogramType(iComponent)]) + ".SU";
        buildAcqMatrixReceiverComp(filenameTmp, receiverSettingsVecTmp, DH);
        nShots[iComponent] = receiverSettingsVecTmp.size();
        acqReceiver.insert(acqReceiver.end(), receiverSettingsVecTmp.begin(), receiverSettingsVecTmp.end());
    }
}

//! \brief Build a source Acquisition Matrix from one SU file
/*!
 \param filename Name of the file
 \param sourceSettingsVec Vector of source settings
 \param DH grid spacing
 */
template <typename ValueType>
void KITGPI::Acquisition::suHandler<ValueType>::buildAcqMatrixSourceComp(std::string const &filename, std::vector<sourceSettings<ValueType>> &sourceSettingsVec, ValueType DH)
{
    std::vector<Segy> header;
    readHeaderSU(filename, header);
    Segy thisHeader;

    IndexType component = getComponentFromName(filename);
    sourceSettingsVec.clear();
    sourceSettingsVec.reserve(header.size());

    for (unsigned int i = 0; i < header.size(); i++) {
        thisHeader = header[i];
        sourceSettingsVec[i].sourceCoords.x = static_cast<IndexType>(thisHeader.sx * common::Math::pow<ValueType>(10, thisHeader.scalco) / DH + 0.5);
        sourceSettingsVec[i].sourceCoords.y = static_cast<IndexType>(thisHeader.sdepth * common::Math::pow<ValueType>(10, thisHeader.scalel) / DH + 0.5);
        sourceSettingsVec[i].sourceCoords.z = static_cast<IndexType>(thisHeader.sy * common::Math::pow<ValueType>(10, thisHeader.scalco) / DH + 0.5);
        sourceSettingsVec[i].sourceType = component;
        sourceSettingsVec[i].waveletType = 3; // each source signal should be read from file
    }
}

//! \brief Build a source Acquisition Matrix from one SU file
/*!
 \param filename Name of the file
 \param acqMat Acquisition Matrix
 \param DH grid spacing
 */
template <typename ValueType>
void KITGPI::Acquisition::suHandler<ValueType>::buildAcqMatrixReceiverComp(std::string const &filename, std::vector<receiverSettings> &receiverSettingsVec, ValueType DH)
{
    std::vector<Segy> header;
    readHeaderSU(filename, header);
    Segy thisHeader;

    IndexType component = getComponentFromName(filename);
    receiverSettingsVec.clear();
    receiverSettingsVec.reserve(header.size());

    for (unsigned int i = 0; i < header.size() / 4; i++) {
        thisHeader = header[i];
        receiverSettingsVec[i].receiverCoords.x = static_cast<IndexType>(thisHeader.gx * common::Math::pow<ValueType>(10, thisHeader.scalco) / DH);
        receiverSettingsVec[i].receiverCoords.y = static_cast<IndexType>(thisHeader.gelev * common::Math::pow<ValueType>(10, thisHeader.scalel) / DH);
        receiverSettingsVec[i].receiverCoords.z = static_cast<IndexType>(thisHeader.gy * common::Math::pow<ValueType>(10, thisHeader.scalco) / DH);
        receiverSettingsVec[i].receiverType = component;
    }
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

template <typename ValueType>
std::vector<KITGPI::Acquisition::sourceSettings<ValueType>> &KITGPI::Acquisition::suHandler<ValueType>::getSourceSettingsVec()
{
    return (acqSource);
}

template <typename ValueType>
std::vector<KITGPI::Acquisition::receiverSettings> &KITGPI::Acquisition::suHandler<ValueType>::getReceiverSettingsVec()
{
    return (acqReceiver);
}

template <typename ValueType>
KITGPI::Acquisition::sourceSettings<ValueType> &KITGPI::Acquisition::suHandler<ValueType>::getSourceSettings(IndexType shotNumber)
{
    return (acqSource[shotNumber]);
}

//! \brief Initialize a Segy struct
/*!
\param tr Segy struct
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
void KITGPI::Acquisition::suHandler<ValueType>::readDataSU(std::string const &filename, scai::lama::DenseMatrix<ValueType> &data, scai::IndexType ns, scai::IndexType ntr)
{
    Segy tr;

    // write su pararllel
    // 1 redistribute data matrix and coordinate vector to block distribution
    auto colDistIn = data.getColDistributionPtr();
    auto comm = data.getRowDistributionPtr()->getCommunicatorPtr();
    auto rowDist = std::make_shared<dmemo::BlockDistribution>(ntr, comm);
    
    
    data.redistribute(rowDist, colDistIn);

    auto numLocalTraces = data.getLocalNumRows();
    scai::hmemo::HArray<char> localBuffer(numLocalTraces * (240 + sizeof(float) * ns));
    
    auto cfile = comm->collectiveFile();
    const char *filetemp = filename.c_str();
    cfile->open(filetemp, "r");
    cfile->readAll(localBuffer,numLocalTraces* (240 + sizeof(float) * ns));
    cfile->close();
    
    scai::hmemo::HArray<float> localTraces(numLocalTraces*ns);
    auto writeLocalData = hmemo::hostWriteAccess(localTraces);
    float *writePointer = writeLocalData.get();

    
    // create local (byte) Harray with type char (1 char = 1 byte)
    // size of array = numLocalTraces*(HeaderSize+Tracessize)
    // get read access and pointer to the buffer
    auto readLocalBuffer = hmemo::hostReadAccess(localBuffer);
    const char *readPointer = readLocalBuffer.get();
    
    for (scai::IndexType localTrace = 0; localTrace < numLocalTraces; localTrace++) {
        readPointer += 240;
        std::memcpy((void *)writePointer, (void *)readPointer, sizeof(float)*ns);
        writePointer += ns;
        readPointer += sizeof(float)*ns; // is float pointer
    }
writeLocalData.release();

data.getLocalStorage().setRawDenseData(numLocalTraces,ns,hmemo::hostReadAccess(localTraces).get());

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
    auto ns = data.getNumColumns();
    auto ntr = data.getNumRows();

    // write su pararllel
    // 1 redistribute data matrix and coordinate vector to block distribution
    auto colDist = data.getColDistributionPtr();
    auto rowDistIn = data.getRowDistributionPtr();
    auto comm = rowDistIn->getCommunicatorPtr();
    auto rowDist = std::make_shared<dmemo::BlockDistribution>(ntr, comm);
    lama::DenseMatrix<float> dataTemp;
    lama::DenseVector<IndexType> coordinatesTemp;
    dataTemp.assignDistribute(data, rowDist, colDist);
    coordinatesTemp.assignDistribute(coordinates1D, rowDist);
    auto readLocalCoordinates = hmemo::hostReadAccess(coordinatesTemp.getLocalValues());

    // 1 get number of local traces
    auto numLocalTraces = dataTemp.getLocalNumRows();
    //read access

    //get read access and pointer to the data
    auto readLocalData = hmemo::hostReadAccess(dataTemp.getLocalStorage().getValues());
    const float *readPointer = readLocalData.get();

    // create local (byte) Harray with type char (1 char = 1 byte)
    // size of array = numLocalTraces*(HeaderSize+Tracessize)
    scai::hmemo::HArray<char> localBuffer(numLocalTraces * (240 + sizeof(float) * ns));
    // get write access and pointer to the buffer
    auto writeLocalBuffer = hmemo::hostWriteAccess(localBuffer);
    char *writePointer = writeLocalBuffer.get();

    // create Segy header
    Segy tr;
    initSegy(tr);

    ValueType xr, yr, zr, x, y, z;
    ValueType XS = 0.0, YS = 0.0, ZS = 0.0;
    const ValueType xshift = 800.0, yshift = 800.0;
    ValueType dtms = (ValueType)(DT * 1000000);

    tr.ntr = ntr; /* number of traces */

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

    // loop over local traces
    for (int localTrace = 0; localTrace < numLocalTraces; localTrace++) {

        IndexType coordinateIndex = readLocalCoordinates[localTrace];
        coord3Drec = modelCoordinates.index2coordinate(coordinateIndex);
        xr = (ValueType)coord3Drec.x;
        yr = (ValueType)coord3Drec.y;
        zr = (ValueType)coord3Drec.z;
        yr = yr * DH;
        xr = xr * DH;
        zr = zr * DH;
        x = xr - XS; // Taking source position as reference point
        y = yr - YS;
        z = zr - ZS;

        tr.tracl = (int)rowDist->local2Global(localTrace) + 1; // trace sequence number within line
        tr.tracr = 1;                                          // trace sequence number within reel
        tr.ep = 1;
        tr.cdp = (int)ntr;
        tr.trid = (short)1;
        tr.offset = (signed int)round(sqrt((XS - xr) * (XS - xr) + (YS - yr) * (YS - yr) + (ZS - zr) * (ZS - zr)) * 1000.0);
        tr.gelev = (signed int)round(yr * 1000.0);
        tr.sdepth = (signed int)round(YS * 1000.0); /* source depth (positive) */
        /* angle between receiver position and reference point
            (sperical coordinate system: swdep=theta, gwdep=phi) */
        tr.gdel = (signed int)round(atan2(-y, z) * 180 * 1000.0 / 3.1415926);
        tr.gwdep = (signed int)round(sqrt(z * z + y * y) * 1000.0);
        tr.swdep = (int)round(((360.0 / (2.0 * 3.1415926)) * atan2(x - xshift, y - yshift)) * 1000.0);
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

        std::memcpy((void *)writePointer, &tr, 240);
        writePointer += 240;
        std::memcpy((void *)writePointer, (void *)readPointer, sizeof(float) * ns);
        writePointer += sizeof(float) * ns;
        readPointer += ns; // is float pointer
    }
    writeLocalBuffer.release();

    //get collective file to write with parall IO
    auto cfile = comm->collectiveFile();
    const char *filetemp = filename.c_str();
    cfile->open(filetemp, "w");
    cfile->writeAll(localBuffer);
    cfile->close();
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
