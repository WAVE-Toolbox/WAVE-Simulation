
#pragma once

#include <scai/lama.hpp>
#include <scai/dmemo/CyclicDistribution.hpp>
#include <scai/lama/DenseVector.hpp>

#include "Sources.hpp"
#include "Receivers.hpp"
#include "Coordinates.hpp"
#include "segy.h"
#include "../Configuration/Configuration.hpp"
#include "math.h"

namespace KITGPI {
    
    namespace Acquisition {
        
        //! Seismogram class
        /*!
         * This class handels a seismogram which consists of several traces.
         */
        template <typename ValueType>
        class Seismogram
        {
            
        public:
            
            //! Default constructor
            Seismogram():numSamples(0),numTracesGlobal(0),numTracesLocal(0),DT(0.0){};
            
            //! Default destructor
            ~Seismogram(){};
            
            void writeToFileRaw(std::string filename); // write only the raw data in text format
            void writeToFileSU(std::string filename, IndexType NX, IndexType NY, IndexType NZ, ValueType DH); // write to Seismic Unix format
            void write(Configuration::Configuration<ValueType>& config);
            
            void readFromFileRaw(std::string filename,dmemo::DistributionPtr distTraces,dmemo::DistributionPtr distSamples);
            
            void allocate(hmemo::ContextPtr ctx, dmemo::DistributionPtr distSeismogram, IndexType NT);
            
            void init(Receivers<ValueType>& receiver, IndexType NT, hmemo::ContextPtr ctx);
            
            void redistribute(dmemo::DistributionPtr distRow,dmemo::DistributionPtr distColumn=NULL);
            void replicate();
            
            void resetData();
            
            IndexType getNumTracesGlobal();
            IndexType getNumTracesLocal();
            IndexType getNumSamples();
            ValueType getDT();
            lama::DenseMatrix<ValueType>& getData();
            lama::DenseVector<ValueType>& getTraceType();
            lama::DenseVector<ValueType>& getCoordinates();
            
            void setDT(ValueType newDT);
            void setSourceCoordinate(IndexType sourceCoord);
            
        private:
            
            IndexType numSamples; //!< Number of samples of one trace
            IndexType numTracesGlobal; //!< Number of global traces
            IndexType numTracesLocal; //!< Number of local traces
            
            /* header information */
            ValueType DT; //!< Temporal sampling in seconds
            lama::DenseVector<ValueType> traceType; //!< Type of trace (receiver or source type)
            lama::DenseVector<ValueType> coordinates; //!< Coordinates of the traces
            IndexType sourceCoordinate; //!< Coordinate of source
            
            /* raw data */
            lama::DenseMatrix<ValueType> data; //!< Raw seismogram data
            
        };
    }
}

//! \brief Write seismogram to disk
/*!
 *
 \param config Configuration class
 */
template <typename ValueType>
void KITGPI::Acquisition::Seismogram<ValueType>::write(Configuration::Configuration<ValueType>& config){
    
    switch(config.getSeismogramFormat()){
        case 1:
            writeToFileRaw(config.getSeismogramFilename());
            break;
        case 2:
            setDT(config.getDT());
            writeToFileSU(config.getSeismogramFilename(), config.getNX(), config.getNY(), config.getNZ(), config.getDH());
            break;
        default:
            COMMON_THROWEXCEPTION(" Unkown SeismogramFormat " )
            break;
    }
    
}


//! \brief Set the temporal sampling DT
/*!
 *
 \param newDT Temporal sampling which will be set to the seismogram
 */
template <typename ValueType>
void KITGPI::Acquisition::Seismogram<ValueType>::setDT(ValueType newDT){
    if(newDT<0){
        COMMON_THROWEXCEPTION(" DT is smaller zero. ")
    }
    DT=newDT;
}

//! \brief Set the source coordinate
/*!
 *
 \param sourceCoord Source coordinate in 1-D format
 */
template <typename ValueType>
void KITGPI::Acquisition::Seismogram<ValueType>::setSourceCoordinate(IndexType sourceCoord){
    sourceCoordinate=sourceCoord;
}


//! \brief Initiate the seismogram by receivers
/*!
 *
 \param receiver Receivers which will record into this seismogram
 \param NT Total number of time steps
 \param ctx Context
 */
template <typename ValueType>
void KITGPI::Acquisition::Seismogram<ValueType>::init(Receivers<ValueType>& receiver, IndexType NT, hmemo::ContextPtr ctx){
    
    /* Allocation */
    dmemo::DistributionPtr dist_traces=receiver.getReceiversDistribution();
    allocate(ctx,dist_traces,NT);
    
    /* set header information */
    lama::DenseVector<ValueType>& coordinates_temp=receiver.getCoordinates();
    lama::DenseVector<ValueType>& traceType_temp=receiver.getReceiversType();
    coordinates=coordinates_temp;
    traceType=traceType_temp;
    
    numTracesLocal=receiver.getNumReceiversLocal();
    numTracesGlobal=receiver.getNumReceiversGlobal();
    numSamples=NT;
}


//! \brief Get reference to Receiver Type
/*!
 *
 * THIS METHOD IS CALLED DURING TIME STEPPING
 * DO NOT WASTE RUNTIME HERE
 *
 */
template <typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Acquisition::Seismogram<ValueType>::getTraceType(){
    return(traceType);
}


//! \brief Get reference to coordinates
/*!
 *
 * THIS METHOD IS CALLED DURING TIME STEPPING
 * DO NOT WASTE RUNTIME HERE
 *
 */
template <typename ValueType>
lama::DenseVector<ValueType>& KITGPI::Acquisition::Seismogram<ValueType>::getCoordinates(){
    return(coordinates);
}


//! \brief Get reference to seismogram data
/*!
 *
 * THIS METHOD IS CALLED DURING TIME STEPPING
 * DO NOT WASTE RUNTIME HERE
 *
 */
template <typename ValueType>
lama::DenseMatrix<ValueType>& KITGPI::Acquisition::Seismogram<ValueType>::getData(){
    return(data);
}


//! \brief Replicate seismogram on all processes
/*!
 * Creates a copy of the seismogram on all processe
 */
template <typename ValueType>
void KITGPI::Acquisition::Seismogram<ValueType>::replicate()
{
    dmemo::DistributionPtr no_dist_numSamples( new scai::dmemo::NoDistribution ( numSamples ) );
    
    dmemo::DistributionPtr no_dist_Traces( new scai::dmemo::NoDistribution ( numTracesGlobal ) );
    
    redistribute(no_dist_Traces,no_dist_numSamples);
}


//! \brief Allocate seismogram
/*!
 * Allocates seismogram based on a given distribution of the traces and the number of samples per trace
 *
 \param ctx Context
 \param distTraces Distribution for traces
 \param NT Total number of samples per trace
 */
template <typename ValueType>
void KITGPI::Acquisition::Seismogram<ValueType>::allocate(hmemo::ContextPtr ctx, dmemo::DistributionPtr distTraces, IndexType NT)
{
    
    numSamples=NT;
    numTracesGlobal=distTraces->getGlobalSize();
    numTracesLocal=distTraces->getLocalSize();
    
    dmemo::DistributionPtr no_dist_NT( new scai::dmemo::NoDistribution ( NT ) );
    
    data.setContextPtr(ctx);
    traceType.setContextPtr(ctx);
    coordinates.setContextPtr(ctx);
    
    data.allocate(distTraces,no_dist_NT);
    traceType.allocate(distTraces);
    coordinates.allocate(distTraces);
}


//! \brief reset seismogram set the seismogram data to zero
/*!
 */
template <typename ValueType>
void KITGPI::Acquisition::Seismogram<ValueType>::resetData()
{
    data.scale(0.0);
}


//! \brief Redistribute seismogram data
/*!
 *
 \param distTraces Distribution of traces
 \param distSamples Distribution of temporal samples
 */
template <typename ValueType>
void KITGPI::Acquisition::Seismogram<ValueType>::redistribute(dmemo::DistributionPtr distTraces,dmemo::DistributionPtr distSamples)
{
    if(distSamples==NULL){
        dmemo::DistributionPtr distSamples( new scai::dmemo::NoDistribution ( numSamples ) );
    }
    
    data.redistribute(distTraces,distSamples);
    traceType.redistribute(distTraces);
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
void KITGPI::Acquisition::Seismogram<ValueType>::readFromFileRaw(std::string filename,dmemo::DistributionPtr distTraces,dmemo::DistributionPtr distSamples)
{
    data.ReadFromFile(filename);
    
    IndexType nrow_temp=data.getNumRows();
    IndexType ncolumn_temp=data.getNumColumns();
    
    numSamples=ncolumn_temp;
    numTracesGlobal=nrow_temp;
    
    if(distTraces==NULL && distSamples==NULL){
        replicate();
    } else {
        redistribute(distTraces,distSamples);
    }
    
}


//! \brief Write a seismogram to disk without header
/*!
 *
 \param filename Filename to write seismogram
 */
template <typename ValueType>
void KITGPI::Acquisition::Seismogram<ValueType>::writeToFileRaw(std::string filename)
{
    if(data.getNumValues()==0) {
        COMMON_THROWEXCEPTION("Seismogramm data is not allocated")
    }
    data.writeToFile(filename);
}


//! \brief Get temporal sampling
template <typename ValueType>
ValueType KITGPI::Acquisition::Seismogram<ValueType>::getDT()
{
    if(DT==0.0){
        COMMON_THROWEXCEPTION("DT is not set for this seismogram")
    }
    return(DT);
}


//! \brief Get number of samples per trace
template <typename ValueType>
IndexType KITGPI::Acquisition::Seismogram<ValueType>::getNumSamples()
{
    if(numSamples==0){
        numSamples=data.getNumColumns();
        if(numSamples==0){
            COMMON_THROWEXCEPTION("Seismogram is not allocated")
        }
    }
    return(numSamples);
}


//! \brief Get number of local traces
template <typename ValueType>
IndexType KITGPI::Acquisition::Seismogram<ValueType>::getNumTracesLocal()
{
    return(numTracesLocal);
}


//! \brief Get number of global traces
template <typename ValueType>
IndexType KITGPI::Acquisition::Seismogram<ValueType>::getNumTracesGlobal()
{
    if(numTracesGlobal==0){
        numTracesGlobal=data.getNumRows();
        if(numTracesGlobal==0){
            COMMON_THROWEXCEPTION("Seismogram is not allocated")
        }
    }
    return(numTracesGlobal);
}

//! \brief Write a seismogram to disk in Seismic Unix (SEG-Y) format
/*!
 *
 \param filename Filename to write seismogram in Seismic Unix (SEG-Y) format
 \param NX Number of grid points in X direction
 \param NY Number of grid points in Y direction
 \param NZ Number of grid points in Z direction
 \param DH Length of space step in meter
 */
template <typename ValueType>
void KITGPI::Acquisition::Seismogram<ValueType>::writeToFileSU(std::string filename,IndexType NX, IndexType NY, IndexType NZ, ValueType DH)
{
    Segy tr;
    /* Define parameters in tr header */
    tr.tracr=0      ;       /* trace sequence number within reel */
    tr.fldr=0       ;       /* field record number */
    tr.tracf=0      ;       /* trace number within field record */
    tr.ep=0         ;       /* energy source point number */
    tr.cdpt=0       ;
    
    /* trace number within CDP ensemble */
    tr.nvs=0        ;   /* number of vertically summed traces (see vscode
                         in bhed structure) */
    tr.nhs=0        ;   /* number of horizontally summed traces (see vscode
                         in bhed structure) */
    tr.duse=0       ;   /* data use:
                         1 = production
                         2 = test */
    tr.selev=0      ; /* source elevation from sea level
                       (above sea level is positive) */
    tr.gdel=0       ; /* datum elevation at receiver group */
    tr.sdel=0       ; /* datum elevation at source */
    tr.gwdep=0      ; /* water depth at receiver group */
    tr.sy=0         ;   /* Y source coordinate */
    tr.gy=0         ;   /* Y group coordinate */
    tr.counit=1     ;   /* coordinate units code:
                         for previous four entries
                         1 = length (meters or feet)
                         2 = seconds of arc (in this case, the
                         X values are longitude and the Y values
                         are latitude, a positive value designates
                         the number of seconds east of Greenwich
                         or north of the equator */
    tr.wevel=0     ;        /* weathering velocity */
    tr.swevel=0    ;        /* subweathering velocity */
    tr.sut=0       ;        /* uphole time at source */
    tr.gut=0       ;        /* uphole time at receiver group */;
    
    tr.sstat=0     ;        /* source static correction */
    tr.gstat=0     ;        /* group static correction */
    tr.tstat=0     ;        /* total static applied */
    tr.laga=0      ; /* lag time A, time in ms between end of 240-
                      byte trace identification header and time
                      break, positive if time break occurs after
                      end of header, time break is defined as
                      the initiation pulse which maybe recorded
                      on an auxiliary trace or as otherwise
                      specified by the recording system */
    tr.lagb=0       ; /* lag time B, time in ms between the time break
                       and the initiation time of the energy source,
                       may be positive or negative */
    tr.delrt=0      ; /* delay recording time, time in ms between
                       initiation time of energy source and time
                       when recording of data samples begins
                       (for deep water work if recording does not
                       start at zero time) */
    tr.muts=0      ; /* mute time--start */
    tr.mute=0      ; /* mute time--end */
    tr.gain=0      ; /* gain type of field instruments code:
                      1 = fixed
                      2 = binary
                      3 = floating point
                      4 ---- N = optional use */
    tr.igc=0       ; /* instrument gain constant */
    tr.igi=0       ; /* instrument early or initial gain */
    tr.corr=0      ; /* correlated:
                      1 = no
                      2 = yes */
    tr.sfs=0       ; /* sweep frequency at start */
    tr.sfe=0       ; /* sweep frequency at end */
    tr.slen=0      ; /* sweep length in ms */
    tr.styp=0      ; /* sweep type code:
                      1 = linear
                      2 = cos-squared
                      3 = other */
    tr.stas=0       ; /* sweep trace length at start in ms */
    tr.stae=0       ; /* sweep trace length at end in ms */
    tr.tatyp=0      ; /* taper type: 1=linear, 2=cos^2, 3=other */
    tr.afilf=0      ; /* alias filter frequency if used */
    tr.afils=0      ; /* alias filter slope */
    tr.nofilf=0     ; /* notch filter frequency if used */
    tr.nofils=0     ; /* notch filter slope */
    tr.lcf=0        ; /* low cut frequency if used */
    tr.hcf=0        ; /* high cut frequncy if used */
    tr.lcs=0        ; /* low cut slope */
    tr.hcs=0        ; /* high cut slope */
    tr.year=0       ; /* year data recorded */
    tr.day=0        ; /* day of year */
    tr.hour=0       ; /* hour of day (24 hour clock) */
    tr.minute=0     ; /* minute of hour */
    tr.sec=0        ; /* second of minute */
    tr.timbas=0     ; /* time basis code:
                       1 = local
                       2 = GMT
                       3 = other */
    tr.trwf=0       ; /* trace weighting factor, defined as 1/2^N
                       volts for the least sigificant bit */
    tr.grnors=0     ; /* geophone group number of roll switch
                       position one */
    tr.grnofr=0     ; /* geophone group number of trace one within
                       original field record */
    tr.grnlof=0     ; /* geophone group number of last trace within
                       original field record */
    tr.gaps=0       ;  /* gap size (total number of groups dropped) */
    tr.otrav=0      ;  /* overtravel taper code:
                        1 = down (or behind)
                        2 = up (or ahead) */
    
    /* local assignments */
    
    tr.f1=0.0;      /* first sample location for non-seismic data */
    
    tr.d2=0.0;      /* sample spacing between traces */
    
    tr.f2=0.0;      /* first trace location */
    
    tr.ungpow=0.0;  /* negative of power used for dynamic
                     range compression */
    tr.unscale=0.0; /* reciprocal of scaling factor to normalize
                     range */
    tr.mark=0     ;
    
    int tracl1;
    lama::Scalar temp;
    double temp3;
    IndexType temp2;
    float xr, yr, zr, x, y, z;
    float XS=0.0, YS=0.0, ZS=0.0;
    const float xshift=800.0, yshift=800.0;
    float dtms=float(DT*1000000);
    
    int ns=int(numSamples);
    int ntr=int(numTracesGlobal);
    tr.ntr=ntr;    /* number of traces */
    std::string tempstr=filename;
    const char*filetemp=tempstr.c_str();
    FILE * pFile;
    pFile=fopen(filetemp,"wb");
    lama::DenseVector<ValueType> tempdata;     
    
    Coordinates<ValueType> coordTransform;
    coordinate3D coord3Dsrc;
    coordinate3D coord3Drec;
    coord3Dsrc=coordTransform.index2coordinate(sourceCoordinate,NX,NY,NZ);  
    YS=coord3Dsrc.y;
    XS=coord3Dsrc.x;
    ZS=coord3Dsrc.z;
    YS=YS*DH;
    XS=XS*DH;        
    ZS=ZS*DH;
    
    for(tracl1=0;tracl1<ntr;tracl1++){
        temp=coordinates.getValue(tracl1);
        temp3=float(temp.getValue<ValueType>());
        temp2=floor(temp3);
        coord3Drec=coordTransform.index2coordinate(temp2,NX,NY,NZ);
        xr=coord3Drec.x;
        yr=coord3Drec.y;
        zr=coord3Drec.z;
        yr=yr*DH;
        xr=xr*DH;
        zr=zr*DH;
        x=xr-XS; // Taking source position as reference point
        y=yr-YS;
        z=zr-ZS;
        
        tr.tracl=tracl1+1;   // trace sequence number within line
        tr.tracr=1;          // trace sequence number within reel
        tr.ep=1;
        tr.cdp=ntr;
        tr.trid=(short)1;  
        tr.offset=(signed int)round(sqrt((XS-xr)*(XS-xr)+(YS-yr)*(YS-yr)+(ZS-zr)*(ZS-zr))*1000.0);
        tr.gelev=(signed int)round(yr*1000.0);
        tr.sdepth=(signed int)round(YS*1000.0);   /* source depth (positive) */
        /* angle between receiver position and reference point
         (sperical coordinate system: swdep=theta, gwdep=phi) */
        tr.gdel=(signed int)round(atan2(-y,z)*180*1000.0/3.1415926);
        tr.gwdep=(signed int)round(sqrt(z*z+y*y)*1000.0);	
        tr.swdep=round(((360.0/(2.0*3.1415926))*atan2(x-xshift,y-yshift))*1000.0);
        tr.scalel=(signed short)-3;
        tr.scalco=(signed short)-3;
        tr.sx=(signed int)round(XS*1000.0);  /* X source coordinate */
        tr.sy=(signed int)round(YS*1000.0);  /* Y source coordinate */

        /* group coordinates */
        tr.gx=(signed int)round(xr*1000.0);
        tr.gy=(signed int)round(yr*1000.0);
        tr.ns=(unsigned short)ns; /* number of samples in this trace */
        tr.dt=(unsigned short)round(dtms);/* sample interval in micro-seconds */
        tr.d1=(float)tr.dt*1.0e-6;        /* sample spacing for non-seismic data */
        
        data.getRow(tempdata,tracl1);
        scai::lama::Scalar tempScalar;	
        for(IndexType sample=0; sample<tempdata.size(); sample++){
            tempScalar=tempdata.getValue(sample);
            tr.data[sample]=float(tempScalar.getValue<ValueType>());   
        }
        
        fwrite(&tr,240,1,pFile);
        fwrite(&tr.data[1],4,ns,pFile);
    }
    fclose(pFile);
}


