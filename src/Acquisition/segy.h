#pragma once

// \cond
//The file segy.h is property of the Colorado School of Mines.
//
//Copyright  2007, Colorado School of Mines,
//All rights reserved.
//
//
//Redistribution and use in source and binary forms, with or
//without modification, are permitted provided that the following
//conditions are met:
//
//*  Redistributions of source code must retain the above copyright
//notice, this list of conditions and the following disclaimer.
//*  Redistributions in binary form must reproduce the above
//copyright notice, this list of conditions and the following
//disclaimer in the documentation and/or other materials provided
//with the distribution.
//*  Neither the name of the Colorado School of Mines nor the names of
//its contributors may be used to endorse or promote products
//derived from this software without specific prior written permission.
//
//Warranty Disclaimer:
//THIS SOFTWARE IS PROVIDED BY THE COLORADO SCHOOL OF MINES AND CONTRIBUTORS
//"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
//LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
//FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
//COLORADO SCHOOL OF MINES OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
//INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
//                                                          BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
//                                                          LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
//CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
//STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
//IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
//POSSIBILITY OF SUCH DAMAGE.
//
//
//Export Restriction Disclaimer:
//We believe that CWP/SU: Seismic Un*x is a low technology product that does
//not appear on the Department of Commerce CCL list of restricted exports.
//Accordingly, we believe that our product meets the qualifications of
//an ECCN (export control classification number) of EAR99 and we believe
//it fits the qualifications of NRR (no restrictions required), and
//is thus not subject to export restrictions of any variety.
//
//Approved Reference Format:
//In publications, please refer to SU as per the following example:
//Cohen, J. K. and Stockwell, Jr. J. W., (200_), CWP/SU: Seismic Un*x
//Release No. __: an open source software  package for seismic
//research and processing,
//Center for Wave Phenomena, Colorado School of Mines.
//
//Articles about SU in peer-reviewed journals:
//Saeki, T., (1999), A guide to Seismic Un*x (SU)(2)---examples of data processing (part 1), data input and preparation of headers, Butsuri-Tansa (Geophysical Exploration), vol. 52, no. 5, 465-477.
//Stockwell, Jr. J. W. (1999), The CWP/SU: Seismic Un*x Package, Computers and Geosciences, May 1999.
//Stockwell, Jr. J. W. (1997), Free Software in Education: A case study of CWP/SU: Seismic Un*x, The Leading Edge, July 1997.
//Templeton, M. E., Gough, C.A., (1998), Web Seismic Un*x: Making seismic reflection processing more accessible, Computers and Geosciences.
//
//Acknowledgements:
//SU stands for CWP/SU:Seismic Un*x, a processing line developed at Colorado
//School of Mines, partially based on Stanford Exploration Project (SEP)
//software.


//! \brief Define the structure for Seismic Unix (SEG-Y) header
#define SU_NFLTS	USHRT_MAX	/* Arbitrary limit on data array size	*/
struct Segy{	/* segy - trace identification header */
    
    signed tracl   :32;	/* trace sequence number within line */
    
    signed tracr   :32;	/* trace sequence number within reel */
    
    signed fldr    :32;	/* field record number */
    
    signed tracf   :32;	/* trace number within field record */
    
    signed ep      :32;	/* energy source point number */
    
    signed cdp     :32;	/* CDP ensemble number */
    
    signed cdpt    :32;	/* trace number within CDP ensemble */
    
    signed trid    :16;	/* trace identification code:
                         1 = seismic data
                         2 = dead
                         3 = dummy
                         4 = time break
                         5 = uphole
                         6 = sweep
                         7 = timing
                         8 = water break
                         9---, N = optional use (N = 32,767)
                         
                         Following are CWP id flags:
                         
                         9 = autocorrelation
                         
                         10 = Fourier transformed - no packing
                         xr[0],xi[0], ..., xr[N-1],xi[N-1]
                         
                         11 = Fourier transformed - unpacked Nyquist
                         xr[0],xi[0],...,xr[N/2],xi[N/2]
                         
                         12 = Fourier transformed - packed Nyquist
                         even N:
                         xr[0],xr[N/2],xr[1],xi[1], ...,
                         xr[N/2 -1],xi[N/2 -1]
                         (note the exceptional second entry)
                         odd N:
                         xr[0],xr[(N-1)/2],xr[1],xi[1], ...,
                         xr[(N-1)/2 -1],xi[(N-1)/2 -1],xi[(N-1)/2]
                         (note the exceptional second & last entries)
                         
                         13 = Complex signal in the time domain
                         xr[0],xi[0], ..., xr[N-1],xi[N-1]
                         
                         14 = Fourier transformed - amplitude/phase
                         a[0],p[0], ..., a[N-1],p[N-1]
                         
                         15 = Complex time signal - amplitude/phase
                         a[0],p[0], ..., a[N-1],p[N-1]
                         
                         16 = Real part of complex trace from 0 to Nyquist
                         
                         17 = Imag part of complex trace from 0 to Nyquist
                         
                         18 = Amplitude of complex trace from 0 to Nyquist
                         
                         19 = Phase of complex trace from 0 to Nyquist
                         
                         21 = Wavenumber time domain (k-t)
                         
                         22 = Wavenumber frequency (k-omega)
                         
                         23 = Envelope of the complex time trace
                         
                         24 = Phase of the complex time trace
                         
                         25 = Frequency of the complex time trace
                         
                         30 = Depth-Range (z-x) traces
                         
                         101 = Seismic data packed to bytes (by supack1)
                         
                         102 = Seismic data packed to 2 bytes (by supack2)
                         */
    
    signed nvs    :16;   /* number of vertically summed traces (see vscode
                          in bhed structure) */
    
    signed nhs    :16;   /* number of horizontally summed traces (see vscode
                          in bhed structure) */
    
    signed duse   :16;   /* data use:
                          1 = production
                          2 = test */
    
    signed offset :32; /* distance from source point to receiver
                        group (negative if opposite to direction
                        in which the line was shot) */
    
    signed gelev  :32; /* receiver group elevation from sea level
                        (above sea level is positive) */
    
    signed selev  :32; /* source elevation from sea level
                        (above sea level is positive) */
    
    signed sdepth :32; /* source depth (positive) */
    
    signed gdel   :32; /* datum elevation at receiver group */
    
    signed sdel   :32; /* datum elevation at source */
    
    signed swdep  :32; /* water depth at source */
    
    signed gwdep  :32; /* water depth at receiver group */
    
    signed scalel :16; /* scale factor for previous 7 entries
                        with value plus or minus 10 to the
                        power 0, 1, 2, 3, or 4 (if positive,
                        multiply, if negative divide) */
    
    signed scalco :16; /* scale factor for next 4 entries
                        with value plus or minus 10 to the
                        power 0, 1, 2, 3, or 4 (if positive,
                        multiply, if negative divide) */
    
    signed  sx    :32;   /* X source coordinate */
    
    signed  sy    :32;   /* Y source coordinate */
    
    signed  gx    :32;   /* X group coordinate */
    
    signed  gy    :32;   /* Y group coordinate */
    
    signed counit :16;   /* coordinate units code:
                          for previous four entries
                          1 = length (meters or feet)
                          2 = seconds of arc (in this case, the
                          X values are longitude and the Y values
                          are latitude, a positive value designates
                          the number of seconds east of Greenwich
                          or north of the equator */
    
    signed wevel  :16;	/* weathering velocity */
    
    signed swevel :16;	/* subweathering velocity */
    
    signed sut    :16;	/* uphole time at source */
    
    signed gut    :16;	/* uphole time at receiver group */
    
    signed sstat  :16;	/* source static correction */
    
    signed gstat  :16;	/* group static correction */
    
    signed tstat  :16;	/* total static applied */
    
    signed laga   :16; /* lag time A, time in ms between end of 240-
                        byte trace identification header and time
                        break, positive if time break occurs after
                        end of header, time break is defined as
                        the initiation pulse which maybe recorded
                        on an auxiliary trace or as otherwise
                        specified by the recording system */
    
    signed lagb   :16; /* lag time B, time in ms between the time break
                        and the initiation time of the energy source,
                        may be positive or negative */
    
    signed delrt  :16; /* delay recording time, time in ms between
                        initiation time of energy source and time
                        when recording of data samples begins
                        (for deep water work if recording does not
                        start at zero time) */
    
    signed muts   :16; /* mute time--start */
    
    signed mute   :16; /* mute time--end */
    
    unsigned ns   :16; /* number of samples in this trace */
    
    unsigned dt   :16; /* sample interval; in micro-seconds */
    
    signed gain   :16; /* gain type of field instruments code:
                        1 = fixed
                        2 = binary
                        3 = floating point
                        4 ---- N = optional use */
    
    signed igc    :16; /* instrument gain constant */
    
    signed igi    :16; /* instrument early or initial gain */
    
    signed corr   :16; /* correlated:
                        1 = no
                        2 = yes */
    
    signed sfs    :16; /* sweep frequency at start */
    
    signed sfe    :16; /* sweep frequency at end */
    
    signed slen   :16; /* sweep length in ms */
    
    signed styp   :16; /* sweep type code:
                        1 = linear
                        2 = cos-squared
                        3 = other */
    
    signed stas   :16; /* sweep trace length at start in ms */
    
    signed stae   :16; /* sweep trace length at end in ms */
    
    signed tatyp  :16; /* taper type: 1=linear, 2=cos^2, 3=other */
    
    signed afilf  :16; /* alias filter frequency if used */
    
    signed afils  :16; /* alias filter slope */
    
    signed nofilf :16; /* notch filter frequency if used */
    
    signed nofils :16; /* notch filter slope */
    
    signed lcf    :16; /* low cut frequency if used */
    
    signed hcf    :16; /* high cut frequncy if used */
    
    signed lcs    :16; /* low cut slope */
    
    signed hcs    :16; /* high cut slope */
    
    signed year   :16; /* year data recorded */
    
    signed day    :16; /* day of year */
    
    signed hour   :16; /* hour of day (24 hour clock) */
    
    signed minute :16; /* minute of hour */
    
    signed sec    :16; /* second of minute */
    
    signed timbas :16; /* time basis code:
                        1 = local
                        2 = GMT
                        3 = other */
    
    signed trwf   :16; /* trace weighting factor, defined as 1/2^N
                        volts for the least sigificant bit */
    
    signed grnors :16; /* geophone group number of roll switch
                        position one */
    
    signed grnofr :16; /* geophone group number of trace one within
                        original field record */
    
    signed grnlof :16; /* geophone group number of last trace within
                        original field record */
    
    signed gaps   :16;  /* gap size (total number of groups dropped) */
    
    signed otrav  :16;  /* overtravel taper code:
                         1 = down (or behind)
                         2 = up (or ahead) */
    
    /* local assignments */
    
    /*        signed pad :32; */ /* double word alignment for Cray 64-bit floats */
    
    float d1;	/* sample spacing for non-seismic data */
    
    float f1;	/* first sample location for non-seismic data */
    
    float d2;	/* sample spacing between traces */
    
    float f2;	/* first trace location */
    
    float ungpow;	/* negative of power used for dynamic
                     range compression */
    
    float unscale;	/* reciprocal of scaling factor to normalize
                     range */
    signed ntr   :32;   /* number of traces */
    
    signed mark  :16;   /* mark selected traces */
    
    signed unass :16;   /* unassigned values */
    
    float  data[SU_NFLTS];
    
};

// \endcond
