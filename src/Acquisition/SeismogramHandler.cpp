#include "SeismogramHandler.hpp"
#include "../Common/HostPrint.hpp"
using namespace scai;

//! \brief Method to write all handled Seismogram to file
/*!
 * This method allows to write all handled Seismogram to file. 
 * The #SeismogramType of each Seismogram will be added to the filename automaticly.
 * 
 *
 \param seismogramFormat =1 MTX: MatrixMaker format, =4 SU: SeismicUnix format
 \param filename base filename of the seismogram
 \param modelCoordinates Coordinate class, which eg. maps 3D coordinates to 1D model indices
 */
template <typename ValueType>
void KITGPI::Acquisition::SeismogramHandler<ValueType>::write(scai::IndexType const seismogramFormat, std::string const &filename, Coordinates<ValueType> const &modelCoordinates)
{
    for (auto &i : seismo) {
        i.write(seismogramFormat, filename, modelCoordinates);
    }
}

//! \brief Method to read all handled Seismograms from file
/*!
 * This method allows to read all handled Seismogram from file. 
 * The #SeismogramType of each Seismogram will be added to the filename automaticly.
 *
 \param seismogramFormat =1 MTX: MatrixMaker format, =4 SU: SeismicUnix format
 \param filename base filename of the seismograms
 \param copyDist Boolean: 0 = read data undistributed (default), data is replicated on each process // 1 = read data with existing distribution of data
 */
template <typename ValueType>
void KITGPI::Acquisition::SeismogramHandler<ValueType>::read(IndexType const seismogramFormat, std::string const &filename, bool copyDist)
{
    for (auto &i : seismo) {
        if (i.getNumTracesGlobal() > 0)
            i.read(seismogramFormat, filename, copyDist);
    }
}

//! \brief Method to normalize Seismogram-traces
/*!
 *
 */
template <typename ValueType>
void KITGPI::Acquisition::SeismogramHandler<ValueType>::normalize(scai::IndexType normalizeTraces)
{
    for (auto &i : seismo) {
        i.normalizeTrace(normalizeTraces);
    }
}

//! \brief Method to integrate the Seismogram-traces
/*!
 *
 */
template <typename ValueType>
void KITGPI::Acquisition::SeismogramHandler<ValueType>::integrate()
{
    for (auto &i : seismo) {
        i.integrateTraces();
    }
}

//! \brief Method to differentiate the Seismogram-traces
/*!
 *
 */
template <typename ValueType>
void KITGPI::Acquisition::SeismogramHandler<ValueType>::differentiate()
{
    for (auto &i : seismo) {
        i.differentiateTraces();
    }
}

/*! \brief Method to filter the Seismogram-traces
 \param freqFilter filter object
*/
template <typename ValueType>
void KITGPI::Acquisition::SeismogramHandler<ValueType>::filter(Filter::Filter<ValueType> const &freqFilter)
{
    for (auto &i : seismo) {
        i.filterTraces(freqFilter);
    }
}

/*! \brief Constructor
 *
 * This constructor will initialize the handled Seismogram.
 * The number of handled Seismogram will be determined by the variable #NUM_ELEMENTS_SEISMOGRAMTYPE.
 * Since #NUM_ELEMENTS_SEISMOGRAMTYPE is equal to the number of different #SeismogramType, the SeismogramHandler will store a Seismogram for each #SeismogramType.
 * Moreover, the #SeismogramType will be set to all handled Seismogram, so each Seismogram knows which #SeismogramType it is.
 */
template <typename ValueType>
KITGPI::Acquisition::SeismogramHandler<ValueType>::SeismogramHandler()
    : seismo(NUM_ELEMENTS_SEISMOGRAMTYPE)
// seismo(4) initializes 4 Seismogram objects by calling the default constructor of the Class Seismogram. Other constructors can be called by e.g. seismo(4,Seismogram<ValueType>(a,b...))
{
    // seismo.shrink_to_fit();
    setTraceType();
    setTraceTypeEM();
}

/*! \brief Setter method for #SeismogramType
 *
 * This method sets the #SeismogramType to each handled Seismogram. The #SeismogramType is assigned in the same order the #SeismogramType are defined in the corresponding enum.
 */
template <typename ValueType>
void KITGPI::Acquisition::SeismogramHandler<ValueType>::setTraceType()
{
    SCAI_ASSERT_DEBUG(static_cast<SeismogramType>(0) == SeismogramType::P, "Cast went wrong");
    SCAI_ASSERT_DEBUG(static_cast<SeismogramType>(1) == SeismogramType::VX, "Cast went wrong");
    SCAI_ASSERT_DEBUG(static_cast<SeismogramType>(2) == SeismogramType::VY, "Cast went wrong");
    SCAI_ASSERT_DEBUG(static_cast<SeismogramType>(3) == SeismogramType::VZ, "Cast went wrong");

    IndexType count = 0;
    for (auto &i : seismo) {
        i.setTraceType(static_cast<SeismogramType>(count));
        ++count;
    }
}

/*! \brief Setter method for #SeismogramTypeEM
 *
 * This method sets the #SeismogramTypeEM to each handled Seismogram. The #SeismogramTypeEM is assigned in the same order the #SeismogramTypeEM are defined in the corresponding enum.
 */
template <typename ValueType>
void KITGPI::Acquisition::SeismogramHandler<ValueType>::setTraceTypeEM()
{

    SCAI_ASSERT_DEBUG(static_cast<SeismogramTypeEM>(0) == SeismogramTypeEM::EZ, "Cast went wrong");
    SCAI_ASSERT_DEBUG(static_cast<SeismogramTypeEM>(1) == SeismogramTypeEM::EX, "Cast went wrong");
    SCAI_ASSERT_DEBUG(static_cast<SeismogramTypeEM>(2) == SeismogramTypeEM::EY, "Cast went wrong");
    SCAI_ASSERT_DEBUG(static_cast<SeismogramTypeEM>(3) == SeismogramTypeEM::HZ, "Cast went wrong");

    IndexType count = 0;
    for (auto &i : seismo) {
        i.setTraceTypeEM(static_cast<SeismogramTypeEM>(count));
        ++count;
    }
}

/*! \brief Method to reset the Seismogram data
 *
 * This method resets the Seismogram data (the content of the traces). The reset will set all values to zero.
 * However, the memory (number of samples * number of traces * number of Seismograms) will stay allocated, only the values are set to zero.
 */
template <typename ValueType>
void KITGPI::Acquisition::SeismogramHandler<ValueType>::resetData()
{
    for (auto &i : seismo) {
        i.resetData();
    }
}

/*! \brief Method to reset all Seismogram (components)
 *
 * This method clears the Seismogram data 
 */
template <typename ValueType>
void KITGPI::Acquisition::SeismogramHandler<ValueType>::resetSeismograms()
{
    for (auto &i : seismo) {
        i.resetSeismogram();
    }
}
//! \brief Getter method for the number of samples
/*!
 * This method returns the number of temporal samples of the specified #SeismogramType.\n
 * See getSeismogram() for information how to use #SeismogramType as input parameter.
 \param type #SeismogramType of the desired Seismogram
 */
template <typename ValueType>
IndexType KITGPI::Acquisition::SeismogramHandler<ValueType>::getNumSamples(SeismogramType type) const
{
    SCAI_ASSERT_ERROR(type >= 0 && type <= NUM_ELEMENTS_SEISMOGRAMTYPE - 1, "SeismogramType unkown");
    return (seismo[type].getNumSamples());
}

//! \brief Getter method for the number of samples
/*!
 * This method returns the number of temporal samples of the specified #SeismogramTypeEM.\n
 * See getSeismogram() for information how to use #SeismogramTypeEM as input parameter.
 \param type #SeismogramTypeEM of the desired Seismogram
 */
template <typename ValueType>
IndexType KITGPI::Acquisition::SeismogramHandler<ValueType>::getNumSamples(SeismogramTypeEM type) const
{
    SCAI_ASSERT_ERROR(type >= 0 && type <= NUM_ELEMENTS_SEISMOGRAMTYPE - 1, "SeismogramTypeEM unkown");
    return (seismo[type].getNumSamples());
}

//! \brief Getter method for the total number of traces
/*!
 * This method returns the total number of global traces, which is the sum of the number of global traces of all handled #Seismogram.\n
 \return Total number of handled traces
 */
template <typename ValueType>
IndexType KITGPI::Acquisition::SeismogramHandler<ValueType>::getNumTracesTotal() const
{
    IndexType sum = 0;
    for (auto &i : seismo) {
        sum += i.getNumTracesGlobal();
    }
    return sum;
}

//! \brief Getter method for the number of global traces
/*!
 * This method returns the number of global traces for the specified #SeismogramType.\n
* See getSeismogram() for information how to use #SeismogramType as input parameter.
 \param type #SeismogramType of the desired Seismogram
 */
template <typename ValueType>
IndexType KITGPI::Acquisition::SeismogramHandler<ValueType>::getNumTracesGlobal(SeismogramType type) const
{
    SCAI_ASSERT_ERROR(type >= 0 && type <= NUM_ELEMENTS_SEISMOGRAMTYPE - 1, "SeismogramType unkown");
    return (seismo[type].getNumTracesGlobal());
}

//! \brief Getter method for the number of global traces
/*!
 * This method returns the number of global traces for the specified #SeismogramTypeEM.\n
* See getSeismogram() for information how to use #SeismogramTypeEM as input parameter.
 \param type #SeismogramTypeEM of the desired Seismogram
 */
template <typename ValueType>
IndexType KITGPI::Acquisition::SeismogramHandler<ValueType>::getNumTracesGlobal(SeismogramTypeEM type) const
{
    SCAI_ASSERT_ERROR(type >= 0 && type <= NUM_ELEMENTS_SEISMOGRAMTYPE - 1, "SeismogramTypeEM unkown");
    return (seismo[type].getNumTracesGlobal());
}

//! \brief Getter method for the Seismogram class
/*!
 * This method returns the Seismogram class for the specified #SeismogramType.
 *
 * **Example Usage:**\n
 * To get const access to the pressure seismogram (SeismogramType::P):\n
 * `const Acquisition::Seismogram pressure = handler.getSeismogram(Acquisition::SeismogramType::P);`\n
 * or to the VX seismogram (SeismogramType::VX):\n
 * `const Acquisition::Seismogram vx = handler.getSeismogram(Acquisition::SeismogramType::VX);`\n
 \param type #SeismogramType of the desired Seismogram
 */
template <typename ValueType>
KITGPI::Acquisition::Seismogram<ValueType> const &KITGPI::Acquisition::SeismogramHandler<ValueType>::getSeismogram(SeismogramType type) const
{
    SCAI_ASSERT_ERROR(type >= 0 && type <= NUM_ELEMENTS_SEISMOGRAMTYPE - 1, "SeismogramType unkown");
    return (seismo[type]);
}

//! \brief Getter method for the Seismogram class
/*!
 * This method returns the Seismogram class for the specified #SeismogramTypeEM.
 *
 * **Example Usage:**\n
 * To get const access to the pressure seismogram (SeismogramTypeEM::EZ):\n
 * `const Acquisition::Seismogram pressure = handler.getSeismogram(Acquisition::SeismogramTypeEM::EZ);`\n
 * or to the HX seismogram (SeismogramTypeEM::HX):\n
 * `const Acquisition::Seismogram hx = handler.getSeismogram(Acquisition::SeismogramTypeEM::HX);`\n
 \param type #SeismogramTypeEM of the desired Seismogram
 */
template <typename ValueType>
KITGPI::Acquisition::Seismogram<ValueType> const &KITGPI::Acquisition::SeismogramHandler<ValueType>::getSeismogram(SeismogramTypeEM type) const
{
    SCAI_ASSERT_ERROR(type >= 0 && type <= NUM_ELEMENTS_SEISMOGRAMTYPE - 1, "SeismogramTypeEM unkown");
    return (seismo[type]);
}

//! \brief Getter method for the Seismogram class
/*!
 * This method returns the Seismogram class for the specified #SeismogramType.
 *
 * **Example Usage:**\n
 * To get const access to the pressure seismogram (SeismogramType::P):\n
 * `Acquisition::Seismogram pressure = handler.getSeismogram(Acquisition::SeismogramType::P);`\n
 * or to the VX seismogram (SeismogramType::VX):\n
 * `Acquisition::Seismogram vx = handler.getSeismogram(Acquisition::SeismogramType::VX);`\n
 \param type #SeismogramType of the desired Seismogram
 */
template <typename ValueType>
KITGPI::Acquisition::Seismogram<ValueType> &KITGPI::Acquisition::SeismogramHandler<ValueType>::getSeismogram(SeismogramType type)
{
    SCAI_ASSERT_ERROR(type >= 0 && type <= NUM_ELEMENTS_SEISMOGRAMTYPE - 1, "SeismogramType unkown");
    return (seismo[type]);
}

//! \brief Getter method for the Seismogram class
/*!
 * This method returns the Seismogram class for the specified #SeismogramTypeEM.
 *
 * **Example Usage:**\n
 * To get const access to the pressure seismogram (SeismogramTypeEM::EZ):\n
 * `Acquisition::Seismogram pressure = handler.getSeismogram(Acquisition::SeismogramTypeEM::EZ);`\n
 * or to the HX seismogram (SeismogramTypeEM::HX):\n
 * `Acquisition::Seismogram hx = handler.getSeismogram(Acquisition::SeismogramTypeEM::HX);`\n
 \param type #SeismogramTypeEM of the desired Seismogram
 */
template <typename ValueType>
KITGPI::Acquisition::Seismogram<ValueType> &KITGPI::Acquisition::SeismogramHandler<ValueType>::getSeismogram(SeismogramTypeEM type)
{
    SCAI_ASSERT_ERROR(type >= 0 && type <= NUM_ELEMENTS_SEISMOGRAMTYPE - 1, "SeismogramTypeEM unkown");
    return (seismo[type]);
}

//! \brief Setter method for the context type
/*!
 * This method sets the context to all handled Seismogram.
 \param ctx Context
 */
template <typename ValueType>
void KITGPI::Acquisition::SeismogramHandler<ValueType>::setContextPtr(scai::hmemo::ContextPtr ctx)
{
    for (auto &i : seismo) {
        i.setContextPtr(ctx);
    }
}

//! \brief Setter method to set DT.
/*!
 *
 * This method sets the temporal sampling DT to all handled Seismogram.
 \param newDT Temporal sampling which will be set to the seismogram
 */
template <typename ValueType>
void KITGPI::Acquisition::SeismogramHandler<ValueType>::setDT(ValueType newDT)
{
    for (auto &i : seismo) {
        i.setDT(newDT);
    }
}

//! \brief Setter method to set source coordinate
/*!
 * This method sets the source coordinate to all handled Seismogram.
 \param sourceCoord Source coordinate in 1-D format
 */
template <typename ValueType>
void KITGPI::Acquisition::SeismogramHandler<ValueType>::setSourceCoordinate(IndexType sourceCoord)
{
    for (auto &i : seismo) {
        i.setSourceCoordinate(sourceCoord);
    }
}

//! \brief Setter method to set matrix for resampling this seismogram.
/*!
 \param rMat Resampling matrix
 */
template <typename ValueType>
void KITGPI::Acquisition::SeismogramHandler<ValueType>::setSeismoDT(ValueType seismoDT)
{
    for (auto &i : seismo) {
        i.setSeismoDT(seismoDT);
    }
}

//! \brief Setter method to set outputInstantaneous.
/*!
 \param instantaneousTraces outputInstantaneous
 */
template <typename ValueType>
void KITGPI::Acquisition::SeismogramHandler<ValueType>::setInstantaneousTrace(scai::IndexType instantaneousTraces)
{
    for (auto &i : seismo) {
        i.setInstantaneousTrace(instantaneousTraces);
    }
}

//! \brief Setter method to set frequencyAGC.
/*!
 \param setFrequencyAGC setFrequencyAGC
 */
template <typename ValueType>
void KITGPI::Acquisition::SeismogramHandler<ValueType>::setFrequencyAGC(ValueType setFrequencyAGC)
{
    for (auto &i : seismo) {
        i.setFrequencyAGC(setFrequencyAGC);
    }
}

//! \brief calculate InverseAGC .
/*!
 */
template <typename ValueType>
void KITGPI::Acquisition::SeismogramHandler<ValueType>::calcInverseAGC()
{
    for (auto &i : seismo) {
        i.calcInverseAGC();
    }
}

//! \brief Setter method to set inverseAGC.
/*!
 \param setFrequencyAGC setFrequencyAGC
 */
template <typename ValueType>
void KITGPI::Acquisition::SeismogramHandler<ValueType>::setInverseAGC(SeismogramHandler<ValueType> seismograms)
{
    for (int i = 0; i < NUM_ELEMENTS_SEISMOGRAMTYPE; i++) {
        seismo[i].setInverseAGC(seismograms.getSeismogram(static_cast<Acquisition::SeismogramType>(i)).getInverseAGC());
    }
}

//! \brief Setter method to set isSeismic
/*!
 */
template <typename ValueType>
void KITGPI::Acquisition::SeismogramHandler<ValueType>::setIsSeismic(bool setIsSeismic)
{
    isSeismic = setIsSeismic;
    for (auto &i : seismo) {
        i.setIsSeismic(setIsSeismic);
    }
}

//! \brief Getter method to get isSeismic
/*!
 */
template <typename ValueType>
bool KITGPI::Acquisition::SeismogramHandler<ValueType>::getIsSeismic() const
{
    return isSeismic;
}

/*! \brief Overloading *= Operation
 *
 \param rhs Scalar factor with which the vectors are multiplied.
 */
template <typename ValueType>
KITGPI::Acquisition::SeismogramHandler<ValueType> KITGPI::Acquisition::SeismogramHandler<ValueType>::operator*=(scai::lama::DenseVector<ValueType> const &rhs)
{
    for (auto &i : seismo) {
        i *= rhs;
    }

    return *this;
}

//! \brief Check seismograms for inf or NaN
/*!
 */
template <typename ValueType>
bool KITGPI::Acquisition::SeismogramHandler<ValueType>::isFinite()
{
    bool isfinite=true;
    for (auto &i : seismo) {
        isfinite = i.isFinite();
            if (isfinite==false){
                break;
            }
    }
    return(isfinite);
}

template class KITGPI::Acquisition::SeismogramHandler<double>;
template class KITGPI::Acquisition::SeismogramHandler<float>;
