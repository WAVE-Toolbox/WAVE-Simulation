#include "SeismogramHandler.hpp"
using namespace scai;

//! \brief Method to write all handled Seismogram to file
/*!
 * This method allows to write all handled Seismogram to file. The Configuration class will be used to determine requiered information,
 * such as filename or header information. The #SeismogramType of each Seismogram will be added to the filename automaticly.
 *
 \param config Configuration class, which will be used to determine filename and header information
 */
template <typename ValueType>
void KITGPI::Acquisition::SeismogramHandler<ValueType>::write(Configuration::Configuration const &config, std::string const &filename) const
{
    for (auto const &i : seismo) {
        i.write(config, filename);
    }
}

//! \brief Method to write normalize Seismogram-traces
/*!
 *
 */
template <typename ValueType>
void KITGPI::Acquisition::SeismogramHandler<ValueType>::normalize()
{
    for (auto &i : seismo) {
        i.getNormalizeTraces();
        i.normalizeTrace();
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

//! \brief Write all handled Seismogram to a raw-file.
/*!
 *
 * This method writes all handled Seismogram to a raw-file (MTX format, without header). The #SeismogramType will be added to the filenames automaticly.
 \param filename Filename of the output file
 */
template <typename ValueType>
void KITGPI::Acquisition::SeismogramHandler<ValueType>::writeToFileRaw(std::string const &filename) const
{
    for (auto const &i : seismo) {
        i.writeToFileRaw(filename);
    }
}

//! \brief Read all handled Seismogram from a raw-file.
/*!
 *
 * This method reads all handled Seismogram from a raw-file (MTX format, without header). The #SeismogramType will be added to the filenames automaticaly.
 \param filename Filename of the input file
 */
template <typename ValueType>
void KITGPI::Acquisition::SeismogramHandler<ValueType>::readFromFileRaw(std::string const &filename, bool copyDist)
{
    for (auto &i : seismo) {
        if (i.getNumTracesGlobal() > 0)
            i.readFromFileRaw(filename, copyDist);
    }
}

//! \brief Read all handled Seismogram from a raw-file.
/*!
 *
 * This method reads all handled Seismogram from a raw-file (MTX format, without header). The #SeismogramType will be added to the filenames automaticaly.
 \param filename Filename of the input file
 */
template <typename ValueType>
void KITGPI::Acquisition::SeismogramHandler<ValueType>::readFromFileRaw(std::string const &filename, scai::dmemo::DistributionPtr distTraces, scai::dmemo::DistributionPtr distSamples)
{
    for (auto &i : seismo) {
        if (i.getNumTracesGlobal() > 0)
            i.readFromFileRaw(filename, distTraces, distSamples);
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

//! \brief Setter methode to set DT.
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

//! \brief Setter methode to set Index for trace-normalization.
/*!
 *
 * This method sets the index for trace-normalization.
 \param normalizeTrace Index for trace-normalization which will normalize the seismogram traces
 */
template <typename ValueType>
void KITGPI::Acquisition::SeismogramHandler<ValueType>::setNormalizeTraces(IndexType normalize)
{
    for (auto &i : seismo) {
        i.setNormalizeTraces(normalize);
    }
}

//! \brief Setter methode to set source coordinate
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

template class KITGPI::Acquisition::SeismogramHandler<double>;
template class KITGPI::Acquisition::SeismogramHandler<float>;
