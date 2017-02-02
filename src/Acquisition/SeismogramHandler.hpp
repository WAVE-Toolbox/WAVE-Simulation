#pragma once

#include "../Configuration/Configuration.hpp"
#include "Acquisition.hpp"
#include "Seismogram.hpp"

namespace KITGPI
{

    namespace Acquisition
    {

        /*! \brief class for a single Seismogram
         *
         * This class handels a single seismogram.
         */
        template <typename ValueType>
        class SeismogramHandler
        {

          public:
            //! \brief Default constructor
            explicit SeismogramHandler();
            //! \brief Default destructor
            ~SeismogramHandler(){};

            void writeToFileRaw(std::string const &filename) const;
            void write(Configuration::Configuration const &config) const;
            void resetData();

            void setSourceCoordinate(IndexType sourceCoord);
            void setDT(ValueType newDT);
            void setContextPtr(hmemo::ContextPtr ctx);

            inline Seismogram<ValueType> const &getSeismogram(SeismogramType type) const;
            inline Seismogram<ValueType> &getSeismogram(SeismogramType type);
            inline IndexType getNumTracesGlobal(SeismogramType type) const;
            inline IndexType getNumSamples(SeismogramType type) const;

          private:
            void setTraceType();

            std::vector<Seismogram<ValueType>> seismo; //!< vector in which the seismogram is stored
        };
    }
}

//! \brief write the seismogram.
/*!
 \param config Configuration
 */
template <typename ValueType>
void KITGPI::Acquisition::SeismogramHandler<ValueType>::write(Configuration::Configuration const &config) const
{
    for (auto const &i : seismo) {
        i.write(config);
    }
}

/*! \brief Default constructor. */
template <typename ValueType>
KITGPI::Acquisition::SeismogramHandler<ValueType>::SeismogramHandler()
    : seismo(NUM_ELEMENTS_SEISMOGRAMTYPE)
{
    seismo.shrink_to_fit();
    setTraceType();
}

/*! \brief Setter method for trace type. */
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

//! \brief Constant methode to write to RAW-file.
/*!
 \param filename Filename of the output file
 */
template <typename ValueType>
void KITGPI::Acquisition::SeismogramHandler<ValueType>::writeToFileRaw(std::string const &filename) const
{
    for (auto const &i : seismo) {
        i.writeToFileRaw(filename);
    }
}

/*! \brief Method to reset data. */
template <typename ValueType>
void KITGPI::Acquisition::SeismogramHandler<ValueType>::resetData()
{
    for (auto &i : seismo) {
        i.resetData();
    }
}

//! \brief Constant getter methode for number of samples.
/*!
 \param type Type of the seismogram
 */
template <typename ValueType>
IndexType KITGPI::Acquisition::SeismogramHandler<ValueType>::getNumSamples(SeismogramType type) const
{
    SCAI_ASSERT_ERROR(type >= 0 && type <= NUM_ELEMENTS_SEISMOGRAMTYPE - 1, "SeismogramType unkown");
    return (seismo[type].getNumSamples());
}

//! \brief Constant getter methode for number of global traces.
/*!
 \param type Type of the seismogram
 */
template <typename ValueType>
IndexType KITGPI::Acquisition::SeismogramHandler<ValueType>::getNumTracesGlobal(SeismogramType type) const
{
    SCAI_ASSERT_ERROR(type >= 0 && type <= NUM_ELEMENTS_SEISMOGRAMTYPE - 1, "SeismogramType unkown");
    return (seismo[type].getNumTracesGlobal());
}

//! \brief Constant getter methode for the seismogram.
/*!
 \param type Type of the seismogram
 */
template <typename ValueType>
KITGPI::Acquisition::Seismogram<ValueType> const &KITGPI::Acquisition::SeismogramHandler<ValueType>::getSeismogram(SeismogramType type) const
{
    SCAI_ASSERT_ERROR(type >= 0 && type <= NUM_ELEMENTS_SEISMOGRAMTYPE - 1, "SeismogramType unkown");
    return (seismo[type]);
}

//! \brief Getter methode for the seismogram.
/*!
 \param type Type of the seismogram
 */
template <typename ValueType>
KITGPI::Acquisition::Seismogram<ValueType> &KITGPI::Acquisition::SeismogramHandler<ValueType>::getSeismogram(SeismogramType type)
{
    SCAI_ASSERT_ERROR(type >= 0 && type <= NUM_ELEMENTS_SEISMOGRAMTYPE - 1, "SeismogramType unkown");
    return (seismo[type]);
}

//! \brief Methode to set the right context-pointer.
/*!
 \param ctx Context
 */
template <typename ValueType>
void KITGPI::Acquisition::SeismogramHandler<ValueType>::setContextPtr(hmemo::ContextPtr ctx)
{
    for (auto &i : seismo) {
        i.setContextPtr(ctx);
    }
}

//! \brief Methode to set DT.
/*!
 \param newDT Temporal sampling which will be set to the seismogram
 */
template <typename ValueType>
void KITGPI::Acquisition::SeismogramHandler<ValueType>::setDT(ValueType newDT)
{
    for (auto &i : seismo) {
        i.setDT(newDT);
    }
}

//! \brief Methode to set Source Coordinate.
/*!
 \param sourceCoord Source coordinate in 1-D format
 */
template <typename ValueType>
void KITGPI::Acquisition::SeismogramHandler<ValueType>::setSourceCoordinate(IndexType sourceCoord)
{
    for (auto &i : seismo) {
        i.setSourceCoordinate(sourceCoord);
    }
}
