#pragma once

namespace KITGPI
{

    /*! \brief Acquisition namespace
     *
     * This namespace handles the seismic acquisition, e.g. Sources, Receivers, Seismogram and SourceSignals.
     */
    namespace Acquisition
    {

        /*! \brief List of seismogram types
         *
         * Types of seismograms, which are used in this programm.
         */
        enum SeismogramType { P,  //!< Pressure
                              VX, //!< Horizontal inline velocity
                              VY, //!< Vertical particle velocity
                              VZ  //!< Horizontal crossline particle velocity
        };

        /*! \brief Map to transfer a SeismogramType to the corresponding string
         * 
         * This map can be used to map a SeismogramType to the corresponding char.
         * This can be usefull to add the SeismogramType to a filename etc.
         */
        static std::map<SeismogramType, const char *> SeismogramTypeString = {
            {P, "p"},
            {VX, "vx"},
            {VY, "vy"},
            {VZ, "vz"}};

        constexpr scai::IndexType NUM_ELEMENTS_SEISMOGRAMTYPE = 4; //!< Number of seismogram types, which are implemented.
    }
}
