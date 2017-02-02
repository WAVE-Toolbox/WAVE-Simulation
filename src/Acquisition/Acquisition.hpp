#pragma once
using namespace scai;

namespace KITGPI
{

    //! \brief Acquisition namespace
    namespace Acquisition
    {

        //! \brief list of seismogram types
        enum SeismogramType { P,
                              VX,
                              VY,
                              VZ };

        //! \brief seismogram types: pressure or velocity
        static std::map<SeismogramType, const char *> SeismogramTypeString = {
            {P, "p"},
            {VX, "vx"},
            {VY, "vy"},
            {VZ, "vz"}};

        constexpr IndexType NUM_ELEMENTS_SEISMOGRAMTYPE = 4; //!< Number of seismogram types
    }
}
