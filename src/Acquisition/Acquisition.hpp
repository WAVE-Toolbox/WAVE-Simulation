
#pragma once

namespace KITGPI {
    
    namespace Acquisition {
        
        enum SeismogramType { P, VX, VY, VZ };
        
        static std::map< SeismogramType, const char * > SeismogramTypeString = {
            {P, "p"},
            {VX, "vx"},
            {VY, "vy"},
            {VZ, "vz"}
        };
        
        constexpr IndexType NUM_ELEMENTS_SEISMOGRAMTYPE=4;
        
    }
}
