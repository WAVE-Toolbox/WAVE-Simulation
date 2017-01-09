#pragma once

#include "Seismogram.hpp"
#include "Acquisition.hpp"

namespace KITGPI {
    
    namespace Acquisition {
        
        template <typename ValueType>
        class SeismogramHandler
        {
            
        public:
            
            explicit SeismogramHandler();
            ~SeismogramHandler(){};
            
            void writeToFileRaw(std::string const& filename) const;
            void resetData();
            
            void setSourceCoordinate(IndexType sourceCoord);
            void setDT(ValueType newDT);
            void setContextPtr(hmemo::ContextPtr ctx);
            
            inline Seismogram<ValueType>const& getSeismogram(SeismogramType type) const;
            inline Seismogram<ValueType>& getSeismogram(SeismogramType type);
            inline IndexType getNumTracesGlobal(SeismogramType type) const;
            inline IndexType getNumSamples(SeismogramType type) const;
            
        private:
            
            void setTraceType();
            
            std::vector<Seismogram<ValueType>> seismo;
            
        };
    }
}

template <typename ValueType>
KITGPI::Acquisition::SeismogramHandler<ValueType>::SeismogramHandler()
:seismo(NUM_ELEMENTS_SEISMOGRAMTYPE)
{
    seismo.shrink_to_fit();
    setTraceType();
}

template <typename ValueType>
void KITGPI::Acquisition::SeismogramHandler<ValueType>::setTraceType()
{
    
    SCAI_ASSERT_DEBUG(static_cast<SeismogramType>(0)==SeismogramType::P, "Cast went wrong");
    SCAI_ASSERT_DEBUG(static_cast<SeismogramType>(1)==SeismogramType::VX, "Cast went wrong");
    SCAI_ASSERT_DEBUG(static_cast<SeismogramType>(2)==SeismogramType::VY, "Cast went wrong");
    SCAI_ASSERT_DEBUG(static_cast<SeismogramType>(3)==SeismogramType::VZ, "Cast went wrong");
    
    IndexType count=0;
    for(auto &i : seismo){
        i.setTraceType(static_cast<SeismogramType>(count));
        ++count;
    }
}

template <typename ValueType>
void KITGPI::Acquisition::SeismogramHandler<ValueType>::writeToFileRaw(std::string const& filename) const
{
    for(auto const& i : seismo){
        i.writeToFileRaw(filename);
    }
}

template <typename ValueType>
void KITGPI::Acquisition::SeismogramHandler<ValueType>::resetData()
{
    for(auto &i : seismo){
        i.resetData();
    }
}

template <typename ValueType>
IndexType KITGPI::Acquisition::SeismogramHandler<ValueType>::getNumSamples(SeismogramType type) const
{
    SCAI_ASSERT_ERROR(type >= 0 && type <= NUM_ELEMENTS_SEISMOGRAMTYPE-1, "SeismogramType unkown");
    return(seismo[type].getNumSamples());
}

template <typename ValueType>
IndexType KITGPI::Acquisition::SeismogramHandler<ValueType>::getNumTracesGlobal(SeismogramType type) const
{
    SCAI_ASSERT_ERROR(type >= 0 && type <= NUM_ELEMENTS_SEISMOGRAMTYPE-1, "SeismogramType unkown");
    return(seismo[type].getNumTracesGlobal());
}

template <typename ValueType>
KITGPI::Acquisition::Seismogram<ValueType>const& KITGPI::Acquisition::SeismogramHandler<ValueType>::getSeismogram(SeismogramType type) const
{
    SCAI_ASSERT_ERROR(type >= 0 && type <= NUM_ELEMENTS_SEISMOGRAMTYPE-1, "SeismogramType unkown");
    return(seismo[type]);
}

template <typename ValueType>
KITGPI::Acquisition::Seismogram<ValueType>& KITGPI::Acquisition::SeismogramHandler<ValueType>::getSeismogram(SeismogramType type)
{
    SCAI_ASSERT_ERROR(type >= 0 && type <= NUM_ELEMENTS_SEISMOGRAMTYPE-1, "SeismogramType unkown");
    return(seismo[type]);
}

template <typename ValueType>
void KITGPI::Acquisition::SeismogramHandler<ValueType>::setContextPtr(hmemo::ContextPtr ctx)
{
    for(auto &i : seismo){
        i.setContextPtr(ctx);
    }
}

template <typename ValueType>
void KITGPI::Acquisition::SeismogramHandler<ValueType>::setDT(ValueType newDT)
{
    for(auto &i : seismo){
        i.setDT(newDT);
    }
}

template <typename ValueType>
void KITGPI::Acquisition::SeismogramHandler<ValueType>::setSourceCoordinate(IndexType sourceCoord)
{
    for(auto &i : seismo){
        i.setSourceCoordinate(sourceCoord);
    }
}
