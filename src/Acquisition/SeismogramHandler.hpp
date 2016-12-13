#pragma once

#include "Seismogram.hpp"

namespace KITGPI {
    
    namespace Acquisition {
        
        
        enum SeismogramType { P, VX, VY, VZ };

        template <typename ValueType>
        class SeismogramHandler
        {
            
        public:
            
            explicit SeismogramHandler();
            ~SeismogramHandler(){};
            
            void resetData();
            
            void setSourceCoordinate(IndexType sourceCoord);
            void setDT(ValueType newDT);
            void setContextPtr(hmemo::ContextPtr ctx);
            
            inline Seismogram<ValueType>const& getSeismogram(SeismogramType type) const;
            inline Seismogram<ValueType>& getSeismogram(SeismogramType type);
            inline Seismogram<ValueType>& getSeismogram(IndexType type);
            inline IndexType getNumTracesGlobal(SeismogramType type) const;
            inline IndexType getNumSamples(SeismogramType type) const;
            
        private:
            
            std::vector<Seismogram<ValueType>> seismo;
            IndexType const NUM_ELEMENTS_ENUM=4;
            
        };
    }
}

template <typename ValueType>
KITGPI::Acquisition::SeismogramHandler<ValueType>::SeismogramHandler()
:seismo(4)
{
    seismo.shrink_to_fit();
}


template <typename ValueType>
void KITGPI::Acquisition::SeismogramHandler<ValueType>::resetData()
{
    for(auto i : seismo){
        i.resetData();
    }
}

template <typename ValueType>
IndexType KITGPI::Acquisition::SeismogramHandler<ValueType>::getNumSamples(SeismogramType type) const
{
    SCAI_ASSERT_ERROR(type >= 0 && type <= NUM_ELEMENTS_ENUM-1, "SeismogramType unkown");
    return(seismo[type].getNumSamples());
}

template <typename ValueType>
IndexType KITGPI::Acquisition::SeismogramHandler<ValueType>::getNumTracesGlobal(SeismogramType type) const
{
    SCAI_ASSERT_ERROR(type >= 0 && type <= NUM_ELEMENTS_ENUM-1, "SeismogramType unkown");
    return(seismo[type].getNumTracesGlobal());
}

template <typename ValueType>
KITGPI::Acquisition::Seismogram<ValueType>& KITGPI::Acquisition::SeismogramHandler<ValueType>::getSeismogram(IndexType type)
{
    SCAI_ASSERT_ERROR(type >= 0 && type <= NUM_ELEMENTS_ENUM-1, "SeismogramType unkown");
    return(seismo[type]);
}

template <typename ValueType>
KITGPI::Acquisition::Seismogram<ValueType>const& KITGPI::Acquisition::SeismogramHandler<ValueType>::getSeismogram(SeismogramType type) const
{
    SCAI_ASSERT_ERROR(type >= 0 && type <= NUM_ELEMENTS_ENUM-1, "SeismogramType unkown");
    return(seismo[type]);
}

template <typename ValueType>
KITGPI::Acquisition::Seismogram<ValueType>& KITGPI::Acquisition::SeismogramHandler<ValueType>::getSeismogram(SeismogramType type)
{
    SCAI_ASSERT_ERROR(type >= 0 && type <= NUM_ELEMENTS_ENUM-1, "SeismogramType unkown");
    return(seismo[type]);
}

template <typename ValueType>
void KITGPI::Acquisition::SeismogramHandler<ValueType>::setContextPtr(hmemo::ContextPtr ctx)
{
    for(auto i : seismo){
        i.setContextPtr(ctx);
    }
}

template <typename ValueType>
void KITGPI::Acquisition::SeismogramHandler<ValueType>::setDT(ValueType newDT)
{
    for(auto i : seismo){
        i.setDT(newDT);
    }
}

template <typename ValueType>
void KITGPI::Acquisition::SeismogramHandler<ValueType>::setSourceCoordinate(IndexType sourceCoord)
{
    for(auto i : seismo){
        i.setSourceCoordinate(sourceCoord);
    }
}
