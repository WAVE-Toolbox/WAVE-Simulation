#pragma once

#include "../Configuration/Configuration.hpp"
#include "../Modelparameter/Modelparameter.hpp"
#include "../Acquisition/Sources.hpp"

namespace KITGPI
{
    //! \brief CheckParameter namespace
    namespace CheckParameter
    {
	template<typename ValueType, typename IndexType> void checkNumericalArtefeactsAndInstabilities(Configuration::Configuration const &config, Modelparameter::Modelparameter<ValueType> &model,scai::dmemo::CommunicatorPtr comm); ;
	template<typename ValueType, typename IndexType> void checkStabilityCriterion(ValueType dt, ValueType dh, ValueType vpMax, std::string dimension, IndexType FDorder, scai::dmemo::CommunicatorPtr comm);
	template<typename ValueType, typename IndexType> void checkNumericalDispersion(ValueType dh, ValueType vsMin, ValueType fc, IndexType spFDo, scai::dmemo::CommunicatorPtr comm);
	
	template<typename ValueType, typename IndexType> void checkAcquisitionGeometry(Configuration::Configuration const &config,scai::dmemo::CommunicatorPtr comm);
	template<typename ValueType, typename IndexType> void checkSources(IndexType NX, IndexType NY, IndexType NZ, std::string sourcefile, scai::dmemo::CommunicatorPtr comm);
	template<typename ValueType, typename IndexType> void checkReceivers(IndexType NX, IndexType NY, IndexType NZ, std::string receiverfile, scai::dmemo::CommunicatorPtr comm);
    }
}

# include "./CheckParameter.tpp"