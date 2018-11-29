#pragma once

#include "../Acquisition/Receivers.hpp"
#include "../Acquisition/Sources.hpp"
#include "../Acquisition/suHandler.hpp"
#include "../Common/Common.hpp"
#include "../Configuration/Configuration.hpp"
#include "../Modelparameter/Modelparameter.hpp"
#include <scai/lama.hpp>

namespace KITGPI
{
    //! \brief CheckParameter namespace
    namespace CheckParameter
    {
        template <typename ValueType>
        void checkNumericalArtefeactsAndInstabilities(Configuration::Configuration const &config, Modelparameter::Modelparameter<ValueType> &model, scai::dmemo::CommunicatorPtr comm);
        template <typename ValueType>
        void checkStabilityCriterion(ValueType dt, ValueType dh, ValueType vpMax, std::string dimension, scai::IndexType FDorder, scai::dmemo::CommunicatorPtr comm);
        template <typename ValueType>
        void checkNumericalDispersion(ValueType dh, ValueType vsMin, ValueType fc, scai::IndexType spFDo, scai::dmemo::CommunicatorPtr comm);
        template <typename ValueType>
        void checkSources(Configuration::Configuration const &config, Acquisition::Sources<ValueType> const &sources, scai::dmemo::CommunicatorPtr comm);
        template <typename ValueType>
        void checkReceivers(Configuration::Configuration const &config, Acquisition::Receivers<ValueType> const &receiver, scai::dmemo::CommunicatorPtr comm);
    }
}

#include "./CheckParameter.tpp"