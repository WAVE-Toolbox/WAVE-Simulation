#include "ModelparameterFactory.hpp"

template <typename ValueType>
typename KITGPI::Modelparameter::Modelparameter<ValueType>::ModelparameterPtr KITGPI::Modelparameter::Factory<ValueType>::Create(std::string type)
{

    // transform to lower cases
    std::transform(type.begin(), type.end(), type.begin(), ::tolower);

    // Assert correctness of input values
    SCAI_ASSERT_ERROR(type.compare("acoustic") == 0 || type.compare("elastic") == 0 || type.compare("visco") == 0 || type.compare("sh") == 0, "Unkown type");

    if (type.compare("acoustic") == 0) {
        return ModelparameterPtr(new Acoustic<ValueType>);
    }
    if (type.compare("elastic") == 0) {
        return ModelparameterPtr(new Elastic<ValueType>);
    }
    if (type.compare("visco") == 0) {
        return ModelparameterPtr(new Viscoelastic<ValueType>);
    }
    if (type.compare("sh") == 0) {
        return ModelparameterPtr(new SH<ValueType>);
    }

    COMMON_THROWEXCEPTION("Reached end of factory without match");
    return nullptr;
}

template class KITGPI::Modelparameter::Factory<float>;
template class KITGPI::Modelparameter::Factory<double>;
