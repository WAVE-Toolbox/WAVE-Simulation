#include "ModelparameterFactory.hpp"

template <typename ValueType>
typename KITGPI::Modelparameter::Modelparameter<ValueType>::ModelparameterPtr KITGPI::Modelparameter::Factory<ValueType>::Create(std::string type)
{

    // transform to lower cases
    std::transform(type.begin(), type.end(), type.begin(), ::tolower);

    // Assert correctness of input values
    SCAI_ASSERT_ERROR(type.compare("acoustic") == 0 || type.compare("elastic") == 0 || type.compare("viscoelastic") == 0 || type.compare("sh") == 0 || type.compare("emem") == 0 || type.compare("tmem") == 0 || type.compare("viscoemem") == 0 || type.compare("viscotmem") == 0, "Unkown type");
    
    if (type.compare("acoustic") == 0) {
        return ModelparameterPtr(new Acoustic<ValueType>);
    }
    if (type.compare("elastic") == 0) {
        return ModelparameterPtr(new Elastic<ValueType>);
    }
    if (type.compare("viscoelastic") == 0) {
        return ModelparameterPtr(new Viscoelastic<ValueType>);
    }
    if (type.compare("sh") == 0) {
        return ModelparameterPtr(new SH<ValueType>);
    }
    if (type.compare("viscosh") == 0) {
        return ModelparameterPtr(new ViscoSH<ValueType>);
    }
    
    if (type.compare("tmem") == 0) {
        return ModelparameterPtr(new TMEM<ValueType>);
    }
    if (type.compare("emem") == 0) {
        return ModelparameterPtr(new EMEM<ValueType>);
    }
    if (type.compare("viscotmem") == 0) {
        return ModelparameterPtr(new ViscoTMEM<ValueType>);
    }
    if (type.compare("viscoemem") == 0) {
        return ModelparameterPtr(new ViscoEMEM<ValueType>);
    }

    return ModelparameterPtr();
}

template class KITGPI::Modelparameter::Factory<float>;
template class KITGPI::Modelparameter::Factory<double>;
