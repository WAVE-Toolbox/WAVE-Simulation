#include "WavefieldsFactory.hpp"

using namespace scai;

template <typename ValueType>
typename KITGPI::Wavefields::WavefieldsEM<ValueType>::WavefieldPtr KITGPI::Wavefields::FactoryEM<ValueType>::Create(std::string dimension, std::string type)
{

    // transform to lower cases
    std::transform(type.begin(), type.end(), type.begin(), ::tolower);
    std::transform(dimension.begin(), dimension.end(), dimension.begin(), ::tolower);

    // Assert correctness of input values
    SCAI_ASSERT_ERROR(dimension.compare("2d") == 0 || dimension.compare("3d") == 0, "Unkown dimension");
    SCAI_ASSERT_ERROR(type.compare("acoustic") == 0 || type.compare("elastic") == 0 || type.compare("viscoelastic") == 0 || type.compare("sh") == 0 || type.compare("emem") == 0 || type.compare("tmem") == 0 || type.compare("viscoemem") == 0 || type.compare("viscotmem") == 0, "Unkown type");

    // 2D
    if (dimension.compare("2d") == 0 && type.compare("tmem") == 0) {
        return WavefieldPtrEM(new FD2Dtmem<ValueType>);
    }
    if (dimension.compare("2d") == 0 && type.compare("emem") == 0) {
        return WavefieldPtrEM(new FD2Demem<ValueType>);
    }
    if (dimension.compare("2d") == 0 && type.compare("viscotmem") == 0) {
        return WavefieldPtrEM(new FD2Dviscotmem<ValueType>);
    }
    if (dimension.compare("2d") == 0 && type.compare("viscoemem") == 0) {
        return WavefieldPtrEM(new FD2Dviscoemem<ValueType>);
    }

    // 3D
    if (dimension.compare("3d") == 0 && type.compare("emem") == 0) {
        return WavefieldPtrEM(new FD3Demem<ValueType>);
    }
    if (dimension.compare("3d") == 0 && type.compare("viscoemem") == 0) {
        return WavefieldPtrEM(new FD3Dviscoemem<ValueType>);
    }

    return WavefieldPtrEM();
}

template class KITGPI::Wavefields::FactoryEM<double>;
template class KITGPI::Wavefields::FactoryEM<float>;
