#include "WavefieldsFactory.hpp"

using namespace scai;

template <typename ValueType>
typename KITGPI::Wavefields::Wavefields<ValueType>::WavefieldPtr KITGPI::Wavefields::Factory<ValueType>::Create(std::string dimension, std::string type)
{

    // transform to lower cases
    std::transform(type.begin(), type.end(), type.begin(), ::tolower);
    std::transform(dimension.begin(), dimension.end(), dimension.begin(), ::tolower);

    // Assert correctness of input values
    SCAI_ASSERT_ERROR(dimension.compare("2d") == 0 || dimension.compare("3d") == 0, "Unkown dimension");
    SCAI_ASSERT_ERROR(type.compare("acoustic") == 0 || type.compare("elastic") == 0 || type.compare("viscoelastic") == 0 || type.compare("sh") == 0 || type.compare("emem") == 0 || type.compare("tmem") == 0 || type.compare("viscoemem") == 0 || type.compare("viscotmem") == 0, "Unkown type");

    // 2D
    if (dimension.compare("2d") == 0 && type.compare("acoustic") == 0) {
        return WavefieldPtr(new FD2Dacoustic<ValueType>);
    }
    if (dimension.compare("2d") == 0 && type.compare("elastic") == 0) {
        return WavefieldPtr(new FD2Delastic<ValueType>);
    }
    if (dimension.compare("2d") == 0 && type.compare("viscoelastic") == 0) {
        return WavefieldPtr(new FD2Dviscoelastic<ValueType>);
    }
    if (dimension.compare("2d") == 0 && type.compare("sh") == 0) {
        return WavefieldPtr(new FD2Dsh<ValueType>);
    }
    if (dimension.compare("2d") == 0 && type.compare("viscosh") == 0) {
        return WavefieldPtr(new FD2Dviscosh<ValueType>);
    }

    // 3D
    if (dimension.compare("3d") == 0 && type.compare("acoustic") == 0) {
        return WavefieldPtr(new FD3Dacoustic<ValueType>);
    }
    if (dimension.compare("3d") == 0 && type.compare("elastic") == 0) {
        return WavefieldPtr(new FD3Delastic<ValueType>);
    }
    if (dimension.compare("3d") == 0 && type.compare("viscoelastic") == 0) {
        return WavefieldPtr(new FD3Dviscoelastic<ValueType>);
    }

    // 2D
    if (dimension.compare("2d") == 0 && type.compare("tmem") == 0) {
        return WavefieldPtr(new FD2Dtmem<ValueType>);
    }
    if (dimension.compare("2d") == 0 && type.compare("emem") == 0) {
        return WavefieldPtr(new FD2Demem<ValueType>);
    }
    if (dimension.compare("2d") == 0 && type.compare("viscotmem") == 0) {
        return WavefieldPtr(new FD2Dviscotmem<ValueType>);
    }
    if (dimension.compare("2d") == 0 && type.compare("viscoemem") == 0) {
        return WavefieldPtr(new FD2Dviscoemem<ValueType>);
    }

    // 3D
    if (dimension.compare("3d") == 0 && type.compare("emem") == 0) {
        return WavefieldPtr(new FD3Demem<ValueType>);
    }
    if (dimension.compare("3d") == 0 && type.compare("viscoemem") == 0) {
        return WavefieldPtr(new FD3Dviscoemem<ValueType>);
    }
    
    return WavefieldPtr();
}

template class KITGPI::Wavefields::Factory<double>;
template class KITGPI::Wavefields::Factory<float>;
