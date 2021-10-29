#include "SourceReceiverImplFactory.hpp"

template <typename ValueType>
typename KITGPI::ForwardSolver::SourceReceiverImpl::SourceReceiverImpl<ValueType>::SourceReceiverImplPtr KITGPI::ForwardSolver::SourceReceiverImpl::Factory<ValueType>::Create(std::string dimension, std::string type, Acquisition::AcquisitionGeometry<ValueType> const &sourceConfig, Acquisition::AcquisitionGeometry<ValueType> &receiverConfig, Wavefields::Wavefields<ValueType> &wavefieldIN)
{
    // transform to lower cases
    std::transform(dimension.begin(), dimension.end(), dimension.begin(), ::tolower);
    std::transform(type.begin(), type.end(), type.begin(), ::tolower);

    // Assert correctness of input values
    SCAI_ASSERT_ERROR(dimension.compare("2d") == 0 || dimension.compare("3d") == 0, "Unkown dimension");
    SCAI_ASSERT_ERROR(type.compare("acoustic") == 0 || type.compare("elastic") == 0 || type.compare("viscoelastic") == 0 || type.compare("sh") == 0 || type.compare("viscosh") == 0 || type.compare("emem") == 0 || type.compare("tmem") == 0 || type.compare("viscoemem") == 0 || type.compare("viscotmem") == 0, "Unkown type");

    // 2D
    if (dimension.compare("2d") == 0 && type.compare("acoustic") == 0) {
        return SourceReceiverImplPtr(new FDTD2Dacoustic<ValueType>(sourceConfig, receiverConfig, wavefieldIN));
    }
    if (dimension.compare("2d") == 0 && (type.compare("elastic") == 0 || type.compare("viscoelastic") == 0)) {
        return SourceReceiverImplPtr(new FDTD2Delastic<ValueType>(sourceConfig, receiverConfig, wavefieldIN));
    }
    if (dimension.compare("2d") == 0 && (type.compare("sh") == 0 || type.compare("viscosh") == 0)) {
        return SourceReceiverImplPtr(new FDTD2Dsh<ValueType>(sourceConfig, receiverConfig, wavefieldIN));
    }

    // 3D
    if (dimension.compare("3d") == 0 && type.compare("acoustic") == 0) {
        return SourceReceiverImplPtr(new FDTD3Dacoustic<ValueType>(sourceConfig, receiverConfig, wavefieldIN));
    }
    if (dimension.compare("3d") == 0 && (type.compare("elastic") == 0 || type.compare("viscoelastic") == 0)) {
        return SourceReceiverImplPtr(new FDTD3Delastic<ValueType>(sourceConfig, receiverConfig, wavefieldIN));
    }

    if (dimension.compare("2d") == 0 && (type.compare("emem") == 0 || type.compare("viscoemem") == 0)) {
        return SourceReceiverImplPtr(new FDTD2Demem<ValueType>(sourceConfig, receiverConfig, wavefieldIN));
    }
    if (dimension.compare("2d") == 0 && (type.compare("tmem") == 0 || type.compare("viscotmem") == 0)) {
        return SourceReceiverImplPtr(new FDTD2Dtmem<ValueType>(sourceConfig, receiverConfig, wavefieldIN));
    }

    // 3D
    if (dimension.compare("3d") == 0 && (type.compare("emem") == 0 || type.compare("viscoemem") == 0)) {
        return SourceReceiverImplPtr(new FDTD3Demem<ValueType>(sourceConfig, receiverConfig, wavefieldIN));
    }
    
    return SourceReceiverImplPtr();
};

template class KITGPI::ForwardSolver::SourceReceiverImpl::Factory<double>;
template class KITGPI::ForwardSolver::SourceReceiverImpl::Factory<float>;
