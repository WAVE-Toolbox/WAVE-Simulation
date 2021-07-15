
#pragma once

#include "SourceReceiverImpl.hpp"

#include "FDTD2Dacoustic.hpp"
#include "FDTD2Delastic.hpp"
#include "FDTD2Dsh.hpp"

#include "FDTD3Dacoustic.hpp"
#include "FDTD3Delastic.hpp"

namespace KITGPI
{

    namespace ForwardSolver
    {
        
        //! \brief SourceReceiverImpl namespace
        namespace SourceReceiverImpl
        {

            //! \brief Factory class.
            template <typename ValueType>
            class Factory
            {
            public:
                //! \brief Declare SourceReceiverImpl pointer
                typedef typename SourceReceiverImpl<ValueType>::SourceReceiverImplPtr SourceReceiverImplPtr;

                Factory() = delete;
                Factory(Factory const &) = delete;
                void operator=(Factory const &) = delete;

                /*! \brief Create the right simmulation with factory methode.
                *
                \param dimension Dimension of the model (2D, 3D)
                \param type Simmulation type (acoustic, elsstic, viscoelastic)
                */
                static SourceReceiverImplPtr Create(std::string dimension, std::string type, Acquisition::AcquisitionGeometry<ValueType> const &sourceConfig, Acquisition::AcquisitionGeometry<ValueType> &receiverConfig, Wavefields::Wavefields<ValueType> &wavefieldIN);
            };
        }/* end namespace SourceReceiverImpl */
    } /* end namespace ForwardSolver */
} /* end namespace KITGPI */
