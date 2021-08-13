
#pragma once

#include "SourceReceiverImpl.hpp"

#include "FDTD2Demem.hpp"
#include "FDTD2Dtmem.hpp"

#include "FDTD3Demem.hpp"

namespace KITGPI
{

    namespace ForwardSolver
    {
        
        //! \brief SourceReceiverImpl namespace
        namespace SourceReceiverImpl
        {

            //! \brief Factory class.
            template <typename ValueType>
            class FactoryEM
            {
            public:
                //! \brief Declare SourceReceiverImpl pointer
                typedef typename SourceReceiverImplEM<ValueType>::SourceReceiverImplPtr SourceReceiverImplPtr;

                FactoryEM() = delete;
                FactoryEM(FactoryEM const &) = delete;
                void operator=(FactoryEM const &) = delete;

                /*! \brief Create the right simmulation with factory methode.
                *
                \param dimension Dimension of the model (2D, 3D)
                \param type Simmulation type (acoustic, elsstic, viscoelastic)
                */
                static SourceReceiverImplPtr Create(std::string dimension, std::string type, Acquisition::AcquisitionGeometry<ValueType> const &sourceConfig, Acquisition::AcquisitionGeometry<ValueType> &receiverConfig, Wavefields::WavefieldsEM<ValueType> &wavefieldIN);
            };
        }/* end namespace SourceReceiverImpl */
    } /* end namespace ForwardSolver */
} /* end namespace KITGPI */
