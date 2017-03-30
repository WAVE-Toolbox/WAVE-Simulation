#include "Wavefields.hpp"

using namespace scai;

/*! \brief Reset a single wavefield to zero.
 \param vector Vector to be reset to 0
 */
template <typename ValueType>
void KITGPI::Wavefields::Wavefields<ValueType>::resetWavefield(scai::lama::DenseVector<ValueType> &vector)
{
    vector.assign(0.0);
}

/*! \brief Intitialisation of a single wavefield vector.
 *
 * This method will set the context, allocate the the wavefield and set the field to zero.
 *
 \param vector Vector to be set
 \param ctx Context pointer
 \param dist Distribution
 */
template <typename ValueType>
void KITGPI::Wavefields::Wavefields<ValueType>::initWavefield(scai::lama::DenseVector<ValueType> &vector, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist)
{
    vector.setContextPtr(ctx);
    vector.allocate(dist);

    resetWavefield(vector);
}

/*! \brief Methode to Write Wavefield for timestep t
 *
 \param vector Vector written to file
 \param type Wavefield-type (acoustic, elastic, viscoelastic)
 \param t Timestep
 */
template <typename ValueType>
void KITGPI::Wavefields::Wavefields<ValueType>::writeWavefield(scai::lama::DenseVector<ValueType> &vector, std::string vectorName, std::string type, IndexType t)
{
    std::string fileName = "wavefields/wavefield" + type + "." + vectorName + "." + std::to_string(t) + ".mtx";
    std::cout << "snapshot for Timestep " << t << "has been written to: " << fileName;

    vector.writeToFile(fileName);
}

//! \brief Getter routine for vX wavefield
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::Wavefields<ValueType>::getVX()
{
    return (VX);
}

//! \brief Getter routine for vY wavefield
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::Wavefields<ValueType>::getVY()
{
    return (VY);
}

//! \brief Getter routine for vZ wavefield
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::Wavefields<ValueType>::getVZ()
{
    return (VZ);
}

//! \brief Getter routine for Sxx wavefield
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::Wavefields<ValueType>::getSxx()
{
    return (Sxx);
}

//! \brief Getter routine for Syy wavefield
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::Wavefields<ValueType>::getSyy()
{
    return (Syy);
}

//! \brief Getter routine for Szz wavefield
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::Wavefields<ValueType>::getSzz()
{
    return (Szz);
}

//! \brief Getter routine for Syz wavefield
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::Wavefields<ValueType>::getSyz()
{
    return (Syz);
}

//! \brief Getter routine for Sxz wavefield
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::Wavefields<ValueType>::getSxz()
{
    return (Sxz);
}

//! \brief Getter routine for Sxy wavefield
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::Wavefields<ValueType>::getSxy()
{
    return (Sxy);
}

//! \brief Getter routine for P
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::Wavefields<ValueType>::getP()
{
    return (P);
}

//! \brief Getter routine for Rxx Relaxation parameter
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::Wavefields<ValueType>::getRxx()
{
    return (Rxx);
}

//! \brief Getter routine for Ryy Relaxation parameter
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::Wavefields<ValueType>::getRyy()
{
    return (Ryy);
}

//! \brief Getter routine for Rzz Relaxation parameter
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::Wavefields<ValueType>::getRzz()
{
    return (Rzz);
}

//! \brief Getter routine for Ryz Relaxation parameter
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::Wavefields<ValueType>::getRyz()
{
    return (Ryz);
}

//! \brief Getter routine for Rxz Relaxation parameter
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::Wavefields<ValueType>::getRxz()
{
    return (Rxz);
}

//! \brief Getter routine for Rxy Relaxation parameter
template <typename ValueType>
scai::lama::DenseVector<ValueType> &KITGPI::Wavefields::Wavefields<ValueType>::getRxy()
{
    return (Rxy);
}

template class KITGPI::Wavefields::Wavefields<float>;
template class KITGPI::Wavefields::Wavefields<double>;
