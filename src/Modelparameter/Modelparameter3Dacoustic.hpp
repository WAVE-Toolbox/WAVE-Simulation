//! Class for Modelparameter for 3-D acoustic simulations (Subsurface properties)
/*!
 This class handels the modelparameter for the 3-D acoustic finite-difference simulation.
 */


#pragma once

#include <scai/lama.hpp>

#include <scai/lama/DenseVector.hpp>
#include <scai/lama/expression/all.hpp>
#include <scai/lama/matrix/all.hpp>
#include <scai/lama/matutils/MatrixCreator.hpp>
#include <scai/lama/norm/L2Norm.hpp>

#include <scai/dmemo/BlockDistribution.hpp>

#include <scai/hmemo/ReadAccess.hpp>
#include <scai/hmemo/WriteAccess.hpp>
#include <scai/hmemo/HArray.hpp>

#include <scai/tracing.hpp>

#include <scai/common/Walltime.hpp>
#include <scai/common/unique_ptr.hpp>
#include <scai/logging.hpp>

#include <iostream>

#include "Modelparameter.hpp"
#include "../Configuration.hpp"

template<typename ValueType>
class Modelparameter3Dacoustic : private Modelparameter<ValueType>
{
public:
    
    //! Default constructor.
    Modelparameter3Dacoustic(){};

    Modelparameter3Dacoustic(Configuration<ValueType> config, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist);
    Modelparameter3Dacoustic(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, lama::Scalar  M, lama::Scalar  rho);
    Modelparameter3Dacoustic(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filenamepi, std::string filenamerho);
    Modelparameter3Dacoustic(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filename);
    
    //! Copy Constructor.
    Modelparameter3Dacoustic(const Modelparameter3Dacoustic& rhs);
    
    //! Destructor, releases all allocated resources.
    ~Modelparameter3Dacoustic(){};

    void init(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, lama::Scalar  M, lama::Scalar  rho);
    void init(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filename);
    void init(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filenamepi, std::string filenamerho);
    
    void write(std::string filenamepi, std::string filenamedensity);
    void write(std::string filename);
        
    lama::DenseVector<ValueType> pi; //!< Vector storing first Lame-Parameter.
    lama::DenseVector<ValueType> density; //!< Vector storing Density.
    
};

/*! \brief Constructor that is using the Configuration class
 *
 \param config Configuration class
 \param ctx Context for the Calculation
 \param dist Distribution
 */
template<typename ValueType>
Modelparameter3Dacoustic<ValueType>::Modelparameter3Dacoustic(Configuration<ValueType> config, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist)
{
    if(config.getReadModel()){
        init(ctx,dist,config.getFilenameModel());
    } else {
        init(ctx,dist,config.getM(),config.getRho());
    }
    
    write(config.getFilenameModel()+".out");

}

/*! \brief Constructor that is generating a homogeneous model
 *
 *  Generates a homogeneous model, which will be initialized by the two given scalar values.
 \param ctx Context
 \param dist Distribution
 \param M First Lame-Parameter given as Scalar
 \param rho Density given as Scalar
 */
template<typename ValueType>
Modelparameter3Dacoustic<ValueType>::Modelparameter3Dacoustic(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, lama::Scalar  M, lama::Scalar  rho)
{
    init(ctx,dist,M,rho);
}


/*! \brief Initialisation that is generating a homogeneous model
 *
 *  Generates a homogeneous model, which will be initialized by the two given scalar values.
 \param ctx Context
 \param dist Distribution
 \param M First Lame-Parameter given as Scalar
 \param rho Density given as Scalar
 */
template<typename ValueType>
void Modelparameter3Dacoustic<ValueType>::init(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, lama::Scalar  M, lama::Scalar  rho)
{
    this->initModelparameter(pi,ctx,dist,M);
    this->initModelparameter(density,ctx,dist,rho);
}


/*! \brief Constructor that is reading models from external files
 *
 *  Reads a model from an external file.
 \param ctx Context
 \param dist Distribution
 \param filenamepi Name of file that will be read for the first Lame-parameter.
 \param filenamerho Name of file that will be read for the Density.
 */
template<typename ValueType>
Modelparameter3Dacoustic<ValueType>::Modelparameter3Dacoustic(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filenamepi, std::string filenamerho)
{
    init(ctx,dist,filenamepi,filenamerho);
}


/*! \brief Initialisation that is reading models from external files
 *
 *  Reads a model from an external file.
 \param ctx Context
 \param dist Distribution
 \param filenamepi Name of file that will be read for the first Lame-parameter.
 \param filenamerho Name of file that will be read for the Density.
 */
template<typename ValueType>
void Modelparameter3Dacoustic<ValueType>::init(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filenamepi, std::string filenamerho)
{
    this->initModelparameter(pi,ctx,dist,filenamepi);
    this->initModelparameter(density,ctx,dist,filenamerho);
}


/*! \brief Constructor that is reading models from external files
 *
 *  Reads a model from an external file.
 \param ctx Context
 \param dist Distribution
 \param filename For the first Lame-parameter ".pi.mtx" is added and for density ".density.mtx" is added.
 */
template<typename ValueType>
Modelparameter3Dacoustic<ValueType>::Modelparameter3Dacoustic(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filename)
{
    init(ctx,dist,filename);
}


/*! \brief Initialisator that is reading models from external files
 *
 *  Reads a model from an external file.
 \param ctx Context
 \param dist Distribution
 \param filename For the first Lame-parameter "model/"+filename+".pi.mtx" is added and for density "model/"+filename+".density.mtx" is added.
 */
template<typename ValueType>
void Modelparameter3Dacoustic<ValueType>::init(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filename)
{
    std::string filenamepi=filename+".pi.mtx";
    std::string filenamedensity=filename+".density.mtx";

    this->initModelparameter(pi,ctx,dist,filenamepi);
    this->initModelparameter(density,ctx,dist,filenamedensity);
}


//! \brief Copy constructor
template<typename ValueType>
Modelparameter3Dacoustic<ValueType>::Modelparameter3Dacoustic(const Modelparameter3Dacoustic& rhs)
{
    pi=rhs.pi.copy();
    density=rhs.density.copy();
}


/*! \brief Write model to an external file
 *
 \param filenamepi Filename for first Lame-Parameter model
 \param filenamedensity Filename for Density model
 */
template<typename ValueType>
void Modelparameter3Dacoustic<ValueType>::write( std::string filenamepi, std::string filenamedensity)
{
    this->writeModelparameter(pi,filenamepi);
    this->writeModelparameter(density,filenamedensity);
};


/*! \brief Write model to an external file
 *
 \param filename Filename to write files. For the first Lame-parameter ".pi.mtx" is added and for density ".density.mtx" is added.
 */
template<typename ValueType>
void Modelparameter3Dacoustic<ValueType>::write(std::string filename)
{
    std::string filenamepi=filename+".pi.mtx";
    std::string filenamedensity=filename+".density.mtx";
    this->writeModelparameter(pi,filenamepi);
    this->writeModelparameter(density,filenamedensity);
};





