//! Class for Modelparameters (Subsurface properties)
/*!
 This class handels the modelparameters for the finite-difference simulation.
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


template<typename ValueType>
class Modelparameter
{
public:
    
    //! Default constructor.
    Modelparameter(){};
    Modelparameter(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, lama::Scalar  M, lama::Scalar  rho);
    Modelparameter(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filename);
    //! Copy Constructor.
    Modelparameter(const Modelparameter& rhs);
    
    //! Destructor, releases all allocated resources.
    ~Modelparameter(){};

    void init(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, lama::Scalar  M, lama::Scalar  rho);
    void init(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filename);
    
    void write(std::string filename);
    
    /// Vector storing first Lame-Parameter.
    lama::DenseVector<ValueType> pi;
    /// Vector storing Density.
    lama::DenseVector<ValueType> density;
    
private:
    
    void allocate(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist);
    void read(std::string filename);
    
};


/*! \brief Constructor that is generating a homogeneous model
 *
 *  Generates a homogeneous model, which will be initialized by the two given scalar values.
 \param M First Lame-Parameter given as Scalar
 \param rho Density given as Scalar
 */
template<typename ValueType>
Modelparameter<ValueType>::Modelparameter(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, lama::Scalar  M, lama::Scalar  rho) {
    init(ctx,dist,M,rho);
}

/*! \brief Constructor that is reading a model from an external file
 *
 *  Reads a model from an external mtx file. The extension _pi.mtx will be added to read in the
 * first Lame Parameter and the extension _density.mtx will be added to read in the density.
 \param filename Name of file that will be read.
 */
template<typename ValueType>
Modelparameter<ValueType>::Modelparameter(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filename) {
    init(ctx,dist,filename);
}

//! \brief Copy constructor
template<typename ValueType>
Modelparameter<ValueType>::Modelparameter(const Modelparameter& rhs){
    pi=rhs.pi.copy();
    density=rhs.density.copy();
}

/*! \brief Init model by a homogeneous model
 *
 *  Generates a homogeneous model, which will be initialized by the two given scalar values.
 \param M First Lame-Parameter given as Scalar
 \param rho Density given as Scalar
 */
template<typename ValueType>
void Modelparameter<ValueType>::init(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, lama::Scalar  M, lama::Scalar  rho) {
    
    allocate(ctx,dist);
    
    pi.assign(M);
    density.assign(rho);
    
    std::string filename="model/out";
    write(filename);
}

/*! \brief Init model by reading a model from an external file
 *
 *  Reads a model from an external mtx file. The extension _pi.mtx will be added to read in the
 * first Lame Parameter and the extension _density.mtx will be added to read in the density.
 \param filename Name of file that will be read.
 */
template<typename ValueType>
void Modelparameter<ValueType>::init(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filename) {
    
    allocate(ctx,dist);
    
    read(filename);
    
    pi.redistribute(dist);
    density.redistribute(dist);
    
    filename=filename+"_out";
    write(filename);
}

/*! \brief Write model to an external file
 *
 *  Write a model to an external mtx file. The extension _pi.mtx will be added to write the
 * first Lame Parameter and the extension _density.mtx will be added to write the density.
 \param filename Name of output file.
 */
template<typename ValueType>
void Modelparameter<ValueType>::write(std::string filename){
    std::string filename_pi=filename+"_pi.mtx";
    std::string filename_density=filename+"_density.mtx";
    
    pi.writeToFile(filename_pi);
    density.writeToFile(filename_density);
};

template<typename ValueType>
void Modelparameter<ValueType>::allocate(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist) {
    pi.setContextPtr(ctx);
    density.setContextPtr(ctx);
    pi.allocate(dist);
    density.allocate(dist);
};

template<typename ValueType>
void Modelparameter<ValueType>::read(std::string filename){
    std::string filename_pi=filename+"_pi.mtx";
    std::string filename_density=filename+"_density.mtx";
    
    pi.readFromFile(filename_pi);
    density.readFromFile(filename_density);
};




