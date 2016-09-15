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
    
    Modelparameter(){};
    Modelparameter(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, lama::Scalar  M, lama::Scalar  rho);
    Modelparameter(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filename);
    
    ~Modelparameter(){}
    
    void init(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, lama::Scalar  M, lama::Scalar  rho);
    void init(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filename);
    
    lama::DenseVector<ValueType> pi;
    lama::DenseVector<ValueType> density;
    
private:
    
    void allocate(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist);
    void write(std::string filename);
    void read(std::string filename);
    
};


template<typename ValueType>
Modelparameter<ValueType>::Modelparameter(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, lama::Scalar  M, lama::Scalar  rho) {
    init(ctx,dist,M,rho);
}


template<typename ValueType>
Modelparameter<ValueType>::Modelparameter(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filename) {
    init(ctx,dist,filename);
}

template<typename ValueType>
void Modelparameter<ValueType>::init(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, lama::Scalar  M, lama::Scalar  rho) {
    
    allocate(ctx,dist);
    
    pi.assign(M);
    density.assign(rho);
    
    write("model/test");
}

template<typename ValueType>
void Modelparameter<ValueType>::init(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filename) {
    
    allocate(ctx,dist);
    
    read(filename);
    
    pi.redistribute(dist);
    density.redistribute(dist);
    
    filename=filename+"_out";
    write(filename);
}

template<typename ValueType>
void Modelparameter<ValueType>::allocate(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist) {
    pi.setContextPtr(ctx);
    density.setContextPtr(ctx);
    pi.allocate(dist);
    density.allocate(dist);
};

template<typename ValueType>
void Modelparameter<ValueType>::write(std::string filename){
    std::string filename_pi=filename+"_pi.mtx";
    std::string filename_density=filename+"_density.mtx";
    
    pi.writeToFile(filename_pi);
    density.writeToFile(filename_density);
};

template<typename ValueType>
void Modelparameter<ValueType>::read(std::string filename){
    std::string filename_pi=filename+"_pi.mtx";
    std::string filename_density=filename+"_density.mtx";
    
    pi.readFromFile(filename_pi);
    density.readFromFile(filename_density);
};