//! Abstract class for a single Modelparameter (Subsurface properties)
/*!
 * This class handels a single modelparameter.
 * As this class is an abstract class, all constructors are protected.
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
    
    //! Default destructor.
    ~Modelparameter(){};

    /*! \brief Abstract initialisation function 
     *
     * Standard initialisation function
     *
     \param ctx Context
     \param dist Distribution
     \param filename filename to read modelparameters (endings will be added by derived classes)
     */
    virtual void init(hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filename)=0;
    
    /*! \brief Abstract write function
     *
     * Standard write function
     *
     \param filename filename to write modelparameters (endings will be added by derived classes)
     */
    virtual void write(std::string filename)=0;
    
    //! \brief Get referenec to density model parameter
    virtual lama::DenseVector<ValueType>* getDensity()=0;
    //! \brief Get referenec to density model parameter
    virtual lama::DenseVector<ValueType>* getInverseDensity()=0;
    //! \brief Get referenec to first Lame model parameter
    virtual lama::DenseVector<ValueType>* getM()=0;
    //! \brief Get referenec to P-wave velocity
    virtual lama::DenseVector<ValueType>* getVelocityP()=0;
    //! \brief Get referenec to S-wave velocity
    virtual lama::DenseVector<ValueType>* getVelocityS()=0;
    
protected:
    void initModelparameter(lama::DenseVector<ValueType>& vector, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, lama::Scalar  value);
    void initModelparameter(lama::DenseVector<ValueType>& vector, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filename);
    
    void writeModelparameter(lama::DenseVector<ValueType>& vector, std::string filename);
    
    
private:
    void allocateModelparameter(lama::DenseVector<ValueType>& vector, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist);

    void readModelparameter(lama::DenseVector<ValueType>& vector, std::string filename);
};


/*! \brief Init a single modelparameter by a constant value
 *
 \param vector Singel modelparameter which will be initialized
 \param ctx Context
 \param dist Distribution
 \param value Value which will be used to initialize the single modelparameter to a homogenoeus model
 */
template<typename ValueType>
void Modelparameter<ValueType>::initModelparameter(lama::DenseVector<ValueType>& vector, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, lama::Scalar  value)
{
    
    allocateModelparameter(vector,ctx,dist);
    
    vector.assign(value);
    
}


/*! \brief Init a single modelparameter by reading a model from an external file
 *
 *  Reads a single model from an external mtx file.
 \param vector Singel modelparameter which will be initialized
 \param ctx Context
 \param dist Distribution
 \param filename Location of external file which will be read in 
 */
template<typename ValueType>
void Modelparameter<ValueType>::initModelparameter(lama::DenseVector<ValueType>& vector, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filename)
{
    
    allocateModelparameter(vector,ctx,dist);
    
    readModelparameter(vector,filename);
    
    vector.redistribute(dist);
    
    writeModelparameter(vector,filename);
}


/*! \brief Write singe modelparameter to an external file
 *
 *  Write a single model to an external file.
 \param vector Single modelparameter which will be written to filename
 \param filename Name of file in which modelparameter will be written
 */
template<typename ValueType>
void Modelparameter<ValueType>::writeModelparameter(lama::DenseVector<ValueType>& vector, std::string filename)
{
    vector.writeToFile(filename);
};


/*! \brief Read a modelparameter from file
 */
template<typename ValueType>
void Modelparameter<ValueType>::readModelparameter(lama::DenseVector<ValueType>& vector, std::string filename)
{
    vector.readFromFile(filename);
};


/*! \brief Allocate a single modelparameter
 */
template<typename ValueType>
void Modelparameter<ValueType>::allocateModelparameter(lama::DenseVector<ValueType>& vector, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist)
{
    vector.setContextPtr(ctx);
    vector.allocate(dist);
};

