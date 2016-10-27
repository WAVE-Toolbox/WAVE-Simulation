

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

namespace KITGPI {
    
    //! \brief Modelparameter namespace
    namespace Modelparameter {
        
        //! Abstract class for a single Modelparameter (Subsurface properties)
        /*!
         * This class handels a single modelparameter.
         * As this class is an abstract class, all constructors are protected.
         */
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
            //virtual void write(std::string filenameOut)=0;
            
            //! \brief Get referenec to density model parameter
            virtual lama::DenseVector<ValueType>& getDensity()=0;
            //! \brief Get referenec to density model parameter
            virtual lama::DenseVector<ValueType>& getInverseDensity()=0;
            //! \brief Get referenec to first Lame model parameter
            virtual lama::DenseVector<ValueType>& getLambda()=0;
            //! \brief Get referenec to second Lame model parameter
            virtual lama::DenseVector<ValueType>& getMu()=0;
            //! \brief Get referenec to P-wave velocity
            virtual lama::DenseVector<ValueType>& getVelocityP()=0;
            //! \brief Get referenec to S-wave velocity
            virtual lama::DenseVector<ValueType>& getVelocityS()=0;
            
        protected:
            void initModelparameter(lama::DenseVector<ValueType>& vector, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, lama::Scalar  value);
            void initModelparameter(lama::DenseVector<ValueType>& vector, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filename);
            
            void writeModelparameter(lama::DenseVector<ValueType>& vector, std::string filename);
            
            void calculateLame(lama::DenseVector<ValueType>& vector, lama::DenseVector<ValueType>& vectorDensity, lama::DenseVector<ValueType>& vectorOut, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filename, std::string filenameDensity, std::string filenameOut);
            void calculateLame(lama::DenseVector<ValueType>& vector, lama::DenseVector<ValueType>& vectorS, lama::DenseVector<ValueType>& vectorDensity, lama::DenseVector<ValueType>& vectorOut, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filename, std::string filenameS, std::string filenameDensity, std::string filenameOut);
            
            void calculateLambda(lama::DenseVector<ValueType>& vector, lama::DenseVector<ValueType>& vectorDensity, lama::DenseVector<ValueType>& vectorOut, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filename, std::string filenameDensity, std::string filenameOut);
            void calculateLambda(lama::DenseVector<ValueType>& vector, lama::DenseVector<ValueType>& vectorS, lama::DenseVector<ValueType>& vectorDensity, lama::DenseVector<ValueType>& vectorOut, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filename, std::string filenameS, std::string filenameDensity, std::string filenameOut);
            void calculateMu(lama::DenseVector<ValueType>& vector, lama::DenseVector<ValueType>& vectorDensity, lama::DenseVector<ValueType>& vectorOut, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filename, std::string filenameDensity, std::string filenameOut);
            
        private:
            void allocateModelparameter(lama::DenseVector<ValueType>& vector, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist);
            
            void readModelparameter(lama::DenseVector<ValueType>& vector, std::string filename);
        };
    }
}

/*! \brief Init a single modelparameter by a constant value
 *
 \param vector Singel modelparameter which will be initialized
 \param ctx Context
 \param dist Distribution
 \param value Value which will be used to initialize the single modelparameter to a homogenoeus model
 */
template<typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::initModelparameter(lama::DenseVector<ValueType>& vector, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, lama::Scalar  value)
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
void KITGPI::Modelparameter::Modelparameter<ValueType>::initModelparameter(lama::DenseVector<ValueType>& vector, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filename)
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
void KITGPI::Modelparameter::Modelparameter<ValueType>::writeModelparameter(lama::DenseVector<ValueType>& vector, std::string filename)
{
    vector.writeToFile(filename);
};


/*! \brief Read a modelparameter from file
 */
template<typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::readModelparameter(lama::DenseVector<ValueType>& vector, std::string filename)
{
    vector.readFromFile(filename);
};


/*! \brief Allocate a single modelparameter
 */
template<typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::allocateModelparameter(lama::DenseVector<ValueType>& vector, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist)
{
    vector.setContextPtr(ctx);
    vector.allocate(dist);
};

/*! \brief Calculate Acoustic Lame-Vector from p-Velocity-Vector
 *  Acoustic:   lambda = rho * vP^2
 *  Elastic:    mu = rho * vS^2
 */
 template<typename ValueType>
 void KITGPI::Modelparameter::Modelparameter<ValueType>::calculateLame(lama::DenseVector<ValueType>& vector, lama::DenseVector<ValueType>& vectorDensity, lama::DenseVector<ValueType>& vectorOut, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filename, std::string filenameDensity, std::string filenameOut)
{
    allocateModelparameter(vector,ctx,dist);
    allocateModelparameter(vectorDensity,ctx,dist);
    allocateModelparameter(vectorOut,ctx,dist);
    
    readModelparameter(vector,filename);
    readModelparameter(vectorDensity,filenameDensity);
    readModelparameter(vectorOut,filenameOut);
    
    vector.redistribute(dist);
    vectorDensity.redistribute(dist);
    vectorOut.redistribute(dist);
    
    vectorOut=vector;
    vectorOut.scale(vector);
    vectorOut.scale(vectorDensity);
    
    writeModelparameter(vector,filename);
    writeModelparameter(vectorDensity,filenameDensity);
    writeModelparameter(vectorOut,filenameOut);
    
};



/*! \brief Calculate Elastic Lame-Vector (Lambda) from VelocityP and VelocityS
 *  Elastic:    lambda = rho * (vP^2 - 2 * vS^2)
 */
template<typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::calculateLame(lama::DenseVector<ValueType>& vector, lama::DenseVector<ValueType>& vectorS, lama::DenseVector<ValueType>& vectorDensity, lama::DenseVector<ValueType>& vectorOut, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filename, std::string filenameS, std::string filenameDensity, std::string filenameOut)
{
    allocateModelparameter(vector,ctx,dist);
    allocateModelparameter(vectorS,ctx,dist);
    allocateModelparameter(vectorDensity,ctx,dist);
    allocateModelparameter(vectorOut,ctx,dist);
    
    readModelparameter(vector,filename);
    readModelparameter(vectorS,filenameS);
    readModelparameter(vectorDensity,filenameDensity);
    readModelparameter(vectorOut,filenameOut);
    
    vector.redistribute(dist);
    vectorS.redistribute(dist);
    vectorDensity.redistribute(dist);
    vectorOut.redistribute(dist);
    
    vectorOut=vectorS;
    vectorOut=(-2.)*vectorOut.scale(vectorS);
    vectorOut.invert();
    vectorOut.scale(vector);
    vectorOut.invert();
    vectorOut+=vector; /* here vector.scale(vector) missing*/
    vectorOut.scale(vector);
    vectorOut=vectorOut.scale(vectorDensity);
    
    writeModelparameter(vector,filename);
    writeModelparameter(vectorS,filenameS);
    writeModelparameter(vectorDensity,filenameDensity);
    writeModelparameter(vectorOut,filenameOut);
    
};


/*! \brief Calculate Lambda (Acoustic)
 */
template<typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::calculateLambda(lama::DenseVector<ValueType>& vector, lama::DenseVector<ValueType>& vectorDensity, lama::DenseVector<ValueType>& vectorOut, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filename, std::string filenameDensity, std::string filenameOut)
{
    calculateLame(vector,vectorDensity,vectorOut,ctx,dist,filename,filenameDensity,filenameOut);
    
};


/*! \brief Calculate Lambda (Elastic)
 */
template<typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::calculateLambda(lama::DenseVector<ValueType>& vector, lama::DenseVector<ValueType>& vectorS, lama::DenseVector<ValueType>& vectorDensity, lama::DenseVector<ValueType>& vectorOut, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filename, std::string filenameS, std::string filenameDensity, std::string filenameOut)
{
    calculateLame(vector,vectorS,vectorDensity,vectorOut,ctx,dist,filename,filenameS,filenameDensity,filenameOut);
    
};


/*! \brief Calculate Mu (Elastic)
 */
template<typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::calculateMu(lama::DenseVector<ValueType>& vector, lama::DenseVector<ValueType>& vectorDensity, lama::DenseVector<ValueType>& vectorOut, hmemo::ContextPtr ctx, dmemo::DistributionPtr dist, std::string filename, std::string filenameDensity, std::string filenameOut)
{
    calculateLame(vector,vectorDensity,vectorOut,ctx,dist,filename,filenameDensity,filenameOut);
    
};







