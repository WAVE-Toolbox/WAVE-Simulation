/*! \brief Class for Receiver handling
 *
 * This class will handle the generation and distribution of the receiver location
 */

#pragma once
#include <scai/lama.hpp>
#include <scai/lama/DenseVector.hpp>
#include <scai/lama/expression/all.hpp>
#include <scai/lama/matrix/all.hpp>

#include <iostream>

#include "Coordinates.hpp"

template<typename ValueType>
class Receiver : private Coordinates<ValueType>
{
    
public:
    Receiver(){}; //!< Default constructor
    ~Receiver(){}; //!< Default deconstructor
    
    // Initialisations of source configurations
    void DiagonalArraySurface(IndexType NX, IndexType NY, IndexType NZ, IndexType OffsetX, IndexType OffsetY, IndexType OffsetZ);
    
    void GetLocalReceiver(dmemo::DistributionPtr dist_wavefield);
    
    IndexType nreceiverglobal=0; //!< Number of receiver global
    IndexType nreceiverlocal=0; //!< Number of receiver local
    
    lama::DenseVector<ValueType> coordinatesglobal; //!< Coordinates of receivers global (1-D coordinates)
    lama::DenseVector<ValueType> coordinateslocal; //!< Coordinates of receivers local (1-D coordinates)
    
};


/*! \brief Creates a diagonal array of receiver in the X-Y plane
 *
 * This method will create an array of receiver in the X-Y plane (basically at the surface). The array will be the diagonal of the 2D-plane. 
 * The whole array will be in the depth of OffsetZ grid points. 
 * This method is temporary and will be removed in the future.
 \param NX Total number of grid points in X
 \param NY Total number of grid points in Y
 \param NZ Total number of grid points in Z (Depth)
 \param OffsetX Distance in X as grid points of start/end of array to start/end of grid
 \param OffsetY Distance in Y as grid points of start/end of array to start/end of grid
 \param OffsetZ Depth of array in grid points
 */
template<typename ValueType>
void Receiver<ValueType>::DiagonalArraySurface(IndexType NX, IndexType NY, IndexType NZ, IndexType OffsetX, IndexType OffsetY, IndexType OffsetZ )
{
    
    // Count Number of global and local receivers
    IndexType i=0;
    IndexType temp;
    for(IndexType x=OffsetX; x<=NX-OffsetX; x++){
        for(IndexType y=OffsetY; y<=NY-OffsetY; y++){
            if(x==y){
                i++;
            }
        }
    }
    
    // Allocate memory
    nreceiverglobal=i;
    coordinatesglobal.allocate(nreceiverglobal);
    
    // Save coordinates into vectors
    i=0;
    for(IndexType x=OffsetX; x<=NX-OffsetX; x++){
        for(IndexType y=OffsetY; y<=NY-OffsetY; y++){
            
            if(x==y){
                temp=this->coordinate2index(x,y,OffsetZ,NX,NY,NZ);
                coordinatesglobal.setValue( i, temp );
                i++;
            }
            
        }
        
    }
    
}

/*! \brief Determine local receiver coordinates from global distribution
 * Determines the local receiver coordinates from the global distribution
 *
 \param dist_wavefield Distribution of the wavefields
 */
template<typename ValueType>
void Receiver<ValueType>::GetLocalReceiver(dmemo::DistributionPtr dist_wavefield)
{
    this->Global2Local(coordinatesglobal,coordinateslocal,dist_wavefield);
    nreceiverlocal=coordinateslocal.size();
    
    std::cout << " Number of global receiver: " << coordinatesglobal.size() << std::endl;
    std::cout << " Number of local receiver: " << nreceiverlocal << std::endl;
}
