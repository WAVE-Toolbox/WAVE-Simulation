#pragma once


namespace KITGPI {
    
    //! \brief Acquisition namespace
    namespace Acquisition {
        
        /*! \brief Struct to save 3-D coordinates
         */
        struct coordinate3D
        {
            IndexType x; //!< x Position in X-direction in grid points (Horizontal 1)
            IndexType y; //!< y Position in Y-direction in grid points (Depth)
            IndexType z; //!< z Position in Z-direction in grid points (Horizontal 2)
        };
        
        /*! \brief Struct to save 2-D coordinates
         */
        struct coordinate2D
        {
            IndexType x; //!< x Position in X-direction in grid points (Horizontal 1)
            IndexType y; //!< y Position in Z-direction in grid points (Depth)
        };
        
        /*! \brief This class manages the transformation of Coordinates
         */
        template <typename ValueType>
        class Coordinates
        {
            
        public:
            
            

            // Coordinate --> Index:
            // Interfaces 3-D
            IndexType coordinate2index(coordinate3D coordinate, IndexType NX, IndexType NY, IndexType NZ);
            IndexType coordinate2index(IndexType X, IndexType Y, IndexType Z, IndexType NX, IndexType NY, IndexType NZ);
            // Interfaces 2-D
            IndexType coordinate2index(coordinate2D coordinate, IndexType NX, IndexType NY);
            IndexType coordinate2index(IndexType X, IndexType Y, IndexType NX, IndexType NY);
            
            // Index --> Coordinate:
            coordinate3D index2coordinate(IndexType coordinate, IndexType NX, IndexType NY, IndexType NZ);
            coordinate2D index2coordinate(IndexType coordinate, IndexType NX, IndexType NY);
            
            bool locatedOnSurface(IndexType coordinate, IndexType NX, IndexType NY, IndexType NZ);
            
        protected:
            
            void Global2Local(lama::DenseVector<ValueType>& coordinatesglobal,hmemo::HArray<IndexType>& coordinateslocal, dmemo::DistributionPtr dist);
            
        private:
            
            // Coordinate --> Index:
            IndexType map3Dcoordinate2index(IndexType X, IndexType Y, IndexType Z, IndexType NX, IndexType NY, IndexType NZ);
            IndexType map2Dcoordinate2index(IndexType X, IndexType Y, IndexType NX, IndexType NY);
            
            // Index --> Coordinate:
            coordinate3D map3Dindex2coordinate(IndexType coordinate, IndexType NX, IndexType NY);
            coordinate2D map2Dindex2coordinate(IndexType coordinate, IndexType NX);
            
        };
        
    }
}
/* ------- */
/* Mapping */
/* ------- */

template <typename ValueType>
bool KITGPI::Acquisition::Coordinates<ValueType>::locatedOnSurface(IndexType coordinate, IndexType NX, IndexType NY, IndexType /*NZ*/){
    coordinate3D result;
    result=map3Dindex2coordinate(coordinate,NX,NY);
    if(result.y==1){
        return(true);
    } else {
        return(false);
    }
}

//! General mapping from 1-D coordinate to 3-D coordinate
template <typename ValueType>
KITGPI::Acquisition::coordinate3D KITGPI::Acquisition::Coordinates<ValueType>::map3Dindex2coordinate(IndexType coordinate, IndexType NX, IndexType NY)
{
    coordinate3D result;
    
    result.z=IndexType(coordinate) / (NX*NY);
    coordinate -= result.z * (NX*NY);
    result.z+=1;
    
    result.y= IndexType(coordinate) / (NX);
    coordinate -= result.y * (NX);
    result.y+=1;
    
    result.x=coordinate;
    result.x+=1;
    
    return(result);
    
}

//! General mapping from 1-D coordinate to 2-D coordinate
template <typename ValueType>
KITGPI::Acquisition::coordinate2D KITGPI::Acquisition::Coordinates<ValueType>::map2Dindex2coordinate(IndexType coordinate, IndexType NX)
{
    
    coordinate2D result;
    
    result.y= IndexType(coordinate) / (NX);
    coordinate -= result.y * (NX);
    result.y+=1;
    
    result.x=coordinate;
    result.x+=1;
    
    return(result);
    
}

//! General mapping from 1-D coordinate to 3-D coordinate
template <typename ValueType>
KITGPI::Acquisition::coordinate3D KITGPI::Acquisition::Coordinates<ValueType>::index2coordinate(IndexType coordinate, IndexType NX, IndexType NY, IndexType /*NZ*/)
{
    return(map3Dindex2coordinate(coordinate,NX,NY));
}

//! General mapping from 1-D coordinate to 2-D coordinate
template <typename ValueType>
KITGPI::Acquisition::coordinate2D KITGPI::Acquisition::Coordinates<ValueType>::index2coordinate(IndexType coordinate, IndexType NX, IndexType /*NY*/)
{
    return(map2Dindex2coordinate(coordinate,NX));
}

//! General mapping from 3-D coordinates to 1-D coordinate
template <typename ValueType>
IndexType KITGPI::Acquisition::Coordinates<ValueType>::map3Dcoordinate2index(IndexType X, IndexType Y, IndexType Z, IndexType NX, IndexType NY, IndexType NZ)
{
    if ( Z > NZ || X > NX || Y > NY || Z < 1 || X < 1 || Y < 1 )
    {
        COMMON_THROWEXCEPTION ( "Could not map from coordinate to indize!" )
        return -100;
    }
    else
    {
        return ( ( X - 1 ) + ( Y - 1 ) * NX + ( Z - 1 ) * NX * NY );
    }
}

//! General mapping from 2-D coordinates to 1-D coordinate
template <typename ValueType>
IndexType KITGPI::Acquisition::Coordinates<ValueType>::map2Dcoordinate2index(IndexType X, IndexType Y, IndexType NX, IndexType NY)
{
    if ( X > NX || Y > NY || Y < 1 || X < 1 )
    {
        COMMON_THROWEXCEPTION ( "Could not map from coordinate to indize!" )
        return -100;
    }
    else
    {
        return ( ( X - 1 ) + ( Y - 1 ) * NX );
    }
}


/* ---------- */
/* Interfaces */
/* ---------- */

/*! \brief Convert 3-D coordinates to 1-D coordinates
 \param X 3-D coordinate in X (Horizontal 1)
 \param Y 3-D coordinate in Y (Depth)
 \param Z 3-D coordinate in Z (Horizontal 2)
 \param NX Total number of grid points in X
 \param NY Total number of grid points in Y
 \param NZ Total number of grid points in Z
 \return 1-D coordinate
 */
template <typename ValueType>
IndexType KITGPI::Acquisition::Coordinates<ValueType>::coordinate2index(IndexType X, IndexType Y, IndexType Z, IndexType NX, IndexType NY, IndexType NZ)
{
    return(map3Dcoordinate2index(X,Y,Z,NX,NY,NZ));
}

/*! \brief Convert 3-D coordinates to 1-D coordinates
 \param coordinate as a coordinate3D struct
 \param NX Total number of grid points in X (Horizontal 1)
 \param NY Total number of grid points in Y (Depth)
 \param NZ Total number of grid points in Z (Horizontal 2)
 \return 1-D coordinate
 */
template <typename ValueType>
IndexType KITGPI::Acquisition::Coordinates<ValueType>::coordinate2index(coordinate3D coordinate, IndexType NX, IndexType NY, IndexType NZ)
{
    return(map3Dcoordinate2index(coordinate.x,coordinate.y,coordinate.z,NX,NY,NZ));
}

/*! \brief Convert 2-D coordinates to 1-D coordinates
 \param X 2-D coordinate in X (Horizontal 1)
 \param Y 2-D coordinate in Y (Depth)
 \param NX Total number of grid points in X (Horizontal 1)
 \param NY Total number of grid points in Y (Depth)
 \return 1-D coordinate
 */
template <typename ValueType>
IndexType KITGPI::Acquisition::Coordinates<ValueType>::coordinate2index(IndexType X, IndexType Y, IndexType NX, IndexType NY)
{
    return(map2Dcoordinate2index(X,Y,NX,NY));
}

/*! \brief Convert 2-D coordinates to 1-D coordinates
 \param coordinate as a coordinate2D struct
 \param NX Total number of grid points in X (Horizontal 1)
 \param NY Total number of grid points in Y (Depth)
 \return 1-D coordinate
 */
template <typename ValueType>
IndexType KITGPI::Acquisition::Coordinates<ValueType>::coordinate2index(coordinate2D coordinate, IndexType NX, IndexType NY)
{
    return(map2Dcoordinate2index(coordinate.x,coordinate.y,NX,NY));
}

/*! \brief Determination of local coordinates based on given global coordinates
 *
 * Calculate the number of coordinates within the local processing unit as well as
 * the coordinates of the local coordinates.
 \param coordinatesglobal DenseVector with global coordinates
 \param localIndices DenseVector with local coordinates
 \param dist Distribution of global grid
 */
template <typename ValueType>
void KITGPI::Acquisition::Coordinates<ValueType>::Global2Local(lama::DenseVector<ValueType>& coordinatesglobal,hmemo::HArray<IndexType>& localIndices, dmemo::DistributionPtr dist)
{
    
    IndexType n_global=coordinatesglobal.size(); // Number of global entries
    
    IndexType coordinatetemp_int;
    scai::lama::Scalar coordinatetemp_scalar=0;
    
    IndexType i=0;
    for(IndexType n=0; n<n_global; n++){
        
        coordinatetemp_scalar = coordinatesglobal.getValue(n);
        coordinatetemp_int=coordinatetemp_scalar.getValue<IndexType>();
        
        if( dist->isLocal(coordinatetemp_int) ) {
            i++;
        }
    }
    
    /* Determine coordinates of local receivers in the global coordinate vector */
    localIndices.resize(i);
    hmemo::WriteAccess<IndexType> write_localIndices(localIndices);
    i=0;
    for(IndexType n=0; n<n_global; n++){
        
        coordinatetemp_scalar = coordinatesglobal.getValue(n);
        coordinatetemp_int=coordinatetemp_scalar.getValue<IndexType>();
        if( dist->isLocal(coordinatetemp_int) ) {
            write_localIndices[i]=n;
            i++;
        }
    }
    
}
