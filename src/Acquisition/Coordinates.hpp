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
            
            IndexType min(){
                IndexType temp=0;
                if(x<y){
                    temp=x;
                } else {
                    temp=y;
                }
                if(z<temp){
                    temp=z;
                }
                return(temp);
            }
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
            
            // Index --> Coordinate:
            coordinate3D index2coordinate(IndexType coordinate, IndexType NX, IndexType NY, IndexType NZ);

            coordinate3D edgeDistance(coordinate3D coordinate, IndexType NX, IndexType NY, IndexType NZ);
            
            bool locatedOnSurface(IndexType coordinate, IndexType NX, IndexType NY, IndexType NZ);
            
        protected:
            
            void Global2Local(lama::DenseVector<ValueType>& coordinatesglobal,hmemo::HArray<IndexType>& coordinateslocal, dmemo::DistributionPtr dist);
            
        private:
            
            // Coordinate --> Index:
            IndexType map3Dcoordinate2index(IndexType X, IndexType Y, IndexType Z, IndexType NX, IndexType NY, IndexType NZ);
            
            // Index --> Coordinate:
            coordinate3D map3Dindex2coordinate(IndexType coordinate, IndexType NX, IndexType NY);
            
            coordinate3D estimateDistanceToEdges3D(IndexType X, IndexType Y, IndexType Z, IndexType NX, IndexType NY, IndexType NZ);

            
        };
        
    }
}
/* ------- */
/* Mapping */
/* ------- */

/*! Returns bool if given coordinate is located on the surface
 *
 \param coordinate 1-D coordinate
 \param NX Number of grid points in X
 \param NY Number of grid points in Y
 \param NZ Number of grid points in Z
 *
 */
template <typename ValueType>
bool KITGPI::Acquisition::Coordinates<ValueType>::locatedOnSurface(IndexType coordinate, IndexType NX, IndexType NY, IndexType /*NZ*/){
    coordinate3D result;
    result=map3Dindex2coordinate(coordinate,NX,NY);
    if(result.y==0){
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
    
    result.y= IndexType(coordinate) / (NX);
    coordinate -= result.y * (NX);
    
    result.x=coordinate;
    
    return(result);
    
}


/*! \brief General mapping from 1-D coordinate to 3-D coordinate
 *
 * Maps a 1-D coordinate back into 3-D coordinates.
 * The 3-D grid starts at 0 and runs to (NX-1), (NY-1) or (NZ-1).
 *
 \param coordinate 1-D coordinate
 \param NX Number of grid points in X
 \param NY Number of grid points in Y
 \param NZ Number of grid points in Z
 *
 */
template <typename ValueType>
KITGPI::Acquisition::coordinate3D KITGPI::Acquisition::Coordinates<ValueType>::index2coordinate(IndexType coordinate, IndexType NX, IndexType NY, IndexType /*NZ*/)
{
    return(map3Dindex2coordinate(coordinate,NX,NY));
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

template <typename ValueType>
KITGPI::Acquisition::coordinate3D KITGPI::Acquisition::Coordinates<ValueType>::estimateDistanceToEdges3D(IndexType X, IndexType Y, IndexType Z, IndexType NX, IndexType NY, IndexType NZ)
{
    coordinate3D distance;
    
    distance.x=!((NX-X)<(X-1))?(X-1):(NX-X);
    distance.y=!((NY-Y)<(Y-1))?(Y-1):(NY-Y);
    distance.z=!((NX-Z)<(Z-1))?(Z-1):(NZ-Z);
    
    return(distance);
    
}


template <typename ValueType>
KITGPI::Acquisition::coordinate3D KITGPI::Acquisition::Coordinates<ValueType>::edgeDistance(coordinate3D coordinate, IndexType NX, IndexType NY, IndexType NZ)
{
    return(estimateDistanceToEdges3D(coordinate.x,coordinate.y,coordinate.z,NX,NY,NZ));
}



