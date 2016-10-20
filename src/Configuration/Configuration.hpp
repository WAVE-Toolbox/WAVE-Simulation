#pragma once

#include <scai/lama/Scalar.hpp>

#include <iostream>
#include <string>
#include <vector>
#include <sstream>

using namespace scai;

//! Namespace of the Geophysical Institute of the Karlsruhe Institute of Technology
namespace KITGPI {
    
    //! Configuration namespace
    namespace Configuration {
        
        //! Class for Configuration of the FD simulation
        /*!
         This class handels the configuration for the finite-difference simulation.
         */
        template<typename ValueType>
        class Configuration
        {
        public:
            
            /*! \brief Default deconstructor
             */
            ~Configuration(){}
            
            Configuration( std::string filename );
            
            void print();
            
            IndexType getNX() { return NX; } ///< Return NX
            IndexType getNY() { return NY; } ///< Return NY (Depth)
            IndexType getNZ() { return NZ; } ///< Return NZ

            ValueType getDH() { return DH; } ///< Return Grid spacing
            
            ValueType getDT() { return DT; } ///< Return Time Step
            ValueType getT() { return T; } ///< Return Total propagation time
            
            IndexType getModelRead() {return ModelRead;} ///< Return Read in Model?
            IndexType getModelWrite() {return ModelWrite;} ///< Return Write Model to file?
            std::string getModelFilename() {return ModelFilename;} ///< Return Filename of Model
            ValueType getVelocityP() { return velocityP; } ///< Return Velocity for homogeneous model
            ValueType getVelocityS() { return velocityS; } ///< Return Velocity for homogeneous model
            ValueType getRho() { return rho; } ///< Return Density for homogeneous model
            
            std::string getSourceFilename() {return SourceFilename;} ///< Return Filename of Source file
            std::string getReceiverFilename() {return ReceiverFilename;} ///< Return Filename of Receiver file
            std::string getSeismogramFilename() {return SeismogramFilename;} ///< Return Filename of Seismogram file
            
            IndexType getN() { return N; } ///< Return N
            
            ValueType getLambda() { return lambda; } ///< Return lambda
            ValueType getMu() { return mu; } ///< Return lambda

            IndexType getNT() { return NT; } ///< Return NT
            
            IndexType getUseCubePartitioning() { return UseCubePartitioning;} ///< Return UseCubePartitioning
            IndexType getProcNX() { return ProcNX;} ///< Return number of cores in X-direction
            IndexType getProcNY() { return ProcNY;} ///< Return number of cores in Y-direction
            IndexType getProcNZ() { return ProcNZ;} ///< Return number of cores in Z-direction

        private:
            
            IndexType NumParameters; ///< Number of parameters in input file
            
            /* read parameters */
            
            // define spatial sampling: number of grid points in direction
            IndexType NX; ///< Grid points horizontal 1
            IndexType NY;///< Grid points depth
            IndexType NZ; ///< Grid points horizontal 2
            
            /// define distance between two grid points in meter
            ValueType DH;
            
            // define temporal sampling
            ValueType DT;  ///< temporal sampling in seconds
            ValueType T;   ///< total simulation time
            
            IndexType ModelRead; ///< Read model from File (1=YES, else=NO)
            IndexType ModelWrite; ///< Write model to File (1=YES, else=NO)
            std::string ModelFilename; ///< Filename to read model
            ValueType velocityP; ///< Density in kilo gramms per cubic meter
            ValueType velocityS; ///< Density in kilo gramms per cubic meter
            ValueType rho;      ///< P-wave velocity in meter per seconds
            
            IndexType UseCubePartitioning; ///< Use cubes for partitioning of the wave fields (1=ON, else=OFF)
            IndexType ProcNX; ///< Number of cores in X-direction
            IndexType ProcNY; ///< Number of cores in Y-direction
            IndexType ProcNZ; ///< Number of cores in Z-direction
            
            std::string SourceFilename; ///< Filename to read source configuration
            std::string ReceiverFilename; ///< Filename to read receiver configuration
            std::string SeismogramFilename; ///< Filename to write seismograms
            
            
            /* calculated parameters */
            
            IndexType N; ///< Number of total grid points NX*NY*NZ
            
            ValueType lambda; ///< First Lame parameter
            ValueType mu; ///<< Second Lame parameter
            IndexType NT; ///< Number of time steps
            
        };
    }
}

/*! \brief Constructor
 *
 \param filename of configuration file
 */
template<typename ValueType>
KITGPI::Configuration::Configuration<ValueType>::Configuration( std::string filename ): NumParameters(19)
{
    // read all lines in file
    
    std::string line;
    std::map<std::string,std::string> map;
    std::ifstream input( filename.c_str() );
    
    while ( std::getline( input, line ) )
    {
        size_t lineEnd = line.size();
        std::string::size_type commentPos1 = line.find_first_of( "#", 0 );
        if ( std::string::npos != commentPos1 )
        {
            std::string::size_type commentPos2 = line.find_first_of( "#", commentPos1 );
            if ( std::string::npos != commentPos2 )
            {
                if( commentPos1 == 0 )
                {
                    continue;
                }
                lineEnd = commentPos1;
            }
        }
        
        std::string::size_type equalPos = line.find_first_of( "=", 0 );
        
        if ( std::string::npos != equalPos )
        {
            // tokenize it  name = val
            std::string name = line.substr( 0, equalPos );
            size_t len = lineEnd - ( equalPos + 1);
            std::string val  = line.substr( equalPos + 1, len);
            map.insert( std::pair<std::string,std::string>( name, val ) );
        }
    }
    input.close();
    
    // check map and assign all values with right "cast" to members
    size_t nArgs = NumParameters;
    if ( map.size() != nArgs )
    {
        std::cout << filename << " does not include a valid configutation with " << nArgs << " arguments." << std::endl;
    }
    
    std::istringstream( map[ "NZ" ] ) >> NZ; // IndexType
    std::istringstream( map[ "NX" ] ) >> NX; // IndexType
    std::istringstream( map[ "NY" ] ) >> NY; // IndexType
    
    std::istringstream( map[ "DH" ] ) >> DH; // ValueType
    
    std::istringstream( map[ "DT" ] ) >> DT; // ValueType
    std::istringstream( map[ "T" ] ) >> T;  // ValueType
    
    std::istringstream( map[ "ModelRead" ] ) >> ModelRead; // IndexType
    std::istringstream( map[ "ModelWrite" ] ) >> ModelWrite; // IndexType
    ModelFilename = std::istringstream( map[ "ModelFilename" ] ).str(); // std::string

    
    std::istringstream( map[ "velocityP" ] ) >> velocityP; // ValueType
    std::istringstream( map[ "velocityS" ] ) >> velocityS; // ValueType

    std::istringstream( map[ "rho" ] ) >> rho; // ValueType
    
    SourceFilename = std::istringstream( map[ "SourceFilename" ] ).str(); // std::string
    ReceiverFilename = std::istringstream( map[ "ReceiverFilename" ] ).str(); // std::string
    SeismogramFilename = std::istringstream( map[ "SeismogramFilename" ] ).str(); // std::string
    
    std::istringstream( map[ "UseCubePartitioning" ] ) >> UseCubePartitioning; // IndexType
    std::istringstream( map[ "ProcNX" ] ) >> ProcNX; // IndexType
    std::istringstream( map[ "ProcNY" ] ) >> ProcNY; // IndexType
    std::istringstream( map[ "ProcNZ" ] ) >> ProcNZ; // IndexType
    
    // calculate other parameters
    
    N = NZ * NX * NY;
    
    lambda =(-2*velocityS*velocityS + velocityP * velocityP ) * rho; // P-wave modulus
    mu=rho*velocityS*velocityS;
    NT = static_cast<IndexType>( ( T / DT ) + 0.5 ); // MATLAB round(T/DT)
    
    
}


/*! \brief Print configuration
 */
template<typename ValueType>
void KITGPI::Configuration::Configuration<ValueType>::print()
{
    IndexType velocity_max = velocityP; // TODO: is velocity a vector?
    double courant = velocity_max * DT / DH;
    
    std::cout << "Configuration:" << std::endl;
    std::cout << "Criteriums:" << std::endl;
    std::cout << "    Courant-number: " << courant << std::endl;
    if ( courant >= 0.8 )
    {
        std::cout << "Simulation will be UNSTABLE" << std::endl;
        std::cout << "Choose smaller DT, eg.: " << DH * 0.3 / velocity_max << std::endl;
        exit(0);
    }
    std::cout << "Modelling-domain:" << std::endl;
    std::cout << "    X: " << DH * NX << " m (Horizontal)" << std::endl;
    std::cout << "    Y: " << DH * NY << " m (Depth)" << std::endl;
    std::cout << "    Z: " << DH * NZ << " m (Horizontal)" << std::endl;

    std::cout << "Acquisition:" << std::endl;
    std::cout << "    Source acquisition will be read in from " << SourceFilename << std::endl;
    std::cout << "    Receiver acquisition will be read in from " << ReceiverFilename << std::endl;
    std::cout << "    Seismograms will be written to " << SeismogramFilename << std::endl;
    std::cout << "Material:" << std::endl;
    if(ModelRead==1) {
        std::cout << "    Model will be read in from disk" << std::endl;
        std::cout << "    First Lame-Parameter: " << ModelFilename << ".M.mtx" << std::endl;
        std::cout << "    Density: " << ModelFilename << ".density.mtx" << std::endl;
    } else {
        std::cout << "    A homogeneous model will be generated" << std::endl;
        std::cout << "    VelocityP:" << velocityP << " m/s" << std::endl;
        std::cout << "    VelocityS:" << velocityS << " m/s" << std::endl;

        std::cout << "    Density:" << rho << " g/cm3" << std::endl;
    }
    if(ModelWrite==1)  std::cout << "    The model will be written to " << ModelFilename+".out*" << std::endl;
    
    std::cout << std::endl;
}

