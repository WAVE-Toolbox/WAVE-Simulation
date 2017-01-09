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
        
	template<typename ValueType>
	  struct PMLVariables{
		  
		ValueType VMaxCPML; ///< Maximum velocity in the CPML-Layer    
		ValueType CenterFrequencyCPML; ///< Center Frequency of the Signal in the CPML-Layer 
		ValueType KMaxCPML; ///< Maximum velocity in the CPML-Layer    
		ValueType NPower; ///< Center Frequency of the Signal in the CPML-Layer 
		
		
	    };  
	    
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
            
            explicit Configuration( std::string const& filename );
            
            void print() const;
            
            bool checkConfigPlausibility();
            
            IndexType getNX() const { return NX; } ///< Return NX
            IndexType getNY() const { return NY; } ///< Return NY (Depth)
            IndexType getNZ() const { return NZ; } ///< Return NZ
            
            ValueType getDH() const { return DH; } ///< Return Grid spacing
            
            ValueType getDT() const { return DT; } ///< Return Time Step
            ValueType getT() const { return T; } ///< Return Total propagation time
            
            IndexType getSpatialFDorder() const
            {
                SCAI_ASSERT_ERROR( spatialFDorder%2 == 0 , " spatialFDorder has to be an integer of 2");
                return spatialFDorder;
            } ///< Return order of spatial FD operator
            
            IndexType getModelRead() const {return ModelRead;} ///< Return Read in Model?
            IndexType getModelWrite() const {return ModelWrite;} ///< Return Write Model to file?
            std::string getModelFilename() const {return ModelFilename;} ///< Return Filename of Model
            IndexType getModelParametrisation() const {return ModelParametrisation;} ///< Return Read Wave-Modulus-Model or Velocity-Model
            
            ValueType getVelocityP() const { return velocityP; } ///< Return Velocity for homogeneous model
            ValueType getVelocityS() const { return velocityS; } ///< Return Velocity for homogeneous model
            ValueType getRho() const { return rho; } ///< Return Density for homogeneous model
            ValueType getTauP() const { return tauP; } ///< Return tauP for homogeneous model
            ValueType getTauS() const { return tauS; } ///< Return tauS for homogeneous model
            
            std::string getSourceFilename() const {return SourceFilename;} ///< Return Filename of Source file
            std::string getReceiverFilename() const {return ReceiverFilename;} ///< Return Filename of Receiver file
            std::string getSeismogramFilename() const {return SeismogramFilename;} ///< Return Filename of Seismogram file
            IndexType getSeismogramFormat() const {return SeismogramFormat;} ///< Return Seismogram Format
            
            IndexType getN() const { return N; } ///< Return N
            
            ValueType getPWaveModulus() const { return pWaveModulus; } ///< Return pWaveModulus
            ValueType getSWaveModulus() const { return sWaveModulus; } ///< Return pWaveModulus
            
            IndexType getNT() const { return NT; } ///< Return NT
            
            IndexType getUseCubePartitioning() const { return UseCubePartitioning;} ///< Return UseCubePartitioning
            IndexType getProcNX() const { return ProcNX;} ///< Return number of cores in X-direction
            IndexType getProcNY() const { return ProcNY;} ///< Return number of cores in Y-direction
            IndexType getProcNZ() const { return ProcNZ;} ///< Return number of cores in Z-direction
            
            IndexType getFreeSurface() const { return FreeSurface;} ///< Return FreeSurface
            
            ValueType getRelaxationFrequency() const { return relaxationFrequency;} ///< Return relaxationFrequency
            IndexType getNumRelaxationMechanisms() const { return numRelaxationMechanisms;} ///< Return number of relaxation mechanism
            
            IndexType getDampingBoundary() const {return DampingBoundary;} ///< Return DampingBoundary
            IndexType getBoundaryWidth() const {return BoundaryWidth;} ///< Return BoundaryWidth
            ValueType getDampingCoeff() const {return DampingCoeff;} ///< Return DampingCoeff
            
            PMLVariables<ValueType> const& getPMLVariables() const {return PMLVar;} ///< Return PML Variables
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
            
            IndexType spatialFDorder;	///< order of spatial FD operator
            
            IndexType ModelRead; ///< Read model from File (1=YES, else=NO)
            IndexType ModelWrite; ///< Write model to File (1=YES, else=NO)
            std::string ModelFilename; ///< Filename to read model
            IndexType ModelParametrisation; ///< Read model from File (1=Wave-Modulus-Model-Vector, 2=Velocity-Vectors)
            ValueType velocityP; ///< Density in kilo gramms per cubic meter
            ValueType velocityS; ///< Density in kilo gramms per cubic meter
            ValueType rho;      ///< P-wave velocity in meter per seconds
            
            ValueType tauP; ///< Tau value for P-waves
            ValueType tauS; ///< Tau value for S-waves
            
            IndexType UseCubePartitioning; ///< Use cubes for partitioning of the wave fields (1=ON, else=OFF)
            IndexType ProcNX; ///< Number of cores in X-direction
            IndexType ProcNY; ///< Number of cores in Y-direction
            IndexType ProcNZ; ///< Number of cores in Z-direction
            
            ValueType relaxationFrequency; ///< Relaxation frequency in Hz for first relaxation mechanism
            ValueType numRelaxationMechanisms; ///< Number of relaxation mechanisms
            
            IndexType FreeSurface; ///< Use the free surface==1 or not==0
            
            IndexType DampingBoundary; ///< Use the Damping Boundary ==1 or not==0
            ValueType DampingCoeff; ///< Damping coefficient
            IndexType BoundaryWidth; ///< Width of damping boundary
            
            PMLVariables<ValueType> PMLVar;
            
            std::string SourceFilename; ///< Filename to read source configuration
            std::string ReceiverFilename; ///< Filename to read receiver configuration
            std::string SeismogramFilename; ///< Filename to write seismograms
            IndexType SeismogramFormat; ///< Seismogram Format 1=MTX 2=SU
            
            /* calculated parameters */
            
            IndexType N; ///< Number of total grid points NX*NY*NZ
            
            ValueType pWaveModulus; ///< P-wave modulus
            ValueType sWaveModulus; ///<< S-wave modulus
            IndexType NT; ///< Number of time steps
            
        };
    }
}

/*! \brief Constructor
 *
 \param filename of configuration file
 */
template<typename ValueType>
KITGPI::Configuration::Configuration<ValueType>::Configuration( std::string const& filename ):NumParameters(34),numRelaxationMechanisms(0)
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
    SCAI_ASSERT(map.size() == nArgs, filename << " does not include a valid configutation with " << nArgs << " arguments.");

    std::istringstream( map[ "NZ" ] ) >> NZ; // IndexType
    std::istringstream( map[ "NX" ] ) >> NX; // IndexType
    std::istringstream( map[ "NY" ] ) >> NY; // IndexType
    
    std::istringstream( map[ "DH" ] ) >> DH; // ValueType
    
    std::istringstream( map[ "DT" ] ) >> DT; // ValueType
    std::istringstream( map[ "T" ] ) >> T;  // ValueType
    
    std::istringstream( map[ "spatialFDorder" ] ) >> spatialFDorder;  // Indextype
    
    std::istringstream( map[ "ModelRead" ] ) >> ModelRead; // IndexType
    std::istringstream( map[ "ModelWrite" ] ) >> ModelWrite; // IndexType
    ModelFilename = std::istringstream( map[ "ModelFilename" ] ).str(); // std::string
    std::istringstream( map[ "ModelParametrisation" ] ) >> ModelParametrisation; // IndexType
    
    
    std::istringstream( map[ "velocityP" ] ) >> velocityP; // ValueType
    std::istringstream( map[ "velocityS" ] ) >> velocityS; // ValueType
    std::istringstream( map[ "rho" ] ) >> rho; // ValueType
    
    std::istringstream( map[ "tauP" ] ) >> tauP; // ValueType
    std::istringstream( map[ "tauS" ] ) >> tauS; // ValueType
    std::istringstream( map[ "relaxationFrequency" ] ) >> relaxationFrequency; // ValueType
    std::istringstream( map[ "numRelaxationMechanisms" ] ) >> numRelaxationMechanisms; // IndexType
    
    SourceFilename = std::istringstream( map[ "SourceFilename" ] ).str(); // std::string
    ReceiverFilename = std::istringstream( map[ "ReceiverFilename" ] ).str(); // std::string
    SeismogramFilename = std::istringstream( map[ "SeismogramFilename" ] ).str(); // std::string
    
    std::istringstream( map[ "SeismogramFormat" ] ) >> SeismogramFormat; // IndexType
    std::istringstream( map[ "UseCubePartitioning" ] ) >> UseCubePartitioning; // IndexType
    std::istringstream( map[ "ProcNX" ] ) >> ProcNX; // IndexType
    std::istringstream( map[ "ProcNY" ] ) >> ProcNY; // IndexType
    std::istringstream( map[ "ProcNZ" ] ) >> ProcNZ; // IndexType
    
    std::istringstream( map[ "FreeSurface" ] ) >> FreeSurface; // IndexType
    
    std::istringstream( map[ "DampingBoundary" ] ) >> DampingBoundary; // IndexType
    std::istringstream( map[ "DampingCoeff" ] ) >> DampingCoeff; // ValueType
    std::istringstream( map[ "BoundaryWidth" ] ) >> BoundaryWidth; // IndexType
    
    std::istringstream( map[ "VMaxCPML" ] ) >> PMLVar.VMaxCPML; // ValueType
    std::istringstream( map[ "CenterFrequencyCPML" ] ) >> PMLVar.CenterFrequencyCPML; // ValueType
    std::istringstream( map[ "NPower" ] ) >> PMLVar.NPower; // ValueType
    std::istringstream( map[ "KMaxCPML" ] ) >> PMLVar.KMaxCPML; // ValueType
    
    
    // calculate other parameters
    
    N = NZ * NX * NY;
    
    sWaveModulus = rho * velocityS * velocityS;
    
    pWaveModulus = rho * velocityP * velocityP; // P-wave modulus
    
    NT = static_cast<IndexType>( ( T / DT ) + 0.5 ); // MATLAB round(T/DT)
    
}


/*! \brief Print configuration
 */
template<typename ValueType>
void KITGPI::Configuration::Configuration<ValueType>::print() const
{
    IndexType velocity_max = velocityP;
    double courant = velocity_max * DT / DH;
    
    std::cout << "Configuration:" << std::endl << std::endl;
    std::cout << "Time Step: DT = " << DT << " s" << std::endl;
    std::cout << "Grid spacing: DH = " << DH << " m" << std::endl;
    std::cout << "Total simulation time: T = " << T << " s" << std::endl;
    std::cout << "Order of spatial FD operator: " << spatialFDorder << std::endl;
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
    if(FreeSurface==1){
        std::cout << "    A free surface will be set atop the model" << std::endl;
    }
    if(DampingBoundary==1){
        std::cout << "    A damping boundary will be placed at the model boundary." << std::endl;
    }
      if(DampingBoundary==2){
        std::cout << "    A CPML boundary will be placed at the model boundary." << std::endl;
    }
    if(numRelaxationMechanisms != 0 ) {
        std::cout << "    Visco-elastic simulation with " << numRelaxationMechanisms<< " Relaxation Mechanisms" << std::endl;
        std::cout << "    The relaxation frequency is " << relaxationFrequency<< " Hz" << std::endl;
    }
    std::cout << "Acquisition:" << std::endl;
    std::cout << "    Source acquisition will be read in from " << SourceFilename << std::endl;
    std::cout << "    Receiver acquisition will be read in from " << ReceiverFilename << std::endl;
    std::cout << "    Seismograms will be written to " << SeismogramFilename << std::endl;
    std::cout << "Material:" << std::endl;
    
    if ( ModelRead==1 )
    {
        switch (ModelParametrisation) {
            case 1:
                std::cout << "    Model will be read in from disk" << std::endl;
                std::cout << "    P-wave modulus: " << ModelFilename << ".pWaveModulus.mtx" << std::endl;
                std::cout << "    Density: " << ModelFilename << ".density.mtx" << std::endl;
                break;
            case 2:
                std::cout << "    Model will be read in from disk" << std::endl;
                std::cout << "    VelocityP-Data: " << ModelFilename << ".vp.mtx" << std::endl;
                std::cout << "    VelocityS-Data " << ModelFilename << ".vs.mtx" << std::endl;
                std::cout << "    Density: " << ModelFilename << ".density.mtx" << std::endl;
                break;
            default:
                COMMON_THROWEXCEPTION(" Unkown ModelParametrisation value! ")
                break;
        }
    } else {
        std::cout << "    A homogeneous model will be generated" << std::endl;
        std::cout << "    VelocityP:" << velocityP << " m/s" << std::endl;
        std::cout << "    VelocityS:" << velocityS << " m/s" << std::endl;
        std::cout << "    Density:" << rho << " g/cm3" << std::endl;
    }
    if(ModelWrite==1)  std::cout << "    The model will be written to " << ModelFilename+".out*" << std::endl;
    
    std::cout << std::endl;
}


/*! \brief Check plausibility of the configuration file
 */
template<typename ValueType>
bool KITGPI::Configuration::Configuration<ValueType>::checkConfigPlausibility()
{
    if ((spatialFDorder % 2 != 0) || (spatialFDorder < 2) || (spatialFDorder > 12))
    {
        return false;
    }
    return true;	// check is positive
}
