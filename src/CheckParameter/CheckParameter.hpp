#pragma once

#include "../Configuration/Configuration.hpp"
#include "../Modelparameter/Modelparameter.hpp"
#include "../Acquisition/Sources.hpp"

namespace KITGPI
{
    //! \brief CheckParameter namespace
    namespace CheckParameter
    {
        //! \brief Class to check different parameters.
        /*!
         * This class checks the stability criterion, numerical dispersion criterion with the parameters specified in the configuration file.
         * It also checks if the input and output files are readable and writeable. 
         */
	template <typename ValueType>
        class CheckParameter
        {
          public:
            //! Default destructor
            ~CheckParameter(){};

            //! Default constructor
	    CheckParameter() = delete; 
	    CheckParameter(Configuration::Configuration const &config, Modelparameter::Modelparameter<ValueType> &model,scai::dmemo::CommunicatorPtr comm);
										       
	    void checkStabilityCriterion(ValueType dt, ValueType dh, ValueType vpMax, std::string dimension, IndexType FDorder );
	    void checkNumericalDispersion(ValueType dh, ValueType vsMin, ValueType fc, IndexType spFDo);
	    
          private:

        };
    }
}