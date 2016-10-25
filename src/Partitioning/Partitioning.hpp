#pragma once

namespace KITGPI {
    
    //! Partitioning namespace
    namespace Partitioning {
        
        //! \brief Abstract class for partitioning
        /*!
         *
         */
        template <typename ValueType>
        class Partitioning
        {
        protected:
            //! Default constructor
            Partitioning(){};
            
            //! Default destructor
            ~Partitioning(){};
            
            //! Getter method for distribution pointer
            virtual dmemo::DistributionPtr getDist()=0;
            
        };
        
    }
}
