#include "CheckParameter.hpp"

//! \brief Constructor, which calls its memberfunctions 
template <typename ValueType>
KITGPI::CheckParameter::CheckParameter<ValueType>::CheckParameter(Configuration::Configuration const &config, Modelparameter::Modelparameter<ValueType> &model,scai::dmemo::CommunicatorPtr comm)
{
    scai::lama::Scalar vMinTmp;
    if (config.get<std::string>("equationType").compare("acoustic") == 0) {
	vMinTmp = model.getVelocityP().min();
    }  
    else {
        vMinTmp = model.getVelocityS().min();
    } 
    scai::lama::Scalar vMaxTmp = model.getVelocityP().max();
    scai::lama::DenseMatrix<ValueType> acquisition_temp;
    scai::lama::DenseVector<ValueType> wavelet_fc; 
    scai::lama::Scalar fcMax;
    acquisition_temp.readFromFile(config.get<std::string>("SourceFilename"));
    acquisition_temp.getColumn(wavelet_fc,7);
    fcMax=wavelet_fc.max();
    
    if (comm->getRank() == MASTERGPI) {
    checkStabilityCriterion(config.get<ValueType>("DT"),config.get<ValueType>("DH"), vMaxTmp.getValue<ValueType>(),config.get<std::string>("dimension"), config.get<IndexType>("spatialFDorder"));
    checkNumericalDispersion(config.get<ValueType>("DH"), vMinTmp.getValue<ValueType>(), fcMax.getValue<ValueType>(), config.get<IndexType>("spatialFDorder"));
    }
}

/*! \brief check Courant-Friedrichs-Lewy-Criterion
 *
 \param dt Temporal sampling interval in seconds. 
 \param dh Spatial sampling interval in meters.
 \param vpMax Maximum P wave velocity.
 \param dimension Dimension specified in config.txt 
 \param spFDo Spatial FD order
 */
template <typename ValueType>
void KITGPI::CheckParameter::CheckParameter<ValueType>::checkStabilityCriterion(ValueType dt, ValueType dh, ValueType vpMax, std::string dimension, IndexType spFDo)
{
    IndexType D;
    if (dimension.compare("2D") == 0) {
	D=2;
    }  
    else if (dimension.compare("3D") == 0) {
	D=3;
    }  
    else {
    SCAI_ASSERT_ERROR(false,"Unknown dimension");
    }
    
    ValueType h;
    switch(spFDo){
    case 2:
	h=1;
	break;
    case 4:
	h=7/6;
	break;
    case 6:
	h=149/120;
	break;
    case 8:
	h=2161/1680;
	break;
    case 10:
	h=53089/40320;
	break;
    case 12:
	h=1187803/887040;
	break;
    default:
	SCAI_ASSERT_ERROR(false,"Unknown spatial FD order")
    }  

    //Assess stability criterion
    SCAI_ASSERT_ERROR(dt<=dh/(h*sqrt(D)*vpMax),"\nCourant-Friedrichs-Lewy-Criterion is not met! \ndt is "<<dt<<" but should be less than dh/(h*sqrt(D)*vpMax="<< dh/(h*sqrt(D)*vpMax) <<"\n\n");
}


/*! \brief check criterion to avoid numerical dispersion
  \param dh Spatial sampling interval in meters.
  \param vMin Minimum wave velocity. (in the acoustic case the minimal P wave velocity is used; in the elastic and viscoelastic cases the minimal S wave velocity is used)
  \param fcMax Maximum center frequency of the sources in hertz.
  \param spFDo Spatial FD order.
 */

template <typename ValueType>
void KITGPI::CheckParameter::CheckParameter<ValueType>::checkNumericalDispersion(ValueType dh, ValueType vMin, ValueType fcMax, IndexType spFDo)
{
    IndexType N;
    switch(spFDo){
    case 2:
	N=12;
	break;
    case 4:
	N=8;
	break;
    case 6:
	N=7;
	break;
    case 8:
	N=6;
	break;
    case 10:
	N=5;
	break;
    case 12:
	N=4;
	break;
    default:
	SCAI_ASSERT_ERROR(false,"Unknown spatial FD order")
    }  
    if(dh>vMin/(2*fcMax*N)){
    std::cout<<"\nCriterion to avoid numerical dispersion is not met! \ndh is "<<dh<<" but should be less than vMin/(2*fcMax*N)="<<vMin/(2*fcMax*N)<<"\n\n";
    }
}


template class KITGPI::CheckParameter::CheckParameter<float>;
template class KITGPI::CheckParameter::CheckParameter<double>;