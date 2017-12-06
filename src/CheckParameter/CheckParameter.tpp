#include "CheckParameter.hpp"

//! \brief Wrapper Function who calls checkStabilityCriterion and checkNumericalDispersion
template<typename ValueType, typename IndexType>
void KITGPI::CheckParameter::checkNumericalArtefeactsAndInstabilities(const KITGPI::Configuration::Configuration& config, Modelparameter::Modelparameter<ValueType>& model, scai::dmemo::CommunicatorPtr comm)
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
    
  KITGPI::CheckParameter::checkStabilityCriterion<ValueType,IndexType>(config.get<ValueType>("DT"),config.get<ValueType>("DH"), vMaxTmp.getValue<ValueType>(),config.get<std::string>("dimension"), config.get<IndexType>("spatialFDorder"), comm);
  KITGPI::CheckParameter::checkNumericalDispersion<ValueType,IndexType>(config.get<ValueType>("DH"), vMinTmp.getValue<ValueType>(), fcMax.getValue<ValueType>(), config.get<IndexType>("spatialFDorder"), comm);
 
}


/*! \brief check Courant-Friedrichs-Lewy-Criterion
 *
 \param dt Temporal sampling interval in seconds. 
 \param dh Spatial sampling interval in meters.
 \param vpMax Maximum P wave velocity.
 \param dimension Dimension specified in config.txt 
 \param spFDo Spatial FD order
 */
template <typename ValueType,typename IndexType>
void KITGPI::CheckParameter::checkStabilityCriterion(ValueType dt, ValueType dh, ValueType vpMax, std::string dimension, IndexType spFDo, scai::dmemo::CommunicatorPtr comm)
{
    if (comm->getRank() == MASTERGPI) {
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
}


/*! \brief check criterion to avoid numerical dispersion
  \param dh Spatial sampling interval in meters.
  \param vMin Minimum wave velocity. (in the acoustic case the minimal P wave velocity is used; in the elastic and viscoelastic cases the minimal S wave velocity is used)
  \param fcMax Maximum center frequency of the sources in hertz.
  \param spFDo Spatial FD order.
 */

template <typename ValueType,typename IndexType>
void KITGPI::CheckParameter::checkNumericalDispersion(ValueType dh, ValueType vMin, ValueType fcMax, IndexType spFDo, scai::dmemo::CommunicatorPtr comm)
{
    if (comm->getRank() == MASTERGPI) {
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
	std::cerr<<"\nCriterion to avoid numerical dispersion is not met! \ndh is "<<dh<<" but should be less than vMin/(2*fcMax*N)="<<vMin/(2*fcMax*N)<<"\n\n";
      }
    }
}

//! \brief Wrapper Function who calls checkSources and checkReceivers
template<typename ValueType, typename IndexType>
void KITGPI::CheckParameter::checkAcquisitionGeometry(Configuration::Configuration const &config, scai::dmemo::CommunicatorPtr comm)
{
    KITGPI::CheckParameter::checkReceivers<ValueType,IndexType>(config.get<IndexType>("NX"), config.get<IndexType>("NY"), config.get<IndexType>("NZ"), config.get<std::string>("ReceiverFilename"),comm);
    KITGPI::CheckParameter::checkSources<ValueType,IndexType>(config.get<IndexType>("NX"), config.get<IndexType>("NY"), config.get<IndexType>("NZ"), config.get<std::string>("SourceFilename"),comm);
}

/*! \brief check if sources are located within the grid
  \param NX Number of gridpoints in x-direction.
  \param NY Number of gridpoints in y-direction.
  \param NZ Number of gridpoints in z-direction.
  \param sourcefile Name of the source file 
 */

template <typename ValueType, typename IndexType>
void KITGPI::CheckParameter::checkSources(IndexType NX, IndexType NY, IndexType NZ, std::string sourcefile, scai::dmemo::CommunicatorPtr comm)
{
    scai::lama::DenseMatrix<ValueType> acquisition_temp;
    acquisition_temp.readFromFile(sourcefile);
    IndexType numRows = acquisition_temp.getNumRows();
    scai::lama::DenseVector<ValueType> row_temp;
    scai::lama::Scalar X_temp;
    scai::lama::Scalar Y_temp;
    scai::lama::Scalar Z_temp;
    
    for (IndexType row_ind = 0; row_ind < numRows; row_ind++) {
	acquisition_temp.getRow(row_temp,row_ind);
	X_temp=row_temp.getValue(0);
	Y_temp=row_temp.getValue(1);
	Z_temp=row_temp.getValue(2);
	if (comm->getRank() == MASTERGPI) {
	  SCAI_ASSERT_ERROR(X_temp.getValue<IndexType>()>=0,"X coordinate (defined as gridpoint) of source #"<<row_ind+1<<" must be positive!")
	  SCAI_ASSERT_ERROR(Y_temp.getValue<IndexType>()>=0,"Y coordinate (defined as gridpoint) of source #"<<row_ind+1<<" must be positive!")
	  SCAI_ASSERT_ERROR(Z_temp.getValue<IndexType>()>=0,"Z coordinate (defined as gridpoint) of source #"<<row_ind+1<<" must be positive!")
	  SCAI_ASSERT_ERROR(X_temp.getValue<IndexType>()<NX,"X coordinate (defined as gridpoint) of source #"<<row_ind+1<<" must be smaller than NX!")
	  SCAI_ASSERT_ERROR(Y_temp.getValue<IndexType>()<NY,"Y coordinate (defined as gridpoint) of source #"<<row_ind+1<<" must be smaller than NY!")
	  SCAI_ASSERT_ERROR(Z_temp.getValue<IndexType>()<NZ,"Z coordinate (defined as gridpoint) of source #"<<row_ind+1<<" must be smaller than NZ!")
	}
    }  
}


/*! \brief check if receivers are located within the grid
  \param NX Number of gridpoints in x-direction.
  \param NY Number of gridpoints in y-direction.
  \param NZ Number of gridpoints in z-direction.
  \param receiverfile Name of the receiver file 
 */

template <typename ValueType, typename IndexType>
void KITGPI::CheckParameter::checkReceivers(IndexType NX, IndexType NY, IndexType NZ, std::string receiverfile, scai::dmemo::CommunicatorPtr comm)
{
    scai::lama::DenseMatrix<ValueType> acquisition_temp;
    acquisition_temp.readFromFile(receiverfile);
    IndexType numRows = acquisition_temp.getNumRows();
    scai::lama::DenseVector<ValueType> row_temp;
    scai::lama::Scalar X_temp;
    scai::lama::Scalar Y_temp;
    scai::lama::Scalar Z_temp;
    
    for (IndexType row_ind = 0; row_ind < numRows; row_ind++) {
	acquisition_temp.getRow(row_temp,row_ind);
	X_temp=row_temp.getValue(0);
	Y_temp=row_temp.getValue(1);
	Z_temp=row_temp.getValue(2);
	if (comm->getRank() == MASTERGPI) {
	  SCAI_ASSERT_ERROR(X_temp.getValue<IndexType>()>=0,"X coordinate (defined as gridpoint) of receiver #"<<row_ind+1<<" must be positive!")
	  SCAI_ASSERT_ERROR(Y_temp.getValue<IndexType>()>=0,"Y coordinate (defined as gridpoint) of receiver #"<<row_ind+1<<" must be positive!")
	  SCAI_ASSERT_ERROR(Z_temp.getValue<IndexType>()>=0,"Z coordinate (defined as gridpoint) of receiver #"<<row_ind+1<<" must be positive!")
	  SCAI_ASSERT_ERROR(X_temp.getValue<IndexType>()<NX,"X coordinate (defined as gridpoint) of receiver #"<<row_ind+1<<" must be smaller than NX!")
	  SCAI_ASSERT_ERROR(Y_temp.getValue<IndexType>()<NY,"Y coordinate (defined as gridpoint) of receiver #"<<row_ind+1<<" must be smaller than NY!")
	  SCAI_ASSERT_ERROR(Z_temp.getValue<IndexType>()<NZ,"Z coordinate (defined as gridpoint) of receiver #"<<row_ind+1<<" must be smaller than NZ!")
	}
    }  
}


