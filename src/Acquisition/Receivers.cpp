#include "Receivers.hpp"
#include "../CheckParameter/CheckParameter.hpp"
#include "../IO/IO.hpp"

/*! \brief Init based on the configuration class and the distribution of the wavefields
 *
 \param allSettings allSettings for receivers
 \param config Configuration class, which is used to derive all required parameters
 \param modelCoordinates Coordinate class, which eg. maps 3D coordinates to 1D model indices
 \param ctx Context
 \param dist_wavefield Distribution of the wavefields
 */
template <typename ValueType>
void KITGPI::Acquisition::Receivers<ValueType>::init(std::vector<receiverSettings> allSettings, Configuration::Configuration const &config, Coordinates<ValueType> const &modelCoordinates, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist_wavefield)
{    
    /* identify seismic and EM wave */
    std::string equationType = config.get<std::string>("equationType");
    isSeismic = Common::checkEquationType<ValueType>(equationType);
    this->getSeismogramHandler().setIsSeismic(isSeismic);
    
    /*reset seismograms. This is necessary when init will be called multiple times*/
    this->getSeismogramHandler().resetSeismograms();
   
    scai::IndexType getNT = static_cast<scai::IndexType>((config.get<ValueType>("T") / config.get<ValueType>("DT")) + 0.5);

    this->setAcquisition(allSettings, modelCoordinates, dist_wavefield, ctx);

    this->initSeismogramHandler(getNT, ctx, dist_wavefield);
    this->getSeismogramHandler().setDT(config.get<ValueType>("DT"));
    this->getSeismogramHandler().setSeismoDT(config.get<ValueType>("seismoDT"));
    this->getSeismogramHandler().setInstantaneousTrace(config.getAndCatch("instantaneousTraces", 0));
}

/*! \brief Init based on the configuration class and the distribution of the wavefields
 *
 \param acquisition_matrix Acquisition matrix (contains eg. seismogram coordinates and seismogram types)
 \param config Configuration class, which is used to derive all required parameters
 \param modelCoordinates Coordinate class, which eg. maps 3D coordinates to 1D model indices
 \param ctx Context
 \param dist_wavefield Distribution of the wavefields
 */
template <typename ValueType>
void KITGPI::Acquisition::Receivers<ValueType>::init(scai::lama::DenseMatrix<ValueType> acquisition_matrix, Configuration::Configuration const &config, Coordinates<ValueType> const &modelCoordinates, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist_wavefield)
{
    std::vector<receiverSettings> allSettings;
    acqMat2settings(acquisition_matrix, allSettings, dist_wavefield);
    
    CheckParameter::checkReceivers(allSettings, modelCoordinates, dist_wavefield->getCommunicatorPtr()); 
    
    init(allSettings, config, modelCoordinates, ctx, dist_wavefield);
}

/*! \brief Init based on the configuration class and the distribution of the wavefields
 *
 \param config Configuration class, which is used to derive all required parameters
 \param modelCoordinates Coordinate class, which eg. maps 3D coordinates to 1D model indices
 \param ctx Context
 \param dist_wavefield Distribution of the wavefields
 */
template <typename ValueType>
void KITGPI::Acquisition::Receivers<ValueType>::init(Configuration::Configuration const &config, Coordinates<ValueType> const &modelCoordinates, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist_wavefield)
{
    std::vector<receiverSettings> allSettings;

    bool useStreamConfig = config.getAndCatch<bool>("useStreamConfig", false);
    SCAI_ASSERT_ERROR(!useStreamConfig, "useStreamConfig cannot be used when useReceiversPerShot = 0!");
    getAcquisitionSettings(config, allSettings);
    
    CheckParameter::checkReceivers(allSettings, modelCoordinates, dist_wavefield->getCommunicatorPtr()); 

    init(allSettings, config, modelCoordinates, ctx, dist_wavefield);
}

/*! \brief Init based on the configuration class and the distribution of the wavefields
 *
 \param config Configuration class, which is used to derive all required parameters
 \param modelCoordinates Coordinate class, which eg. maps 3D coordinates to 1D model indices
 \param ctx Context
 \param dist_wavefield Distribution of the wavefields
 */
template <typename ValueType>
void KITGPI::Acquisition::Receivers<ValueType>::init(Configuration::Configuration const &config, Coordinates<ValueType> const &modelCoordinates, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist_wavefield, scai::IndexType shotNumber, std::vector<sourceSettings<ValueType>> sourceSettingsEncode)
{
    std::vector<receiverSettings> allSettings;
    
    bool useStreamConfig = config.getAndCatch<bool>("useStreamConfig", false);
    if (useStreamConfig) {
        Configuration::Configuration configBig(config.get<std::string>("streamConfigFilename"));
        Acquisition::Coordinates<ValueType> modelCoordinatesBig(configBig);
        std::vector<Acquisition::sourceSettings<ValueType>> sourceSettingsBig;  
        Acquisition::Sources<ValueType> sources;
        ValueType shotIncr = 0; 
        sources.getAcquisitionSettings(config, shotIncr); // to get numshots
        sourceSettingsBig = sources.getSourceSettings();
        std::vector<scai::IndexType> uniqueShotNos;
        Acquisition::calcuniqueShotNo(uniqueShotNos, sourceSettingsBig);
        IndexType numshots = uniqueShotNos.size();
        shotIncr = config.getAndCatch("shotIncr", 0.0);
        sources.getAcquisitionSettings(config, shotIncr); // to get shotIndsIncr
        sourceSettingsBig = sources.getSourceSettings();
        Acquisition::calcuniqueShotNo(uniqueShotNos, sourceSettingsBig);
        
        std::vector<receiverSettings> allSettingsBig;
        if (configBig.get<scai::IndexType>("useReceiversPerShot") == 0) {
            HOST_PRINT(dist_wavefield->getCommunicatorPtr(), "useStreamConfig = 1 is not possible when useReceiversPerShot = 0!");
        } else if (configBig.get<scai::IndexType>("useReceiversPerShot") == 1) {
            getAcquisitionSettings(config, allSettingsBig, shotNumber);
        } else if (configBig.get<scai::IndexType>("useReceiversPerShot") == 2) {
            std::vector<IndexType> shotIndsIncr = sources.getShotIndsIncr();
            getAcquisitionSettings(config, allSettingsBig, shotNumber, numshots, shotIndsIncr, sourceSettingsEncode);
        }
        
        IndexType shotIndPerShot = 0;
        if (sourceSettingsEncode.size() == 0) {
            Acquisition::getuniqueShotInd(shotIndPerShot, uniqueShotNos, shotNumber);
        } else {
            Acquisition::getuniqueShotInd(shotIndPerShot, sourceSettingsEncode, shotNumber);
        }
        std::vector<Acquisition::coordinate3D> cutCoordinates;
        Acquisition::getCutCoord(config, cutCoordinates, sourceSettingsBig, modelCoordinates, modelCoordinatesBig);
        Acquisition::getSettingsPerShot<ValueType>(allSettings, allSettingsBig, cutCoordinates.at(shotIndPerShot), modelCoordinates, config.get<IndexType>("BoundaryWidth"));
    } else {
        if (config.get<scai::IndexType>("useReceiversPerShot") == 0) {
            getAcquisitionSettings(config, allSettings);
        } else if (config.get<scai::IndexType>("useReceiversPerShot") == 1) {
            getAcquisitionSettings(config, allSettings, shotNumber);
        } else if (config.get<scai::IndexType>("useReceiversPerShot") == 2) {
            std::vector<Acquisition::sourceSettings<ValueType>> sourceSettings;  
            Acquisition::Sources<ValueType> sources;
            ValueType shotIncr = 0; 
            sources.getAcquisitionSettings(config, shotIncr); // to get numshots
            sourceSettings = sources.getSourceSettings();
            std::vector<scai::IndexType> uniqueShotNos;
            Acquisition::calcuniqueShotNo(uniqueShotNos, sourceSettings);
            IndexType numshots = uniqueShotNos.size();
            shotIncr = config.getAndCatch("shotIncr", 0.0);
            sources.getAcquisitionSettings(config, shotIncr); // to get shotIndsIncr
            std::vector<IndexType> shotIndsIncr = sources.getShotIndsIncr();
            getAcquisitionSettings(config, allSettings, shotNumber, numshots, shotIndsIncr, sourceSettingsEncode);
        }
    }
    
    try {
        CheckParameter::checkReceivers(allSettings, modelCoordinates, dist_wavefield->getCommunicatorPtr());
    } catch (scai::common::Exception &e) {
        HOST_PRINT(dist_wavefield->getCommunicatorPtr(), "Message from receiver.cpp init: Error while reading receiver settings for shot " << shotNumber);
        COMMON_THROWEXCEPTION("e.what()");
    }
    
    init(allSettings, config, modelCoordinates, ctx, dist_wavefield);
}

/*! \brief Init based on the configuration class and the distribution of the wavefields
 *
 \param config Configuration class, which is used to derive all required parameters
 \param modelCoordinates Coordinate class, which eg. maps 3D coordinates to 1D model indices
 \param ctx Context
 \param dist_wavefield Distribution of the wavefields
 */
template <typename ValueType>
void KITGPI::Acquisition::Receivers<ValueType>::initWholeSpace(Configuration::Configuration const &config, Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist_wavefield, scai::lama::DenseVector<scai::IndexType> const &seismogramTypes)
{
    std::vector<receiverSettings> allSettings;
    receiverSettings thisSettings;
    std::vector<scai::IndexType> uniqueReceiverTypes;
            
    uniqueReceiverTypes.push_back(seismogramTypes[0]);
    for (int i=0; i<seismogramTypes.size(); i++) { 
        int count = 0;
        for (unsigned index=0; index < uniqueReceiverTypes.size(); index++) { 
            if (uniqueReceiverTypes[index] != seismogramTypes[i]) {
                count++;
            }
        }
        if (count != 0)
            uniqueReceiverTypes.push_back(seismogramTypes[i]);
    }  
    for (unsigned i=0; i<uniqueReceiverTypes.size(); i++) {
        for (int index=0; index < dist_wavefield->getGlobalSize(); index++) {        
            Acquisition::coordinate3D coordinate = modelCoordinates.index2coordinate(index);
            thisSettings.receiverCoords = coordinate;
            thisSettings.receiverType = uniqueReceiverTypes[i];            
            allSettings.push_back(thisSettings);
        }
    }        
    
    CheckParameter::checkReceivers(allSettings, modelCoordinates, dist_wavefield->getCommunicatorPtr());
        
    init(allSettings, config, modelCoordinates, ctx, dist_wavefield);
}

template <typename ValueType>
void KITGPI::Acquisition::Receivers<ValueType>::checkRequiredNumParameter(scai::IndexType numParameterCheck)
{
    /* Check if number of parameters is supported */
    if (numParameterCheck != 4) {
        COMMON_THROWEXCEPTION("Receivers acquisition file has an unkown format ")
    }
}

/*! \brief Gets the Acquisition Matrix
 *
 * Uses configuration to determine if sources are initialized by txt or SU and then get the Acquisition Matrix
 *
 \param config Configuration
 \param acqMat Acquisition Matrix
 */
template <typename ValueType>
void KITGPI::Acquisition::Receivers<ValueType>::getAcquisitionSettings(Configuration::Configuration const &config, std::vector<receiverSettings> &allSettings)
{
    if (config.get<bool>("initReceiverFromSU")) {
        su.readAllSettingsFromSU(allSettings, config.get<std::string>("ReceiverFilename"), config.get<ValueType>("DH"));
    } else if (config.get<scai::IndexType>("useReceiversPerShot") == 0 || config.get<scai::IndexType>("useReceiversPerShot") == 2) {
        readAllSettings(allSettings, config.get<std::string>("ReceiverFilename") + ".txt");
    }
}

/*! \brief Gets the Acquisition Matrix
 *
 * Uses configuration to determine if sources are initialized by txt or SU and then get the Acquisition Matrix
 *
 \param config Configuration
 \param acqMat Acquisition Matrix
 */
template <typename ValueType>
void KITGPI::Acquisition::Receivers<ValueType>::getAcquisitionSettings(Configuration::Configuration const &config, std::vector<receiverSettings> &allSettings, scai::IndexType shotNumber)
{
    bool useStreamConfig = config.getAndCatch<bool>("useStreamConfig", false);
    std::string filenameTmp;
    if (!useStreamConfig) {
        filenameTmp = config.get<std::string>("ReceiverFilename") + ".shot_" + std::to_string(shotNumber);
    } else {
        Configuration::Configuration configBig(config.get<std::string>("streamConfigFilename"));
        filenameTmp = configBig.get<std::string>("ReceiverFilename") + ".shot_" + std::to_string(shotNumber);
    }
    if (config.get<bool>("initReceiverFromSU")) {
        su.readAllSettingsFromSU(allSettings, filenameTmp, config.get<ValueType>("DH"));
    } else {
        readAllSettings(allSettings, filenameTmp + ".txt");
    }
}

/*! \brief Gets the Acquisition Matrix
 *
 * Uses configuration to determine if sources are initialized by txt or SU and then get the Acquisition Matrix
 *
 \param config Configuration
 \param allSettings receiverSettings
 \param numshots numshots is the number of shots not encoded
 \param shotIndsIncr selected shot indices, shotIndsIncr.size() <= numshots
 \param sourceSettingsEncode sourceSettings for encoded receivers
 */
template <typename ValueType>
void KITGPI::Acquisition::Receivers<ValueType>::getAcquisitionSettings(Configuration::Configuration const &config, std::vector<receiverSettings> &allSettings, scai::IndexType shotNumber, scai::IndexType numshots, std::vector<IndexType> shotIndsIncr, std::vector<sourceSettings<ValueType>> sourceSettingsEncode)
{
    std::vector<receiverSettings> receiverSettings;
    getAcquisitionSettings(config, receiverSettings);
    bool useStreamConfig = config.getAndCatch<bool>("useStreamConfig", false);
    std::string filenameTmp;
    if (!useStreamConfig) {
        filenameTmp = config.get<std::string>("ReceiverFilename") + ".mark";
    } else {
        Configuration::Configuration configBig(config.get<std::string>("streamConfigFilename"));
        filenameTmp = configBig.get<std::string>("ReceiverFilename") + ".mark";
    }
    scai::lama::CSRSparseMatrix<ValueType> receiverMarkMatrix;
    scai::IndexType numrecs = receiverSettings.size();
    receiverMarkMatrix.allocate(numshots, numrecs+1); 
    IO::readMatrix<ValueType>(receiverMarkMatrix, filenameTmp, 1); // .mtx file
        
    scai::IndexType numshotsIncr = shotIndsIncr.size();
    // config is configBig when useStreamConfig != 0, so we introduce sourceSettingsEncode.size() == 0 to ensure that useSourceEncode != 0 defined in configBig does not affect useStreamConfig != 0 defined in stream config.
    if (sourceSettingsEncode.size() == 0) {// normal shots
        for (scai::IndexType shotInd = 0; shotInd < numshotsIncr; shotInd++) {
            receiverMarkMatrix.getRow(receiverMarkVector, shotIndsIncr[shotInd]);
            if (receiverMarkVector[0] == shotNumber) {
                break;
            } else if (shotIndsIncr[shotInd] == numshots - 1) {
                SCAI_ASSERT(receiverMarkVector[0] == shotNumber, "receiverMarkVector[0] != shotNumber");
            }
        }  
    } else {// in case of supershot
        scai::lama::SparseVector<ValueType> receiverMarkRow(numrecs+1, 0);
        receiverMarkVector = receiverMarkRow;
        for (scai::IndexType shotInd = 0; shotInd < numshotsIncr; shotInd++) {
            if (std::abs(sourceSettingsEncode[shotInd].sourceNo) == shotNumber) {
                receiverMarkMatrix.getRow(receiverMarkRow, shotIndsIncr[shotInd]); 
                receiverMarkVector += receiverMarkRow;
            }
        }       
    }  
    receiverMarkVector.unaryOp(receiverMarkVector, common::UnaryOp::SIGN);
    receiverMarkVector[0] = shotNumber;
    allSettings.clear();
    for (scai::IndexType irec = 0; irec < numrecs; irec++) {
        if (receiverMarkVector[irec+1] != 0) {
            allSettings.push_back(receiverSettings[irec]);
        }
    }     
}

/*! \brief write receiverMarkVector
 *
 \param config Configuration
 \param allSettings receiverSettings
 \param dist_wavefield dist_wavefield
 */
template <typename ValueType>
void KITGPI::Acquisition::Receivers<ValueType>::writeReceiverMark(Configuration::Configuration const &config, IndexType shotNumber, IndexType stage, IndexType iteration)
{
    IndexType useSourceEncode = config.getAndCatch("useSourceEncode", 0);
    if (useSourceEncode != 0) {
        bool useStreamConfig = config.getAndCatch<bool>("useStreamConfig", false);
        std::string filenameTmp;
        if (!useStreamConfig) {
            filenameTmp = config.get<std::string>("ReceiverFilename");
        } else {
            Configuration::Configuration configBig(config.get<std::string>("streamConfigFilename"));
            filenameTmp = configBig.get<std::string>("ReceiverFilename");
        }
        if (stage != 0)
            filenameTmp += ".stage_" + std::to_string(stage) + ".It_" + std::to_string(iteration);            
        
        scai::lama::DenseVector<ValueType> temp;
        temp = receiverMarkVector;
        IO::writeVector(temp, filenameTmp + ".shot_" + std::to_string(shotNumber) + ".mark", 1);
    }
}

/*! \brief Gets the Acquisition Matrix
 *
 * Uses configuration to determine if sources are initialized by txt or SU and then get the Acquisition Matrix
 *
 \param config Configuration
 \param allSettings receiverSettings
 \param dist_wavefield dist_wavefield
 */
template <typename ValueType>
void KITGPI::Acquisition::Receivers<ValueType>::acqMat2settings(scai::lama::DenseMatrix<ValueType> &acqMat, std::vector<receiverSettings> &allSettings, scai::dmemo::DistributionPtr dist_wavefield)
{
    scai::IndexType nRows = acqMat.getNumRows();
    allSettings.reserve(nRows);
    if (dist_wavefield->getCommunicator().getRank() == 0) {
        auto read_acquisition_temp_HA = hostReadAccess(acqMat.getLocalStorage().getData());
        for (scai::IndexType row = 0; row < nRows; row++) {
            allSettings[row].receiverCoords.x = read_acquisition_temp_HA[row * 9 + 0];
            allSettings[row].receiverCoords.y = read_acquisition_temp_HA[row * 9 + 1];
            allSettings[row].receiverCoords.z = read_acquisition_temp_HA[row * 9 + 2];
            allSettings[row].receiverType = read_acquisition_temp_HA[row * 9 + 3];
        }
    }
}

/*! \brief Decode the supershot to individual shots (encodeType = 2 or 3) or encode individual shots to the supershot (encodeType = 0 or 1)
 *
 \param config Configuration
 \param numshots numshots is the number of shots not encoded
 */
template <typename ValueType>
void KITGPI::Acquisition::Receivers<ValueType>::encode(Configuration::Configuration const &config, std::string const &filename, scai::IndexType shotNumber, std::vector<sourceSettings<ValueType>> sourceSettingsEncode, scai::IndexType encodeType)
{
    IndexType useSourceEncode = config.getAndCatch("useSourceEncode", 0);
    if (useSourceEncode != 0) {
        SCAI_ASSERT_ERROR(config.get<scai::IndexType>("useReceiversPerShot") == 2, "useReceiversPerShot != 2"); 
        std::vector<Acquisition::sourceSettings<ValueType>> sourceSettings;  
        Acquisition::Sources<ValueType> sources;
        ValueType shotIncr = 0; 
        sources.getAcquisitionSettings(config, shotIncr); // to get numshots
        sourceSettings = sources.getSourceSettings();
        std::vector<scai::IndexType> uniqueShotNos;
        Acquisition::calcuniqueShotNo(uniqueShotNos, sourceSettings);
        IndexType numshots = uniqueShotNos.size();
        shotIncr = config.getAndCatch("shotIncr", 0.0);
        sources.getAcquisitionSettings(config, shotIncr); // to get shotIndsIncr
        std::vector<IndexType> shotIndsIncr = sources.getShotIndsIncr();
        
        std::vector<receiverSettings> receiverSettings;
        getAcquisitionSettings(config, receiverSettings);
        bool useStreamConfig = config.getAndCatch<bool>("useStreamConfig", false);
        std::string filenameTmp;
        if (!useStreamConfig) {
            filenameTmp = config.get<std::string>("ReceiverFilename") + ".mark";
        } else {
            Configuration::Configuration configBig(config.get<std::string>("streamConfigFilename"));
            filenameTmp = configBig.get<std::string>("ReceiverFilename") + ".mark";
        }
        scai::lama::CSRSparseMatrix<ValueType> receiverMarkMatrix;
        scai::IndexType numrecs = receiverSettings.size();
        receiverMarkMatrix.allocate(numshots, numrecs+1);        
        IO::readMatrix<ValueType>(receiverMarkMatrix, filenameTmp, 1); // .mtx file
        
        std::vector<IndexType> receiverTypes;
        IndexType count = 0;
        IndexType countEncode = 0;
        for (scai::IndexType irec = 0; irec < numrecs; irec++) {
            if (receiverMarkVector[irec+1] != 0) { // to count receiverTypes (0,1,2,3)
                if (count == 0) {
                    receiverTypes.push_back(receiverSettings[irec].getType());
                } else {                    
                    if (std::find(receiverTypes.begin(), receiverTypes.end(), receiverSettings[irec].receiverType) != receiverTypes.end()) {
                        // receiverType already included
                    } else { 
                        receiverTypes.push_back(receiverSettings[irec].getType());
                    }
                }
                count++;
            }
        }  
        IndexType numrecTypes = receiverTypes.size();  
        scai::lama::DenseVector<ValueType> tempRow;
        IndexType seismoFormat = config.get<scai::IndexType>("SeismogramFormat");
        
        IndexType numshotsIncr = shotIndsIncr.size();
        scai::lama::SparseVector<ValueType> receiverMarkRow(numrecs+1, 0);        
        IndexType gradientDomain = config.getAndCatch("gradientDomain", 0);
        Filter::Filter<ValueType> freqFilter;
        if (gradientDomain != 0) {
            IndexType tStepEnd = static_cast<IndexType>((config.get<ValueType>("T") / config.get<ValueType>("DT")) + 0.5);
            freqFilter.init(config.get<ValueType>("DT"), tStepEnd);
        }
        
        for (scai::IndexType i = 0; i < numrecTypes; i++) {
            scai::lama::DenseMatrix<ValueType> dataSingle;
            scai::lama::DenseMatrix<ValueType> &data = this->getSeismogramHandler().getSeismogram(static_cast<SeismogramType>(receiverTypes[i]-1)).getData();
            std::vector<scai::lama::DenseMatrix<ValueType>> &dataDecode = this->getSeismogramHandler().getSeismogram(static_cast<SeismogramType>(receiverTypes[i]-1)).getDataDecode();
            IndexType countDecode = 0;
            if (encodeType == 0 || encodeType == 1) { // encode 2->1
                data.scale(0);
            } else {
                dataDecode.clear();
                for (scai::IndexType shotInd = 0; shotInd < numshotsIncr; shotInd++) {
                    if (std::abs(sourceSettingsEncode[shotInd].sourceNo) == shotNumber) {
                        scai::lama::DenseMatrix<ValueType> temp;
                        dataDecode.push_back(temp);
                    }
                }
            }
            if (numshots == numrecs) { // common offset data
                dataSingle.allocate(1, data.getNumColumns());
                countEncode = 0;
                if (encodeType == 0 || encodeType == 1) { // encode 1->2->1, similar with encode in multi offset data
                    if (encodeType == 1) {
                        if (this->getSeismogramHandler().getIsSeismic())
                            filenameTmp = filename + "." + SeismogramTypeString[static_cast<SeismogramType>(receiverTypes[i]-1)];
                        else
                            filenameTmp = filename + "." + SeismogramTypeStringEM[static_cast<SeismogramTypeEM>(receiverTypes[i]-1)];
                    }
                    for (scai::IndexType shotInd = 0; shotInd < numshotsIncr; shotInd++) {
                        if (std::abs(sourceSettingsEncode[shotInd].sourceNo) == shotNumber && receiverMarkVector[shotIndsIncr[shotInd]+1] != 0 && receiverSettings[shotIndsIncr[shotInd]].getType() == receiverTypes[i]) {
                            if (encodeType == 1) {
                                hmemo::HArray<ValueType> localsignal = IO::readMatrix<ValueType>(filenameTmp, shotIndsIncr[shotInd], seismoFormat);
                                dataSingle.setLocalRow(localsignal, 0, scai::common::BinaryOp::COPY);
                            } else {
                                dataSingle = dataDecode[countDecode];
                            }
                            receiverMarkMatrix.getRow(receiverMarkRow, shotIndsIncr[shotInd]); 
                            if (receiverMarkRow[shotIndsIncr[shotInd]+1] != 0) {
                                dataSingle.getRow(tempRow, 0);                      
                                if (sourceSettingsEncode[shotInd].amp < 0)
                                    tempRow *= -1;
                                if (gradientDomain != 0) {
                                    freqFilter.calc("ideal", "bp", 1, sourceSettingsEncode[shotInd].fc);
                                    freqFilter.apply(tempRow);
                                }
                                data.setRow(tempRow, countEncode, scai::common::BinaryOp::ADD);
                            }
                            countEncode++;
                            countDecode++;
                        }
                    }
                } else if (encodeType == 2 || encodeType == 3) { // decode 1->2->1, similar with decode in multi offset data
                    for (scai::IndexType shotInd = 0; shotInd < numshotsIncr; shotInd++) {
                        if (receiverMarkVector[shotIndsIncr[shotInd]+1] != 0 && receiverSettings[shotIndsIncr[shotInd]].getType() == receiverTypes[i]) {
                            receiverMarkMatrix.getRow(receiverMarkRow, shotIndsIncr[shotInd]); 
                            if (receiverMarkRow[shotIndsIncr[shotInd]+1] != 0) {
                                data.getRow(tempRow, countEncode); 
                                if (sourceSettingsEncode[shotInd].amp < 0)
                                    tempRow *= -1;
                                if (gradientDomain != 0) {
                                    freqFilter.calc("ideal", "bp", 1, sourceSettingsEncode[shotInd].fc);
                                    freqFilter.apply(tempRow);
                                }
                                dataSingle.setRow(tempRow, 0, scai::common::BinaryOp::COPY);
                            }
                            dataDecode[countDecode] = dataSingle;
                            countEncode++;
                            countDecode++;
                        }
                    }  
                }
            } else { // multi offset data
                for (scai::IndexType shotInd = 0; shotInd < numshotsIncr; shotInd++) {
                    if (std::abs(sourceSettingsEncode[shotInd].sourceNo) == shotNumber) {
                        receiverMarkMatrix.getRow(receiverMarkRow, shotIndsIncr[shotInd]); 
                        IndexType sourceNo = receiverMarkRow[0];
                        receiverMarkRow[0] = 0;
                        IndexType numrecSingle = receiverMarkRow.sum();
                        count = 0;
                        countEncode = 0;
                        dataSingle.allocate(numrecSingle, data.getNumColumns());
                        if (encodeType == 0 || encodeType == 1) { // encode 2->1
                            if (encodeType == 1) {
                                if (this->getSeismogramHandler().getIsSeismic())
                                    filenameTmp = filename + ".shot_" + std::to_string(sourceNo) + "." + SeismogramTypeString[static_cast<SeismogramType>(receiverTypes[i]-1)];
                                else
                                    filenameTmp = filename + ".shot_" + std::to_string(sourceNo)  + "." + SeismogramTypeStringEM[static_cast<SeismogramTypeEM>(receiverTypes[i]-1)];
                                
                                IO::readMatrix(dataSingle, filenameTmp, seismoFormat);  
                            } else {
                                dataSingle = dataDecode[countDecode];
                            }
                            if (sourceSettingsEncode[shotInd].amp < 0)
                                dataSingle *= -1;
                            if (gradientDomain != 0) {
                                freqFilter.calc("ideal", "bp", 1, sourceSettingsEncode[shotInd].fc);
                                freqFilter.apply(dataSingle);
                            }
                            for (scai::IndexType irec = 0; irec < numrecs; irec++) {
                                if (receiverMarkVector[irec+1] != 0 && receiverSettings[irec].getType() == receiverTypes[i]) {
                                    if (receiverMarkRow[irec+1] != 0) {
                                        dataSingle.getRow(tempRow, count); 
                                        data.setRow(tempRow, countEncode, scai::common::BinaryOp::ADD);
                                        count++;
                                    }
                                    countEncode++;
                                }
                            }  
                        } else if (encodeType == 2 || encodeType == 3) { // decode 1->2
                            for (scai::IndexType irec = 0; irec < numrecs; irec++) {
                                if (receiverMarkVector[irec+1] != 0 && receiverSettings[irec].getType() == receiverTypes[i]) {
                                    if (receiverMarkRow[irec+1] != 0) {
                                        data.getRow(tempRow, countEncode); 
                                        dataSingle.setRow(tempRow, count, scai::common::BinaryOp::COPY);
                                        count++;
                                    }
                                    countEncode++;
                                }
                            }  
                            if (sourceSettingsEncode[shotInd].amp < 0)
                                dataSingle *= -1;
                            if (gradientDomain != 0) {
                                freqFilter.calc("ideal", "bp", 1, sourceSettingsEncode[shotInd].fc);
                                freqFilter.apply(dataSingle);
                            }
                            dataDecode[countDecode] = dataSingle;
                            if (encodeType == 3) {
                                if (this->getSeismogramHandler().getIsSeismic())
                                    filenameTmp = filename + ".shot_" + std::to_string(sourceNo) + "." + SeismogramTypeString[static_cast<SeismogramType>(receiverTypes[i]-1)];
                                else
                                    filenameTmp = filename + ".shot_" + std::to_string(sourceNo)  + "." + SeismogramTypeStringEM[static_cast<SeismogramTypeEM>(receiverTypes[i]-1)];
                                
                                IO::writeMatrix(dataSingle, filenameTmp, seismoFormat);
                            }
                        }
                        countDecode++;
                    }
                }
            }  
        }
    }
}

/*! \brief Decode the supershot to individual shots (encodeType = 0 or 1)
 *
 \param config Configuration
 \param numshots numshots is the number of shots not encoded
 */
template <typename ValueType>
void KITGPI::Acquisition::Receivers<ValueType>::decode(Configuration::Configuration const &config, std::string const &filename, scai::IndexType shotNumber, std::vector<sourceSettings<ValueType>> sourceSettingsEncode, scai::IndexType encodeType)
{
    encodeType += 2; // 2 or 3
    encode(config, filename, shotNumber, sourceSettingsEncode, encodeType);
}

/*! \brief Gets the size of modelPerShot
 *
 * Uses configuration to determine the size of modelPerShot
 *
 \param config Configuration
 */
template <typename ValueType>
void KITGPI::Acquisition::Receivers<ValueType>::getModelPerShotSize(scai::dmemo::CommunicatorPtr commAll, Configuration::Configuration const &config, ValueType &NXPerShot, scai::IndexType &numShotPerSuperShot)
{
    NXPerShot = config.get<ValueType>("NX");
    numShotPerSuperShot = 1;
    if (config.get<scai::IndexType>("useReceiversPerShot") > 0) {
        double start_t, end_t;             /* For timing */
        start_t = common::Walltime::get();
        HOST_PRINT(commAll, "Start calculating the size of modelPerShot\n");
        
        bool useStreamConfig = config.getAndCatch("useStreamConfig", false);
        IndexType useSourceEncode = config.getAndCatch("useSourceEncode", 0);
        Acquisition::Sources<ValueType> sources;
        std::vector<sourceSettings<ValueType>> sourceSettingsBig; 
        Acquisition::Coordinates<ValueType> modelCoordinates(config);     
        std::vector<receiverSettings> receiverSettings; 
        
        ValueType shotIncr = 0;        
        sources.getAcquisitionSettings(config, shotIncr); // to get numshots
        sourceSettingsBig = sources.getSourceSettings(); 
        std::vector<scai::IndexType> uniqueShotNos;
        Acquisition::calcuniqueShotNo(uniqueShotNos, sourceSettingsBig);
        IndexType numshots = uniqueShotNos.size();
        shotIncr = config.getAndCatch("shotIncr", 0.0);
        sources.getAcquisitionSettings(config, shotIncr); // to get shotIndsIncr
        sourceSettingsBig = sources.getSourceSettings();
        std::vector<IndexType> shotIndsIncr = sources.getShotIndsIncr();
        
        IndexType seedtime = 0;
        sources.calcSourceSettingsEncode(commAll, config, seedtime);
        std::vector<sourceSettings<ValueType>> sourceSettingsEncode; 
        std::vector<Acquisition::coordinate3D> cutCoordinates;
        IndexType shotNumber;
        IndexType minXsrc = sourceSettingsBig[0].sourceCoords.x;
        IndexType maxXsrc = sourceSettingsBig[0].sourceCoords.x;
        IndexType minXrec = sourceSettingsBig[0].sourceCoords.x;
        IndexType maxXrec = sourceSettingsBig[0].sourceCoords.x;
        if (useStreamConfig) {
            Configuration::Configuration configBig(config.get<std::string>("streamConfigFilename"));
            Acquisition::Coordinates<ValueType> modelCoordinatesBig(configBig); 
            Acquisition::getCutCoord(config, cutCoordinates, sourceSettingsBig, modelCoordinates, modelCoordinatesBig);
            minXsrc -= cutCoordinates.at(0).x;
            maxXsrc -= cutCoordinates.at(0).x;
            minXrec -= cutCoordinates.at(0).x;
            maxXrec -= cutCoordinates.at(0).x;
        }
        IndexType shotSkip = 1;
        if (sourceSettingsBig.size()/10 > 1)   // we select 10 shots to determinate the modelPerShot size.
            shotSkip = sourceSettingsBig.size()/10;
        for (unsigned shotInd = 0; shotInd < sourceSettingsBig.size(); shotInd += shotSkip) {  
            if (useSourceEncode == 0) {
                shotNumber = std::abs(sourceSettingsBig[shotInd].sourceNo);
                if (useStreamConfig)
                    sourceSettingsBig[shotInd].sourceCoords.x -= cutCoordinates.at(shotInd).x;
                if (minXsrc > sourceSettingsBig[shotInd].sourceCoords.x)
                    minXsrc = sourceSettingsBig[shotInd].sourceCoords.x;
                if (maxXsrc < sourceSettingsBig[shotInd].sourceCoords.x)
                    maxXsrc = sourceSettingsBig[shotInd].sourceCoords.x;
            } else {
                sourceSettingsEncode = sources.getSourceSettingsEncode();
                shotNumber = std::abs(sourceSettingsEncode[shotInd].sourceNo);
                if (useStreamConfig)
                    sourceSettingsEncode[shotInd].sourceCoords.x -= cutCoordinates.at(shotInd).x;
                if (minXsrc > sourceSettingsEncode[shotInd].sourceCoords.x)
                    minXsrc = sourceSettingsEncode[shotInd].sourceCoords.x;
                if (maxXsrc < sourceSettingsEncode[shotInd].sourceCoords.x)
                    maxXsrc = sourceSettingsEncode[shotInd].sourceCoords.x;
            }
            if (config.get<scai::IndexType>("useReceiversPerShot") == 1) {
                getAcquisitionSettings(config, receiverSettings, shotNumber);
            } else if (config.get<scai::IndexType>("useReceiversPerShot") == 2) {
                getAcquisitionSettings(config, receiverSettings, shotNumber, numshots, shotIndsIncr, sourceSettingsEncode);
            }
            
            for (unsigned i = 0; i < receiverSettings.size(); i++) {  
                if (useStreamConfig)
                    receiverSettings[i].receiverCoords.x -= cutCoordinates.at(shotInd).x;
                if (minXrec > receiverSettings[i].receiverCoords.x)
                    minXrec = receiverSettings[i].receiverCoords.x;
                if (maxXrec < receiverSettings[i].receiverCoords.x)
                    maxXrec = receiverSettings[i].receiverCoords.x;
            }
        }
            
        if (maxXrec < maxXsrc)
            maxXrec = maxXsrc;
        if (minXrec > minXsrc)
            minXrec = minXsrc;
        
        IndexType BoundaryWidth = config.get<IndexType>("BoundaryWidth");    
        SCAI_ASSERT_ERROR(minXrec >= BoundaryWidth, "The minimum source or receiver x index is " << minXrec);
        if (NXPerShot < maxXrec + BoundaryWidth + 1) {
            NXPerShot = maxXrec + BoundaryWidth + 1;
            HOST_PRINT(commAll, "NXPerShot = " << NXPerShot << "\n");
        }

        IndexType numShotDomains = config.get<IndexType>("NumShotDomains");
        Common::checkNumShotDomains(numShotDomains, commAll);
        if (useSourceEncode != 0) {
            IndexType numshotsIncr = shotIndsIncr.size();
            numShotPerSuperShot = ceil(ValueType(numshotsIncr) / numShotDomains);
            HOST_PRINT(commAll, "numShotPerSuperShot = " << numShotPerSuperShot << "\n");
        }
        
        // to get numTracesGlobal for checking COP data
        Acquisition::Coordinates<ValueType> modelCoordinates2; 
        modelCoordinates2.init(config, 1, NXPerShot);
        IndexType shotDomain = numShotDomains; 
        dmemo::CommunicatorPtr commShot = commAll->split(shotDomain);
        common::Grid3D grid(config.get<IndexType>("NY"), config.get<IndexType>("NZ"), NXPerShot);
        dmemo::DistributionPtr dist = std::make_shared<dmemo::GridDistribution>(grid, commShot);
        
        CheckParameter::checkReceivers(receiverSettings, modelCoordinates2, dist->getCommunicatorPtr()); 

        hmemo::ContextPtr ctx = hmemo::Context::getContextPtr();
        init(receiverSettings, config, modelCoordinates2, ctx, dist);
        
        end_t = common::Walltime::get();
        HOST_PRINT(commAll, "Finish calculating the size of modelPerShot in: " << end_t - start_t << " sec.\n");
    }
}

template class KITGPI::Acquisition::Receivers<double>;
template class KITGPI::Acquisition::Receivers<float>;
