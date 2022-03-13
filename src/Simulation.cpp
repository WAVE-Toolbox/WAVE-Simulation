#include <scai/common/Settings.hpp>
#include <scai/common/Walltime.hpp>
#include <scai/dmemo/CommunicatorStack.hpp>
#include <scai/dmemo/GenBlockDistribution.hpp>
#include <scai/dmemo/GridDistribution.hpp>
#include <scai/lama.hpp>
#include <scai/tracing.hpp>
#include <scai/common/macros/assert.hpp>

#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>

#include "Configuration/Configuration.hpp"
#include "Configuration/ValueType.hpp"

#include "Acquisition/Receivers.hpp"
#include "Acquisition/Sources.hpp"

#include "Modelparameter/ModelparameterFactory.hpp"
#include "Wavefields/WavefieldsFactory.hpp"
#include "ForwardSolver/Derivatives/DerivativesFactory.hpp"
#include "ForwardSolver/ForwardSolverFactory.hpp"

#include "CheckParameter/CheckParameter.hpp"
#include "Common/HostPrint.hpp"
#include "Common/Common.hpp"
#include "Partitioning/Partitioning.hpp"
#include <scai/lama/io/PartitionIO.hpp>

using namespace scai;
using namespace KITGPI;

extern bool verbose; // global variable definition

int main(int argc, const char *argv[])
{
    // parse command line arguments to be set as environment variables, e.g.
    // --SCAI_CONTEXT=CUDA --SCAI_SETTINGS=domains.txt

    common::Settings::parseArgs(argc, argv);

    double start_t, end_t;             /* For timing */
    double globalStart_t, globalEnd_t; /* For timing */
    globalStart_t = common::Walltime::get();
    IndexType seedtime = (int)time(0);

    if (argc != 2) {
        std::cout << "\n\nNo configuration file given!\n\n"
                  << std::endl;
        return (2);
    }

    /* --------------------------------------- */
    /* Read configuration from file            */
    /* --------------------------------------- */
    Configuration::Configuration config(argv[1]);
    verbose = config.get<bool>("verbose");
    
    Configuration::Configuration configBig;
    bool useStreamConfig = config.getAndCatch("useStreamConfig", false);
    
    if (useStreamConfig){
        configBig.readFromFile(config.get<std::string>("streamConfigFilename"), true);    
    }
    
    std::string dimension = config.get<std::string>("dimension");
    std::string equationType = config.get<std::string>("equationType");
    IndexType numRelaxationMechanisms = config.get<IndexType>("numRelaxationMechanisms");

    /* inter node communicator */
    dmemo::CommunicatorPtr commAll = dmemo::Communicator::getCommunicatorPtr(); // default communicator, set by environment variable SCAI_COMMUNICATOR
    common::Settings::setRank(commAll->getNodeRank());

    HOST_PRINT(commAll, "\n WAVE-Simulation " << dimension << " " << equationType << " - LAMA Version\n");
    HOST_PRINT(commAll, "", "  - Running on " << commAll->getSize() << " mpi processes -\n\n");

    if (commAll->getRank() == MASTERGPI) {
        config.print();
    }

    SCAI_REGION("WAVE-Simulation.main")

    std::string settingsFilename; // filename for processor specific settings
    if (common::Settings::getEnvironment(settingsFilename, "SCAI_SETTINGS")) {
        // each processor reads line of settings file that matches its node name and node rank
        common::Settings::readSettingsFile(settingsFilename.c_str(), commAll->getNodeName(), commAll->getNodeRank());
    }
        
    /* --------------------------------------- */
    /* coordinate mapping (3D<->1D)            */
    /* --------------------------------------- */
    Acquisition::Receivers<ValueType> receivers; 
    ValueType NXPerShot;
    IndexType numShotPerSuperShot;
    receivers.getModelPerShotSize(commAll, config, NXPerShot, numShotPerSuperShot);
    Acquisition::Coordinates<ValueType> modelCoordinates(config, 1, NXPerShot);
    Acquisition::Coordinates<ValueType> modelCoordinatesBig;

    if (config.get<bool>("useVariableGrid")) {
        CheckParameter::checkVariableGrid(config, commAll, modelCoordinates);
        for (int layer = 0; layer < modelCoordinates.getNumLayers(); layer++) {
            HOST_PRINT(commAll, "\n Number of gridpoints in layer: " << layer << " = " << modelCoordinates.getNGridpoints(layer));
        }
        auto numGridpointsRegular = config.get<IndexType>("NX") * config.get<IndexType>("NY") * config.get<IndexType>("NZ");
        HOST_PRINT(commAll, "\n Number of gripoints total: " << modelCoordinates.getNGridpoints());
        HOST_PRINT(commAll, "\n Percentage of gridpoints of the underlying regular grid given by NX*NY*NZ: " << (float)modelCoordinates.getNGridpoints() / numGridpointsRegular * 100 << "% \n\n");
    }

    /* ------------------------------------------------- */
    /* context and communicator for shot parallelisation */
    /* ------------------------------------------------- */
    /* execution context */
    hmemo::ContextPtr ctx = hmemo::Context::getContextPtr(); // default context, set by environment variable SCAI_CONTEXT

    IndexType shotDomain = Partitioning::getShotDomain(config, commAll); // will contain the domain to which this processor belongs

    // Build subsets of processors for the shots
    dmemo::CommunicatorPtr commShot = commAll->split(shotDomain);
    dmemo::CommunicatorPtr commInterShot = commAll->split(commShot->getRank());
    SCAI_DMEMO_TASK(commShot)

    /* --------------------------------------- */
    /* Distribution                            */
    /* --------------------------------------- */
    dmemo::DistributionPtr dist = nullptr;
    dmemo::DistributionPtr distBig = nullptr;
    IndexType configPartitioning = config.get<IndexType>("partitioning");
    if (useStreamConfig) {
        SCAI_ASSERT_ERROR(configPartitioning == 1, "partitioning != 1 when useStreamConfig"); 
    }
    switch (configPartitioning) {
    case 0:
    case 2:
    case 3:
        //Block distribution = starting distribution for graph partitioner
        dist = std::make_shared<dmemo::BlockDistribution>(modelCoordinates.getNGridpoints(), commShot);
        break;
    case 1:
        SCAI_ASSERT(!config.get<bool>("useVariableGrid"), "Grid distribution is not available for the variable grid");
        dist = Partitioning::gridPartition<ValueType>(config, commShot, NXPerShot);
        if (useStreamConfig) {
            modelCoordinatesBig.init(configBig);
            distBig = Partitioning::gridPartition<ValueType>(configBig, commShot);
        }
        break;
    default:
        COMMON_THROWEXCEPTION("unknown partitioning method = " << configPartitioning);
    }

    if (config.get<bool>("writeCoordinate") && shotDomain == 0) {
        // every shotdomain owns the same coordinates
        modelCoordinates.writeCoordinates(dist, ctx, config.get<std::string>("coordinateFilename"), config.get<IndexType>("FileFormat"));
    }

    /* --------------------------------------- */
    /* Factories                               */
    /* --------------------------------------- */
    ForwardSolver::Derivatives::Derivatives<ValueType>::DerivativesPtr derivatives(ForwardSolver::Derivatives::Factory<ValueType>::Create(dimension));
    
    Modelparameter::Modelparameter<ValueType>::ModelparameterPtr model(Modelparameter::Factory<ValueType>::Create(equationType));
    Wavefields::Wavefields<ValueType>::WavefieldPtr wavefields(Wavefields::Factory<ValueType>::Create(dimension, equationType));
    ForwardSolver::ForwardSolver<ValueType>::ForwardSolverPtr solver(ForwardSolver::Factory<ValueType>::Create(dimension, equationType));
    
    /* --------------------------------------- */
    /* Memory estimation                       */
    /* --------------------------------------- */
    HOST_PRINT(commAll, " ========== "<< dimension <<" "<< equationType <<" Memory Estimation: ===========\n\n")
    ValueType memDerivatives = derivatives->estimateMemory(config, dist, modelCoordinates);
    IndexType numShotDomains = config.get<IndexType>("NumShotDomains"); // total number of shot domains
    Common::checkNumShotDomains(numShotDomains, commAll);
    ValueType memWavefileds = wavefields->estimateMemory(dist, numRelaxationMechanisms);
    ValueType memModel = model->estimateMemory(dist);
    ValueType memSolver = solver->estimateMemory(config, dist, modelCoordinates);
    ValueType memTotal = memDerivatives + memWavefileds + memModel + memSolver;

    HOST_PRINT(commAll, " -  Derivative Matrices \t" << memDerivatives << " MB\n");
    HOST_PRINT(commAll, " -  Wavefield vectors \t\t" << memWavefileds << " MB\n");
    HOST_PRINT(commAll, " -  Model Vectors \t\t" << memModel << " MB\n");
    HOST_PRINT(commAll, " -  Boundary Condition Vectors \t" << memSolver << " MB\n");
    HOST_PRINT(commAll, "\n Memory Usage (total / per partition): \n " << memTotal << " / " << memTotal / dist->getNumPartitions() << " MB ");
    if (numShotDomains > 1)
        HOST_PRINT(commAll, "\n Total Memory Usage (" << numShotDomains << " shot domains): \n " << memTotal * numShotDomains << " MB  ");

    HOST_PRINT(commAll, "\n\n ===========================================================\n")

    /* --------------------------------------- */
    /* Call partitioner                        */
    /* --------------------------------------- */
    if (configPartitioning == 2) {
        SCAI_REGION("WAVE-Simulation.partitioningGEO")
        start_t = common::Walltime::get();
        dist = Partitioning::graphPartition(config, ctx, commShot, dist, *derivatives, modelCoordinates);
        end_t = common::Walltime::get();
        HOST_PRINT(commAll, "", "Finished Geographer graph partitioning in " << end_t - start_t << " sec.\n\n");
    }

    if (configPartitioning == 3) {
        SCAI_REGION("WAVE-Simulation.partitioningMetis")
        start_t = common::Walltime::get();
        dist = Partitioning::metisPartition(config, ctx, commShot, dist, *derivatives, modelCoordinates);
        end_t = common::Walltime::get();
        HOST_PRINT(commAll, "", "Finished ParMetis graph partitioning in " << end_t - start_t << " sec.\n\n");
    }

    bool writePartition = config.getAndCatch("writePartition", false);

    if (writePartition) {
        scai::lama::DenseVector<IndexType> partition(dist, commShot->getRank());
        IO::writeVector(partition, config.get<std::string>("partitionFilename") + std::to_string(shotDomain), config.get<IndexType>("fileFormat"));
    }

    /* --------------------------------------- */
    /* Calculate derivative matrizes           */
    /* --------------------------------------- */
    SCAI_REGION("WAVE-Simulation.initMatrices")
    start_t = common::Walltime::get();
    derivatives->init(dist, ctx, modelCoordinates, commShot);
    end_t = common::Walltime::get();
    HOST_PRINT(commAll, "", "Finished initializing matrices in " << end_t - start_t << " sec.\n\n");

    /* --------------------------------------- */
    /* Wavefields                              */
    /* --------------------------------------- */
    SCAI_REGION("WAVE-Simulation.initWavefields")
    start_t = common::Walltime::get();
    wavefields->init(ctx, dist, numRelaxationMechanisms);
    end_t = common::Walltime::get();
    HOST_PRINT(commAll, "", "Finished initializing wavefield in " << end_t - start_t << " sec.\n\n");

    /* --------------------------------------- */
    /* Acquisition geometry                    */
    /* --------------------------------------- */
    start_t = common::Walltime::get();
    
    Acquisition::Sources<ValueType> sources;
    
    std::vector<Acquisition::sourceSettings<ValueType>> sourceSettings; 
    std::vector<Acquisition::sourceSettings<ValueType>> sourceSettingsEncode; 
    std::vector<Acquisition::coordinate3D> cutCoordinates;
    ValueType shotIncr = config.getAndCatch("shotIncr", 0.0);
    sources.getAcquisitionSettings(config, shotIncr);
    if (useStreamConfig) {
        std::vector<Acquisition::sourceSettings<ValueType>> sourceSettingsBig;
        sourceSettingsBig = sources.getSourceSettings(); 
        Acquisition::getCutCoord(config, cutCoordinates, sourceSettingsBig, modelCoordinates, modelCoordinatesBig);
        Acquisition::getSettingsPerShot(sourceSettings, sourceSettingsBig, cutCoordinates);
        sources.setSourceSettings(sourceSettings); // for useSourceEncode
    } else {
        sourceSettings = sources.getSourceSettings(); 
    }
    CheckParameter::checkSources(sourceSettings, modelCoordinates, commAll);
    // calculate vector with unique shot numbers and get number of shots
    std::vector<IndexType> uniqueShotNos;
    std::vector<IndexType> uniqueShotNosEncode;
    Acquisition::calcuniqueShotNo(uniqueShotNos, sourceSettings);
    IndexType numshots = 1;
    IndexType useSourceEncode = config.getAndCatch("useSourceEncode", 0);
    sources.calcSourceSettingsEncode(commAll, config, seedtime);
    if (useSourceEncode == 0) {
        numshots = uniqueShotNos.size();
    } else {
        sourceSettingsEncode = sources.getSourceSettingsEncode();
        Acquisition::calcuniqueShotNo(uniqueShotNosEncode, sourceSettingsEncode);
        numshots = numShotDomains;
    }
    SCAI_ASSERT_ERROR(numshots >= numShotDomains, "numshots < numShotDomains");
    sources.writeShotIndIncr(commAll, config, uniqueShotNos);
    sources.writeSourceFC(commAll, config); 
    sources.writeSourceEncode(commAll, config); 
    Acquisition::writeCutCoordToFile(commAll, config, cutCoordinates, uniqueShotNos, NXPerShot);    
    
    if (config.get<IndexType>("useReceiversPerShot") == 0) {
        receivers.init(config, modelCoordinates, ctx, dist);
    }
    
    end_t = common::Walltime::get();
    HOST_PRINT(commAll, "", "Finished initializing Acquisition in " << end_t - start_t << " sec.\n\n");
    
    /* --------------------------------------- */
    /* Modelparameter                          */
    /* --------------------------------------- */
    SCAI_REGION("WAVE-Simulation.initModel")
    start_t = common::Walltime::get();
    Modelparameter::Modelparameter<ValueType>::ModelparameterPtr modelPerShot(Modelparameter::Factory<ValueType>::Create(equationType)); 
    if (!useStreamConfig) {
        model->init(config, ctx, dist, modelCoordinates);
    } else { 
        model->init(configBig, ctx, distBig, modelCoordinatesBig);   
    }
    end_t = common::Walltime::get();
    HOST_PRINT(commAll, "", "Finished initializing model in " << end_t - start_t << " sec.\n\n");

    /* --------------------------------------- */
    /* Forward solver                          */
    /* --------------------------------------- */
    SCAI_REGION("WAVE-Simulation.initForwardSolver")
    start_t = common::Walltime::get();
    ValueType DT = config.get<ValueType>("DT");
    IndexType tStepEnd = Common::time2index(config.get<ValueType>("T"), DT);
    if (!useStreamConfig) {
        solver->initForwardSolver(config, *derivatives, *wavefields, *model, modelCoordinates, ctx, DT);
    }
    end_t = common::Walltime::get();
    HOST_PRINT(commAll, "", "Finished initializing forward solver in " << end_t - start_t << " sec.\n\n");
    
    /* --------------------------------------- */
    /* Hilbert handler                         */
    /* --------------------------------------- */
    IndexType snapType = config.get<IndexType>("snapType");
    Hilbert::HilbertFFT<ValueType> hilbertHandlerTime;
    IndexType decomposition = config.getAndCatch("decomposeWavefieldType", 0);
    if (decomposition != 0) {
        IndexType kernelSize = Common::calcNextPowTwo<ValueType>(tStepEnd - 1);  
        hilbertHandlerTime.setCoefficientLength(kernelSize);
        hilbertHandlerTime.calcHilbertCoefficient(); 
        snapType = decomposition + 3;
    }
    
    /* --------------------------------------- */
    /* Wavefield additional                        */
    /* --------------------------------------- */
    Wavefields::Wavefields<ValueType>::WavefieldPtr wavefieldsTemp = Wavefields::Factory<ValueType>::Create(dimension, equationType);
    
    /* --------------------------------------- */
    /* Preparation                             */
    /* --------------------------------------- */
    if (!useStreamConfig) {            
        model->prepareForModelling(modelCoordinates, ctx, dist, commShot); 
        solver->prepareForModelling(*model, DT);
    }
    
    IndexType maxcount = 1;    
    std::vector<IndexType> shotHistory(numshots, 0);
    std::shared_ptr<const dmemo::BlockDistribution> shotDist;
    if (config.getAndCatch("useRandomSource", 0) != 0) {  
        shotDist = dmemo::blockDistribution(numShotDomains, commInterShot);
    } else {
        shotDist = dmemo::blockDistribution(numshots, commInterShot);
    }
    IndexType numRand = numshots / numShotDomains;  
    if (decomposition != 0) {
        numRand = 2;
    }
    
    double tInit = common::Walltime::get();
    HOST_PRINT(commAll, "\nFinished all initialization in " << tInit - globalStart_t << " sec.\n");
        
    /* --------------------------------------- */
    /* Loop over shots                         */
    /* --------------------------------------- */
    for (IndexType randInd = 0; randInd < numRand; randInd++) { 
        sources.calcUniqueShotInds(commAll, config, shotHistory, maxcount, seedtime);
        std::vector<IndexType> uniqueShotInds = sources.getUniqueShotInds();
        IndexType shotNumber;
        IndexType shotIndTrue = 0;
        for (IndexType shotInd = shotDist->lb(); shotInd < shotDist->ub(); shotInd++) {
            SCAI_REGION("WAVE-Simulation.shotLoop")
            shotIndTrue = uniqueShotInds[shotInd];
            
            /* Update Source */
            std::vector<Acquisition::sourceSettings<ValueType>> sourceSettingsShot;
            if (useSourceEncode == 0) {
                shotNumber = uniqueShotNos[shotIndTrue];
                Acquisition::createSettingsForShot(sourceSettingsShot, sourceSettings, shotNumber);
            } else {
                shotNumber = uniqueShotNosEncode[shotIndTrue];
                Acquisition::createSettingsForShot(sourceSettingsShot, sourceSettingsEncode, shotNumber);
            }                    
            sources.init(sourceSettingsShot, config, modelCoordinates, ctx, dist);

            if (!useStreamConfig) {
                CheckParameter::checkNumericalArtefactsAndInstabilities<ValueType>(config, sourceSettingsShot, *model, modelCoordinates, shotNumber);
            } else {
                HOST_PRINT(commShot, "Shot number " << shotNumber << " (" << "domain " << commInterShot->getRank() << ", index " << shotIndTrue + 1 << " of " << numshots << "): Switch to model subset\n");
                IndexType shotIndPerShot = shotIndTrue;
                if (useSourceEncode == 3) {
                    Acquisition::getuniqueShotInd(shotIndPerShot, sourceSettingsEncode, shotNumber);
                }
                model->getModelPerShot(*modelPerShot, dist, modelCoordinates, modelCoordinatesBig, cutCoordinates.at(shotIndPerShot));
                modelPerShot->prepareForModelling(modelCoordinates, ctx, dist, commShot); 
                modelPerShot->write((config.get<std::string>("ModelFilename") + ".shot_" + std::to_string(shotNumber)), config.get<IndexType>("FileFormat"));
                solver->initForwardSolver(config, *derivatives, *wavefields, *modelPerShot, modelCoordinates, ctx, DT);
                solver->prepareForModelling(*modelPerShot, DT);
                
                CheckParameter::checkNumericalArtefactsAndInstabilities<ValueType>(config, sourceSettingsShot, *modelPerShot, modelCoordinates, shotNumber);
            }

            if (randInd == 1 && decomposition != 0) {
                lama::DenseMatrix<ValueType> sourcesignalHilbert = sources.getsourcesignal();
                hilbertHandlerTime.hilbert(sourcesignalHilbert);
                sources.setsourcesignal(sourcesignalHilbert);
            }
            bool writeSource = config.getAndCatch("writeSource", false);
            if (writeSource) {
                if (randInd == 1 && decomposition != 0) {
                    sources.getSeismogramHandler().write(config.get<IndexType>("SeismogramFormat"), config.get<std::string>("writeSourceFilename") + ".shot_" + std::to_string(shotNumber) + ".Hilbert", modelCoordinates);
                } else {
                    sources.getSeismogramHandler().write(config.get<IndexType>("SeismogramFormat"), config.get<std::string>("writeSourceFilename") + ".shot_" + std::to_string(shotNumber), modelCoordinates);
                }
            }

            if (config.get<IndexType>("useReceiversPerShot") != 0) {
                receivers.init(config, modelCoordinates, ctx, dist, shotNumber, sourceSettingsEncode);
            }

            if (randInd == 1 && decomposition != 0) {
                HOST_PRINT(commShot, "Start time stepping for shot number " << shotNumber << " (" << "domain " << shotDomain << ", index " << shotIndTrue + 1 << " of " << numshots << ") Hilbert\n", "\nTotal Number of time steps: " << tStepEnd << "\n");
            } else {
                HOST_PRINT(commShot, "Start time stepping for shot number " << shotNumber << " (" << "domain " << shotDomain << ", index " << shotIndTrue + 1 << " of " << numshots << ")\n", "\nTotal Number of time steps: " << tStepEnd << "\n");
            }
            
            start_t = common::Walltime::get();
            wavefields->resetWavefields();

            double start_t2 = 0.0, end_t2 = 0.0;

            /* --------------------------------------- */
            /* Loop over time steps                    */
            /* --------------------------------------- */
            ValueType DTinv = 1 / config.get<ValueType>("DT");
            lama::DenseVector<ValueType> compensation;
            if (!useStreamConfig) {
                if (config.getAndCatch("compensation", 0))
                    compensation = model->getCompensation(DT, 1);
                for (IndexType tStep = 0; tStep < tStepEnd; tStep++) {

                    SCAI_REGION("WAVE-Simulation.timeLoop")
                    if ((tStep - 1) % 100 == 0) {
                        start_t2 = common::Walltime::get();
                    }
                    *wavefieldsTemp = *wavefields;

                    solver->run(receivers, sources, *model, *wavefields, *derivatives, tStep);
                    
                    if (config.getAndCatch("compensation", 0))
                        *wavefields *= compensation;
                    
                    if (randInd == 0 && decomposition != 0) {
                        *wavefieldsTemp -= *wavefields;
                        *wavefieldsTemp *= -DTinv;
                        wavefields->decompose(decomposition, *wavefieldsTemp, *derivatives);
                    }

                    if (tStep % 100 == 0 && tStep != 0) {
                        end_t2 = common::Walltime::get();
                        HOST_PRINT(commShot, "", "Calculated " << tStep << " time steps" << " in shot  " << shotNumber << " at t = " << end_t2 - globalStart_t << "\nLast 100 timesteps calculated in " << end_t2 - start_t2 << " sec. - Estimated runtime (Simulation/total): " << (int)((tStepEnd / 100) * (end_t2 - start_t2)) << " / " << (int)((tStepEnd / 100) * (end_t2 - start_t2) + tInit) << " sec.\n\n");
                    }

                    if (snapType > 0 && tStep >= Common::time2index(config.get<ValueType>("tFirstSnapshot"), DT) && tStep <= Common::time2index(config.get<ValueType>("tlastSnapshot"), DT) && (tStep - Common::time2index(config.get<ValueType>("tFirstSnapshot"), DT)) % Common::time2index(config.get<ValueType>("tincSnapshot"), DT) == 0) {
                        if (randInd == 1 && decomposition != 0) {
                            wavefields->write(1, config.get<std::string>("WavefieldFileName") + ".shot_" + std::to_string(shotNumber) + ".HilbertT", tStep, *derivatives, *model, config.get<IndexType>("FileFormat"));
                            wavefields->write(2, config.get<std::string>("WavefieldFileName") + ".shot_" + std::to_string(shotNumber) + ".HilbertT", tStep, *derivatives, *model, config.get<IndexType>("FileFormat"));
                        } else {
                            wavefields->write(snapType, config.get<std::string>("WavefieldFileName") + ".shot_" + std::to_string(shotNumber), tStep, *derivatives, *model, config.get<IndexType>("FileFormat"));
                        }
                    }
                }
            } else {                
                if (config.getAndCatch("compensation", 0))
                    compensation = modelPerShot->getCompensation(DT, 1);
                for (IndexType tStep = 0; tStep < tStepEnd; tStep++) {

                    SCAI_REGION("WAVE-Simulation.timeLoop")
                    if ((tStep - 1) % 100 == 0) {
                        start_t2 = common::Walltime::get();
                    }
                    *wavefieldsTemp = *wavefields;

                    solver->run(receivers, sources, *modelPerShot, *wavefields, *derivatives, tStep);
                    
                    if (config.getAndCatch("compensation", 0))
                        *wavefields *= compensation;
                    
                    if (randInd == 0 && decomposition != 0) {
                        *wavefieldsTemp -= *wavefields;
                        *wavefieldsTemp *= -DTinv;
                        wavefields->decompose(decomposition, *wavefieldsTemp, *derivatives);
                    }

                    if (tStep % 100 == 0 && tStep != 0) {
                        end_t2 = common::Walltime::get();
                        HOST_PRINT(commShot, "", "Calculated " << tStep << " time steps" << " in shot  " << shotNumber << " at t = " << end_t2 - globalStart_t << "\nLast 100 timesteps calculated in " << end_t2 - start_t2 << " sec. - Estimated runtime (Simulation/total): " << (int)((tStepEnd / 100) * (end_t2 - start_t2)) << " / " << (int)((tStepEnd / 100) * (end_t2 - start_t2) + tInit) << " sec.\n\n");
                    }

                    if (snapType > 0 && tStep >= Common::time2index(config.get<ValueType>("tFirstSnapshot"), DT) && tStep <= Common::time2index(config.get<ValueType>("tlastSnapshot"), DT) && (tStep - Common::time2index(config.get<ValueType>("tFirstSnapshot"), DT)) % Common::time2index(config.get<ValueType>("tincSnapshot"), DT) == 0) {
                        if (randInd == 1 && decomposition != 0) {
                            wavefields->write(1, config.get<std::string>("WavefieldFileName") + ".shot_" + std::to_string(shotNumber) + ".HilbertT", tStep, *derivatives, *modelPerShot, config.get<IndexType>("FileFormat"));
                            wavefields->write(2, config.get<std::string>("WavefieldFileName") + ".shot_" + std::to_string(shotNumber) + ".HilbertT", tStep, *derivatives, *modelPerShot, config.get<IndexType>("FileFormat"));
                        } else {
                            wavefields->write(snapType, config.get<std::string>("WavefieldFileName") + ".shot_" + std::to_string(shotNumber), tStep, *derivatives, *modelPerShot, config.get<IndexType>("FileFormat"));
                        }
                    }
                }
            }
            solver->resetCPML();
            end_t = common::Walltime::get();
            HOST_PRINT(commShot, "Finished time stepping for shot number: " << shotNumber << " in " << end_t - start_t << " sec.\n", "");
            
            // check wavefield and seismogram for NaNs or infinite values
            SCAI_ASSERT_ERROR(commShot->all(wavefields->isFinite(dist)) && commShot->all(receivers.getSeismogramHandler().isFinite()),"Infinite or NaN value in seismogram or/and velocity wavefield!") // if all processors return isfinite=true, everything is finite
            
            if (config.get<IndexType>("normalizeTraces") == 3) {
                receivers.getSeismogramHandler().setFrequencyAGC(config.get<ValueType>("CenterFrequencyCPML"));
                receivers.getSeismogramHandler().calcInverseAGC();
                receivers.getSeismogramHandler().write(5, config.get<std::string>("SeismogramFilename") + ".shot_" + std::to_string(shotNumber), modelCoordinates);
            }
            receivers.getSeismogramHandler().normalize(config.get<IndexType>("normalizeTraces"));

            if (randInd == 1 && decomposition != 0) { 
                receivers.getSeismogramHandler().write(config.get<IndexType>("SeismogramFormat"), config.get<std::string>("SeismogramFilename") + ".shot_" + std::to_string(shotNumber) + ".Hilbert", modelCoordinates);
            } else {
                receivers.getSeismogramHandler().write(config.get<IndexType>("SeismogramFormat"), config.get<std::string>("SeismogramFilename") + ".shot_" + std::to_string(shotNumber), modelCoordinates);
                receivers.decode(config, config.get<std::string>("SeismogramFilename"), shotNumber, sourceSettingsEncode, 1);
                receivers.writeReceiverMark(config, shotNumber);
            }                
        }
                   
        if (config.getAndCatch("useRandomSource", 0) == 0 && decomposition == 0) 
            break;
    }
    globalEnd_t = common::Walltime::get();

    commAll->synchronize();

    HOST_PRINT(commAll, "\nTotal runtime of WAVE-Simulation: " << globalEnd_t - globalStart_t << " sec.\nWAVE-Simulation finished!\n\n");
    return 0;
}
