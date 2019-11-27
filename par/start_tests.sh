#!/bin/bash

# Setup test environment
export OMP_NUM_THREADS=1
export NUM_MPI_PROCESSES=4
export SCAI_UNSUPPORTED=IGNORE
export SCAI_TRACE=OFF

BINDIR="./../build/bin"

SOFI_EXE="${BINDIR}/Simulation"
UNITTEST_EXE="${BINDIR}/Test_unit"
INTEGRATIONTEST_EXE="${BINDIR}/Test_integration"
MODEL_EXE="${BINDIR}//tools/TwoLayer"

export SCAI_LOG=OFF

# Run unit tests
${UNITTEST_EXE}
if [ "$?" != "0" ]; then
	echo "Test failed ! "
	exit 1
fi

# Run simulation tests
rm -rf ci/*.ci.*

mpirun -np ${NUM_MPI_PROCESSES} ${SOFI_EXE} ci/configuration_ci.2D.sh.txt
${INTEGRATIONTEST_EXE} ci/configuration_ci.2D.sh.txt
if [ "$?" != "0" ]; then
	echo "Test failed ! "
	exit
fi

mpirun -np ${NUM_MPI_PROCESSES} ${SOFI_EXE} ci/configuration_ci.2D.acoustic.txt
${INTEGRATIONTEST_EXE} ci/configuration_ci.2D.acoustic.txt
if [ "$?" != "0" ]; then
	echo "Test failed ! "
	exit
fi


mpirun -np ${NUM_MPI_PROCESSES} ${SOFI_EXE} ci/configuration_ci.2D.acoustic.varGrid.txt
${INTEGRATIONTEST_EXE} ci/configuration_ci.2D.acoustic.varGrid.txt
if [ "$?" != "0" ]; then
	echo "Test failed ! "
	exit
fi


mpirun -np ${NUM_MPI_PROCESSES} ${SOFI_EXE} ci/configuration_ci.3D.acoustic.varGrid.txt
${INTEGRATIONTEST_EXE} ci/configuration_ci.3D.acoustic.varGrid.txt
if [ "$?" != "0" ]; then
	echo "Test failed ! "
	exit
fi


${MODEL_EXE} ci/configuration_ci.3D.acoustic.txt || { echo "Model creation failed !" ; exit 1; } 
mpirun -np ${NUM_MPI_PROCESSES} ${SOFI_EXE} ci/configuration_ci.3D.acoustic.txt
${INTEGRATIONTEST_EXE} ci/configuration_ci.3D.acoustic.txt
if [ "$?" != "0" ]; then
	echo "Test failed ! "
	exit
fi

${MODEL_EXE} ci/configuration_ci.2D.elastic.txt || { echo "Model creation failed !"; exit 1; }
mpirun -np ${NUM_MPI_PROCESSES} ${SOFI_EXE} ci/configuration_ci.2D.elastic.txt
${INTEGRATIONTEST_EXE} ci/configuration_ci.2D.elastic.txt
if [ "$?" != "0" ]; then
	echo "Test failed ! "
	exit
fi

${MODEL_EXE} ci/configuration_ci.3D.elastic.txt || { echo "Model creation failed !"; exit 1; }
mpirun -np ${NUM_MPI_PROCESSES} ${SOFI_EXE} ci/configuration_ci.3D.elastic.txt
${INTEGRATIONTEST_EXE} ci/configuration_ci.3D.elastic.txt
if [ "$?" != "0" ]; then
	echo "Test failed ! "
	exit
fi

${MODEL_EXE} ci/configuration_ci.2D.visco.txt || { echo "Model creation failed !"; exit 1; }
mpirun -np ${NUM_MPI_PROCESSES} ${SOFI_EXE} ci/configuration_ci.2D.visco.txt
${INTEGRATIONTEST_EXE} ci/configuration_ci.2D.visco.txt
if [ "$?" != "0" ]; then
	echo "Test failed ! "
	exit
fi


${MODEL_EXE} ci/configuration_ci.3D.visco.txt || { echo "Model creation failed !"; exit 1; }
mpirun -np ${NUM_MPI_PROCESSES} ${SOFI_EXE} ci/configuration_ci.3D.visco.txt
${INTEGRATIONTEST_EXE} ci/configuration_ci.3D.visco.txt
if [ "$?" != "0" ]; then
	echo "Test failed ! "
	exit
fi

echo "All tests passed ! "
