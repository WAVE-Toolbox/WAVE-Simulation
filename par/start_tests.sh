#!/bin/bash

# Setup test environment
export OMP_NUM_THREADS=2
export NUM_MPI_PROCESSES=4
export SCAI_UNSUPPORTED=IGNORE

SOFI_EXE="./../build/bin/SOFI"
UNITTEST_EXE="./../build/bin/Test_unit"
INTEGRATIONTEST_EXE="./../build/bin/Test_integration"

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
exit
if [ "$?" != "0" ]; then
	echo "Test failed ! "
	return
fi

mpirun -np ${NUM_MPI_PROCESSES} ${SOFI_EXE} ci/configuration_ci.2D.acoustic.txt
${INTEGRATIONTEST_EXE} ci/configuration_ci.2D.acoustic.txt
if [ "$?" != "0" ]; then
	echo "Test failed ! "
	return
fi

mpirun -np ${NUM_MPI_PROCESSES} ${SOFI_EXE} ci/configuration_ci.3D.acoustic.txt
${INTEGRATIONTEST_EXE} ci/configuration_ci.3D.acoustic.txt
if [ "$?" != "0" ]; then
	echo "Test failed ! "
	return
fi

mpirun -np ${NUM_MPI_PROCESSES} ${SOFI_EXE} ci/configuration_ci.2D.elastic.txt
${INTEGRATIONTEST_EXE} ci/configuration_ci.2D.elastic.txt
if [ "$?" != "0" ]; then
	echo "Test failed ! "
	return
fi

mpirun -np ${NUM_MPI_PROCESSES} ${SOFI_EXE} ci/configuration_ci.3D.elastic.txt
${INTEGRATIONTEST_EXE} ci/configuration_ci.3D.elastic.txt
if [ "$?" != "0" ]; then
	echo "Test failed ! "
	return
fi

mpirun -np ${NUM_MPI_PROCESSES} ${SOFI_EXE} ci/configuration_ci.2D.visco.txt
${INTEGRATIONTEST_EXE} ci/configuration_ci.2D.visco.txt
if [ "$?" != "0" ]; then
	echo "Test failed ! "
	return
fi

mpirun -np ${NUM_MPI_PROCESSES} ${SOFI_EXE} ci/configuration_ci.3D.visco.txt
${INTEGRATIONTEST_EXE} ci/configuration_ci.3D.visco.txt
if [ "$?" != "0" ]; then
	echo "Test failed ! "
	return
fi
