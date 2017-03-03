#!/bin/bash

# Setup test environment
export OMP_NUM_THREADS=2
export NUM_MPI_PROCESSES=2
export SCAI_UNSUPPORTED=IGNORE

# Run unit tests
./../bin/Tests/Test_unit
if [ "$?" != "0" ]; then
	echo "Test failed ! "
	exit 1
fi

# Run simulation tests
rm -rf ci/*.ci.*

mpirun -np ${NUM_MPI_PROCESSES} ./../bin/SOFI ci/configuration_ci.2D.acoustic.txt
./../bin/Tests/Test_CompareSeismogram ci/configuration_ci.2D.acoustic.txt
if [ "$?" != "0" ]; then
	echo "Test failed ! "
	return
fi

mpirun -np ${NUM_MPI_PROCESSES} ./../bin/SOFI ci/configuration_ci.3D.acoustic.txt
./../bin/Tests/Test_CompareSeismogram ci/configuration_ci.3D.acoustic.txt
if [ "$?" != "0" ]; then
	echo "Test failed ! "
	return
fi

mpirun -np ${NUM_MPI_PROCESSES} ./../bin/SOFI ci/configuration_ci.2D.elastic.txt
./../bin/Tests/Test_CompareSeismogram ci/configuration_ci.2D.elastic.txt
if [ "$?" != "0" ]; then
	echo "Test failed ! "
	return
fi

mpirun -np ${NUM_MPI_PROCESSES} ./../bin/SOFI ci/configuration_ci.3D.elastic.txt
./../bin/Tests/Test_CompareSeismogram ci/configuration_ci.3D.elastic.txt
if [ "$?" != "0" ]; then
	echo "Test failed ! "
	return
fi

mpirun -np ${NUM_MPI_PROCESSES} ./../bin/SOFI ci/configuration_ci.2D.visco.txt
./../bin/Tests/Test_CompareSeismogram ci/configuration_ci.2D.visco.txt
if [ "$?" != "0" ]; then
	echo "Test failed ! "
	return
fi

mpirun -np ${NUM_MPI_PROCESSES} ./../bin/SOFI ci/configuration_ci.3D.visco.txt
./../bin/Tests/Test_CompareSeismogram ci/configuration_ci.3D.visco.txt
if [ "$?" != "0" ]; then
	echo "Test failed ! "
	return
fi