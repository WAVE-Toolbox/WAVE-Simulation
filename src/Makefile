CC=g++
CFLAGS=-std=c++11 -Wall -fopenmp
IFLAGS=-I${SCAI_ROOT}/include -L${SCAI_ROOT}/lib -DSCAI_ASSERT_LEVEL_ERROR -DSCAI_LOG_LEVEL_ERROR -DSCAI_TRACE_OFF
LDFLAGS=-lscai_lama -lscai_solver -lscai_dmemo -lscai_sparsekernel -lscai_utilskernel -lscai_blaskernel -lscai_kregistry -lscai_hmemo -lscai_tasking -lscai_tracing -lscai_logging -lscai_common

install: FDSimulation

FDSimulation: FDSimulation.cpp
	$(CC) $(CFLAGS) -o ../bin/FDSimulation FDSimulation.cpp $(IFLAGS) $(LDFLAGS)