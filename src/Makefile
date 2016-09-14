install: FDSimulation

FDSimulation: FDSimulation.cpp
	g++ -std=c++11 -o ../bin/FDSimulation FDSimulation.cpp -I${SCAI_ROOT}/include -L${SCAI_ROOT}/lib -lscai_lama -DSCAI_ASSERT_LEVEL_ERROR -DSCAI_LOG_LEVEL_ERROR -DSCAI_TRACE_OFF -lscai_solver -lscai_lama -lscai_dmemo -lscai_sparsekernel -lscai_utilskernel -lscai_blaskernel -lscai_kregistry -lscai_hmemo -lscai_tasking -lscai_tracing -lscai_logging -lscai_common -Wall -fopenmp