CC=g++
CFLAGS=-std=c++11 -fopenmp -Wall -Wextra -Werror
IFLAGS=-I${SCAI_ROOT}/include -L${SCAI_ROOT}/lib -DSCAI_ASSERT_LEVEL_ERROR -DSCAI_LOG_LEVEL_ERROR -DSCAI_TRACE_OFF
LDFLAGS=-lscai_lama -lscai_solver -lscai_dmemo -lscai_sparsekernel -lscai_utilskernel -lscai_blaskernel -lscai_kregistry -lscai_hmemo -lscai_tasking -lscai_tracing -lscai_logging -lscai_common


install: clean SOFI tests

SOFI: SOFI3D SOFI2D


SOFI3D: SOFI3Dacoustic SOFI3Delastic SOFI3Dvisco

SOFI3Dacoustic: SOFI3Dacoustic.cpp
	$(CC) $(CFLAGS) -o ../bin/SOFI3Dacoustic SOFI3Dacoustic.cpp $(IFLAGS) $(LDFLAGS)

SOFI3Delastic: SOFI3Delastic.cpp
	$(CC) $(CFLAGS) -o ../bin/SOFI3Delastic SOFI3Delastic.cpp $(IFLAGS) $(LDFLAGS)

SOFI3Dvisco: SOFI3Dvisco.cpp
	$(CC) $(CFLAGS) -o ../bin/SOFI3Dvisco SOFI3Dvisco.cpp $(IFLAGS) $(LDFLAGS)


SOFI2D: SOFI2Dacoustic  SOFI2Delastic SOFI2Dvisco

SOFI2Dacoustic: SOFI2Dacoustic.cpp
	$(CC) $(CFLAGS) -o ../bin/SOFI2Dacoustic SOFI2Dacoustic.cpp $(IFLAGS) $(LDFLAGS)

SOFI2Delastic: SOFI2Delastic.cpp
	$(CC) $(CFLAGS) -o ../bin/SOFI2Delastic SOFI2Delastic.cpp $(IFLAGS) $(LDFLAGS)

SOFI2Dvisco: SOFI2Dvisco.cpp
	$(CC) $(CFLAGS) -o ../bin/SOFI2Dvisco SOFI2Dvisco.cpp $(IFLAGS) $(LDFLAGS)


tests: Test_SOFI3Dacoustic Test_SOFI3Delastic Test_SOFI3Dvisco

Test_SOFI3Dacoustic: Tests/Test_SOFI3Dacoustic.cpp
	$(CC) $(CFLAGS) -o ../bin/Tests/Test_SOFI3Dacoustic Tests/Test_SOFI3Dacoustic.cpp $(IFLAGS) $(LDFLAGS)

Test_SOFI3Delastic: Tests/Test_SOFI3Delastic.cpp
	$(CC) $(CFLAGS) -o ../bin/Tests/Test_SOFI3Delastic Tests/Test_SOFI3Delastic.cpp $(IFLAGS) $(LDFLAGS)

Test_SOFI3Dvisco: Tests/Test_SOFI3Dvisco.cpp
	$(CC) $(CFLAGS) -o ../bin/Tests/Test_SOFI3Dvisco Tests/Test_SOFI3Dvisco.cpp $(IFLAGS) $(LDFLAGS)


clean:
	rm -rf ../bin/*
	mkdir ../bin/Tests
