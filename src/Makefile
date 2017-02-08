CC=g++
CFLAGS=-std=c++11 -fopenmp -Wall -Wextra -Werror
IFLAGS=-I${SCAI_ROOT}/include -DSCAI_ASSERT_LEVEL_DEBUG -DSCAI_LOG_LEVEL_DEBUG -DSCAI_TRACE_OFF
LDFLAGS= -L${SCAI_ROOT}/lib -lscai_lama -lscai_solver -lscai_dmemo -lscai_sparsekernel -lscai_utilskernel -lscai_blaskernel -lscai_kregistry -lscai_hmemo -lscai_tasking -lscai_tracing -lscai_logging -lscai_common

# Flags for googletest framework
CODECOVERAGE=OFF
IFLAGSTEST = -I${GTEST_DIR}/include -L${GTEST_DIR}/
LDFLAGSTEST = -lgtest -lgtest_main
CFLAGSTEST = -isystem -pthread 
ifeq ($(CODECOVERAGE),ON)
	CFLAGSTEST += -fprofile-arcs -ftest-coverage -fPIC -O0
endif



install: clean SOFI

.cpp.o:
	$(CC) $(CFLAGS) -c  $< -o $@ $(IFLAGS)

subdirs := $(wildcard */)
SOFI_SCR := $(wildcard $(addsuffix */*.cpp,$(subdirs))) $(wildcard **/*.cpp)
SOFIUnitTest_SCR := $(wildcard $(addsuffix */*UnitTest.cpp,$(subdirs))) $(wildcard */*UnitTest.cpp)
SOFI_SCR := $(filter-out $(SOFIUnitTest_SCR) $(wildcard **/*Test_CompareSeismogram.cpp), $(SOFI_SCR))

SOFI_OBJ = $(SOFI_SCR:%.cpp=%.o)

SOFI: SOFI.cpp $(SOFI_OBJ)
	$(CC) $(CFLAGS) $(SOFI_OBJ) -o ../bin/SOFI SOFI.cpp $(IFLAGS) $(LDFLAGS)


SOFIexplicit: SOFI3D SOFI2D

SOFI3D: SOFI3Dacoustic SOFI3Delastic SOFI3Dvisco

SOFI3Dacoustic: SOFIexplicit/SOFI3Dacoustic.cpp
	$(CC) $(CFLAGS) -o ../bin/SOFI3Dacoustic SOFIexplicit/SOFI3Dacoustic.cpp $(IFLAGS) $(LDFLAGS)

SOFI3Delastic: SOFIexplicit/SOFI3Delastic.cpp
	$(CC) $(CFLAGS) -o ../bin/SOFI3Delastic SOFIexplicit/SOFI3Delastic.cpp $(IFLAGS) $(LDFLAGS)

SOFI3Dvisco: SOFIexplicit/SOFI3Dvisco.cpp
	$(CC) $(CFLAGS) -o ../bin/SOFI3Dvisco SOFIexplicit/SOFI3Dvisco.cpp $(IFLAGS) $(LDFLAGS)


SOFI2D: SOFI2Dacoustic  SOFI2Delastic SOFI2Dvisco

SOFI2Dacoustic: SOFIexplicit/SOFI2Dacoustic.cpp
	$(CC) $(CFLAGS) -o ../bin/SOFI2Dacoustic SOFIexplicit/SOFI2Dacoustic.cpp $(IFLAGS) $(LDFLAGS)

SOFI2Delastic: SOFIexplicit/SOFI2Delastic.cpp
	$(CC) $(CFLAGS) -o ../bin/SOFI2Delastic SOFIexplicit/SOFI2Delastic.cpp $(IFLAGS) $(LDFLAGS)

SOFI2Dvisco: SOFIexplicit/SOFI2Dvisco.cpp
	$(CC) $(CFLAGS) -o ../bin/SOFI2Dvisco SOFIexplicit/SOFI2Dvisco.cpp $(IFLAGS) $(LDFLAGS)



tests: integrationTest unitTest

UNITTESTPATH := $(shell find $(shell pwd) -name '*UnitTest.cpp')

unitTest: $(UNITTESTPATH)
	$(CC) $(CFLAGS) $(CFLAGSTEST) -o ../bin/Tests/Test_unit $(UNITTESTPATH) $(IFLAGS) $(IFLAGSTEST) $(LDFLAGS) $(LDFLAGSTEST)

integrationTest: Tests/Test_CompareSeismogram.cpp
	$(CC) $(CFLAGS) -o ../bin/Tests/Test_CompareSeismogram Tests/Test_CompareSeismogram.cpp $(IFLAGS) $(LDFLAGS)



clean:
	rm -rf ../bin/*
	mkdir ../bin/Tests
