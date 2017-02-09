CC=g++
CFLAGS=-std=c++11 -fopenmp -Wall -Wextra -Werror
IFLAGS=-I${SCAI_ROOT}/include -DSCAI_ASSERT_LEVEL_DEBUG -DSCAI_LOG_LEVEL_DEBUG -DSCAI_TRACE_OFF
LDFLAGS= -L${SCAI_ROOT}/lib -lscai_lama -lscai_solver -lscai_dmemo -lscai_sparsekernel -lscai_utilskernel -lscai_blaskernel -lscai_kregistry -lscai_hmemo -lscai_tasking -lscai_tracing -lscai_logging -lscai_common

# Flags for googletest framework
CODECOVERAGE=OFF
IFLAGSTEST = -I${GTEST_DIR}/include
LDFLAGSTEST = -L${GTEST_DIR}/ -lgtest -lgtest_main
CFLAGSTEST = -isystem -pthread 


all: SOFI tests

install: SOFI

tests: integrationTest unitTest

.cpp.o:
	$(CC) $(CFLAGS) -c  $< -o $@ $(IFLAGS)

SUBDIRS := $(wildcard */)
SOFI_LIB_CPP := $(wildcard $(addsuffix */*.cpp,$(SUBDIRS))) $(wildcard **/*.cpp)
SOFI_MAIN_CPP := SOFI.cpp

SOFI_UNITTEST_SCR := $(wildcard $(addsuffix */*UnitTest.cpp,$(SUBDIRS))) $(wildcard */*UnitTest.cpp)
SOFI_UNITTEST_OBJ = $(SOFI_UNITTEST_SCR:%.cpp=%.o)

SOFI_LIB_SCR := $(filter-out $(SOFI_UNITTEST_SCR) $(wildcard **/*Test_CompareSeismogram.cpp), $(SOFI_LIB_CPP))
SOFI_LIB_OBJ = $(SOFI_LIB_SCR:%.cpp=%.o)
SOFI_MAIN_OBJ = $(SOFI_MAIN_CPP:%.cpp=%.o)

SOFI: $(SOFI_MAIN_OBJ) $(SOFI_LIB_OBJ)
	$(CC) $(CFLAGS) $(SOFI_MAIN_OBJ) $(SOFI_LIB_OBJ) -o ../bin/SOFI $(IFLAGS) $(LDFLAGS)

$(SOFI_UNITTEST_OBJ): CFLAGS+=$(CFLAGSTEST)
$(SOFI_UNITTEST_OBJ): LDFLAGS +=$(LDFLAGSTEST)
$(SOFI_UNITTEST_OBJ): IFLAGS +=$(IFLAGSTEST)

unitTest: $(SOFI_UNITTEST_OBJ) $(SOFI_LIB_OBJ)
	$(CC) $(CFLAGS) $(CFLAGSTEST) $(SOFI_UNITTEST_OBJ) $(SOFI_LIB_OBJ) -o ../bin/Tests/Test_unit $(IFLAGS) $(IFLAGSTEST) $(LDFLAGS) $(LDFLAGSTEST)

integrationTest: Tests/Test_CompareSeismogram.cpp
	mkdir -p ../bin/Tests/
	$(CC) $(CFLAGS) $(SOFI_LIB_OBJ) -o ../bin/Tests/Test_CompareSeismogram Tests/Test_CompareSeismogram.cpp $(IFLAGS) $(LDFLAGS)

clean:
	rm -rf $(SOFI_LIB_OBJ)
	rm -rf $(SOFI_MAIN_OBJ)
	rm -rf $(SOFI_UNITTEST_OBJ)
	rm -rf ../bin/*
	mkdir ../bin/Tests
