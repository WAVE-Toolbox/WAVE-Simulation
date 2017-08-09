CC=g++
CFLAGS=-std=c++11 -fopenmp -Wall -Wextra -Werror
IFLAGS=-I${SCAI_ROOT}/include -DSCAI_ASSERT_LEVEL_DEBUG -DSCAI_LOG_LEVEL_DEBUG -DSCAI_TRACE_OFF
LDFLAGS= -L${SCAI_ROOT}/lib -lscai_lama -lscai_solver -lscai_dmemo -lscai_sparsekernel -lscai_utilskernel -lscai_blaskernel -lscai_kregistry -lscai_hmemo -lscai_tasking -lscai_tracing -lscai_logging -lscai_common

# Flags for googletest framework used for unit tests
IFLAGSTEST = -I${GTEST_DIR}/include
LDFLAGSTEST = -L${GTEST_DIR}/ -lgtest -lgtest_main
CFLAGSTEST = -isystem -pthread 

######################
# Install directory
######################
BINDIR=../bin

LIB_TARGET=$(BINDIR)/libSOFI.a
SOFI_TARGET=$(BINDIR)/SOFI
SOFI_UNITTEST_TARGET=$(BINDIR)/Tests/Test_unit
SOFI_INTEGRATIONTEST_TARGET=$(BINDIR)/Tests/Test_CompareSeismogram

################################################################
################################################################
# No input parameters are required below this line
################################################################
################################################################

######################
# General targets
######################
all: sofi tests

install: clean sofi 

tests: integrationTest unitTest

lib: $(LIB_TARGET)

sofi: $(SOFI_TARGET)

unitTest: $(SOFI_UNITTEST_TARGET)

integrationTest: $(SOFI_INTEGRATIONTEST_TARGET)

.PHONY : clean
clean: clean-bin clean-obj

######################
# Pattern rule
######################
.cpp.o:
	$(CC) $(CFLAGS) -c  $< -o $@ $(IFLAGS)

SUBDIRS := $(wildcard */)

######################
# Locate CPP files for compilation
######################
SOFI_LIB_CPP := $(wildcard $(addsuffix */*.cpp,$(SUBDIRS))) $(wildcard **/*.cpp)
SOFI_UNITTEST_CPP := $(wildcard $(addsuffix */*UnitTest.cpp,$(SUBDIRS))) $(wildcard */*UnitTest.cpp)
SOFI_INTEGRATIONTEST_CPP := $(wildcard **/*Test_CompareSeismogram.cpp)

######################
# Static SOFI library 
######################
SOFI_LIB_SCR := $(filter-out $(SOFI_UNITTEST_CPP) $(SOFI_INTEGRATIONTEST_CPP), $(SOFI_LIB_CPP))
SOFI_LIB_OBJ := $(SOFI_LIB_SCR:%.cpp=%.o)

$(LIB_TARGET): $(SOFI_LIB_OBJ)
	ar rcuvs $@ $^
	ranlib $@

#####################
# SOFI main program
#####################
SOFI_MAIN_CPP := SOFI.cpp
SOFI_MAIN_OBJ := $(SOFI_MAIN_CPP:%.cpp=%.o)

$(SOFI_TARGET): $(SOFI_MAIN_OBJ) $(LIB_TARGET)
	$(CC) $(CFLAGS) $^ -o $@ $(IFLAGS) $(LDFLAGS)

######################
# Unit tests
######################
SOFI_UNITTEST_OBJ = $(SOFI_UNITTEST_CPP:%.cpp=%.o)
$(SOFI_UNITTEST_OBJ): CFLAGS+=$(CFLAGSTEST)
$(SOFI_UNITTEST_OBJ): LDFLAGS +=$(LDFLAGSTEST)
$(SOFI_UNITTEST_OBJ): IFLAGS +=$(IFLAGSTEST)

$(SOFI_UNITTEST_TARGET): $(SOFI_UNITTEST_OBJ) $(LIB_TARGET)
	mkdir -p $(BINDIR)/Tests/
	$(CC) $(CFLAGS) $(CFLAGSTEST) $^ -o $@ $(IFLAGS) $(IFLAGSTEST) $(LDFLAGS) $(LDFLAGSTEST)

######################
# Integration tests
######################
SOFI_INTEGRATIONTEST_OBJ := $(SOFI_INTEGRATIONTEST_CPP:%.cpp=%.o)

$(SOFI_INTEGRATIONTEST_TARGET): $(SOFI_INTEGRATIONTEST_OBJ) $(LIB_TARGET)
	mkdir -p $(BINDIR)/Tests/
	$(CC) $(CFLAGS) $^ -o $@ $(IFLAGS) $(LDFLAGS)

######################
# Cleaning
######################
.PHONY : clean-bin
clean-bin:
	rm -rf $(BINDIR)/SOFI
	rm -rf $(BINDIR)/libSOFI.a
	rm -rf $(BINDIR)/Tests/Test_unit
	rm -rf $(BINDIR)/Tests/Test_CompareSeismogram

.PHONY : clean-obj
clean-obj:
	rm -rf $(SOFI_LIB_OBJ) $(SOFI_MAIN_OBJ) $(SOFI_UNITTEST_OBJ) $(SOFI_INTEGRATIONTEST_OBJ)
