#!/bin/bash
CSHELL=sh

# --- MACROS
# define program name
EXE= infoli.x

# define C source files
SRC= infoli.c
UTILS= utils.c
OMPSS_SRC=ompss/infoli.c ompss/simulation.c ompss/infoli_log.c ompss/infoli_perf.c

# define C header files
HDR= infoli.h
OMPSS_HEADERS= ompss/infoli.h ompss/infoli_log.h ompss/infoli_perf.h

# define C object files
OBJ= infoli.o
UTILS_O= utils.o
UTILS_PHI_O= utils_phi.o
OMPSS_OBJ=${OMPSS_SRC:.c=.o}

# common directives
COMP_FLAGS= -c -Wall -O3 -Iutils -xCORE-AVX512 -qopt-report-phase=vec -qopt-report-phase=ipo -opt-subscript-in-range -lm
COMPSS_FLAGS= -std=gnu99 -c -Wall -O3 -g -ggdb -Iutils -xCORE-AVX512 -qopt-report-phase=vec -qopt-report=5 -qopt-subscript-in-range -qopt-report-phase=ipo -lm -lrt
UTILS_FLAGS= -c -Wall -O3 -lm
LINK_FLAGS= -O3 -lm

# regression check directive
REG_TEST_FLAGS= -DG_CAL_FROM_FILE

# profiling directive
PROFILING= -DPROFILING

# mic directives
MIC= -mmic

# define mic power measurement files
MICPWRDIR=../lib/
MICPWRH= ${MICPWRDIR}/power/mic_power.h
MICPWRSO= mic_power
IMFSO= imf
INTLCSO= intlc
IRNGSO= irng
SVMLSO= svml

# version directives
OMP_FLAGS= -openmp
HYBRID_FLAGS= -mt_mpi -openmp
OMPSS_FLAGS= -qopenmp

# subdirectories
OMP_DIR= openmp/
OMPSS_DIR= ompss/
HYBRID_DIR= hybrid/
UTILS_DIR= utils/
REG_DIR= regression_test/

# bin directory
BIN_DIR= ../bin/

# --- TARGETS
all: hybrid_xeon hybrid_phi omp_xeon omp_phi ompss_xeon


omp_xeon: ${OMP_DIR}${SRC} ${OMP_DIR}${HDR} ${UTILS_DIR}${UTILS_O}
	@echo #
	@echo "Building For Pure OpenMP Version on Xeon"
	icc ${OMP_DIR}${SRC} -o ${OMP_DIR}${OBJ} ${OMP_FLAGS} ${COMP_FLAGS}
	icc ${OMP_DIR}${OBJ} ${UTILS_DIR}${UTILS_O} -o ${BIN_DIR}${EXE} ${OMP_FLAGS} ${LINK_FLAGS}

ompss_xeon: ${OMPSS_OBJ} ${OMPSS_HDR} ${UTILS_DIR}${UTILS_O}
	@echo #
	@echo "Building For Pure OpenMP Version on Xeon"
	icc ${OMPSS_OBJ} ${UTILS_DIR}${UTILS_O} -o ${BIN_DIR}${EXE} ${OMPSS_FLAGS} ${LINK_FLAGS}

%.o: %.c
	icc -o $@ $< ${OMPSS_FLAGS} ${COMPSS_FLAGS}

omp_phi: ${OMP_DIR}${SRC} ${OMP_DIR}${HDR} ${UTILS_DIR}${UTILS_PHI_O}
	@echo #
	@echo "Building For Pure OpenMP Version on Phi"
	icc ${OMP_DIR}${SRC} -o ${OMP_DIR}${OBJ} ${MIC} ${OMP_FLAGS} ${COMP_FLAGS} 
	icc ${OMP_DIR}${OBJ} ${UTILS_DIR}${UTILS_PHI_O} -o ${BIN_DIR}${EXE} ${MIC} ${OMP_FLAGS} ${LINK_FLAGS} 

omp_phi_pow: ${OMP_DIR}${SRC} ${OMP_DIR}${HDR} ${UTILS_DIR}${UTILS_PHI_O}
	@echo #
	@echo "Building For Pure OpenMP Version on Phi"
	icc ${OMP_DIR}${SRC} -o ${OMP_DIR}${OBJ} ${MIC} ${OMP_FLAGS} ${COMP_FLAGS} -I${MICPWRDIR}
	icc ${OMP_DIR}${OBJ} ${UTILS_DIR}${UTILS_PHI_O} -o ${BIN_DIR}${EXE} ${MIC} ${OMP_FLAGS} ${LINK_FLAGS} -L${MICPWRDIR} \
	-l${MICPWRSO} -l${IMFSO} -l${INTLCSO} -l${IRNGSO} -l${SVMLSO}

omp_test: ${OMP_DIR}${SRC} ${OMP_DIR}${HDR} ${UTILS_DIR}${UTILS_O}
	@echo #
	@echo "Regression Test Build For OpenMP - will test on x86"
	icc ${OMP_DIR}${SRC} -o ${OMP_DIR}${OBJ} ${OMP_FLAGS} ${REG_TEST_FLAGS} ${COMP_FLAGS}
	icc ${OMP_DIR}${OBJ} ${UTILS_DIR}${UTILS_O} -o ${BIN_DIR}${EXE} ${OMP_FLAGS} ${REG_TEST_FLAGS} ${LINK_FLAGS}

hybrid_xeon: ${HYBRID_DIR}${SRC} ${HYBRID_DIR}${HDR} ${UTILS_DIR}${UTILS_O}
	@echo #
	@echo "Building For Hybrid Version on Xeon"
	mpiicc ${HYBRID_DIR}${SRC} -o ${HYBRID_DIR}${OBJ} ${HYBRID_FLAGS} ${COMP_FLAGS}
	mpiicc ${HYBRID_DIR}${OBJ} ${UTILS_DIR}${UTILS_O} -o ${BIN_DIR}${EXE} ${HYBRID_FLAGS} ${LINK_FLAGS}

hybrid_xeon_prof: ${HYBRID_DIR}${SRC} ${HYBRID_DIR}${HDR} ${UTILS_DIR}${UTILS_O}
	@echo #
	@echo "Building For Hybrid Version on Xeon with Profiling Info Enabled"
	mpiicc ${HYBRID_DIR}${SRC} -o ${HYBRID_DIR}${OBJ} ${HYBRID_FLAGS} ${PROFILING} ${COMP_FLAGS}
	mpiicc ${HYBRID_DIR}${OBJ} ${UTILS_DIR}${UTILS_O} -o ${BIN_DIR}${EXE} ${HYBRID_FLAGS} ${PROFLING} ${LINK_FLAGS}

hybrid_phi: ${HYBRID_DIR}${SRC} ${HYBRID_DIR}${HDR} ${UTILS_DIR}${UTILS_PHI_O}
	@echo #
	@echo "Building For Hybrid Version on Phi"
	mpiicc ${HYBRID_DIR}${SRC} -o ${HYBRID_DIR}${OBJ} ${MIC} ${HYBRID_FLAGS} ${COMP_FLAGS}
	mpiicc ${HYBRID_DIR}${OBJ} ${UTILS_DIR}${UTILS_PHI_O} -o ${BIN_DIR}${EXE} ${MIC} ${HYBRID_FLAGS} ${LINK_FLAGS}

hybrid_test: ${HYBRID_DIR}${SRC} ${HYBRID_DIR}${HDR} ${UTILS_DIR}${UTILS_O}
	@echo #
	@echo "Regression Test Build For Hybrid - will test on x86"
	mpiicc ${HYBRID_DIR}${SRC} -o ${HYBRID_DIR}${OBJ} ${HYBRID_FLAGS} ${REG_TEST_FLAGS} ${COMP_FLAGS}
	mpiicc ${HYBRID_DIR}${OBJ} ${UTILS_DIR}${UTILS_O} -o ${BIN_DIR}${EXE} ${HYBRID_FLAGS} ${REG_TEST_FLAGS} ${LINK_FLAGS}

${UTILS_DIR}${UTILS_O}: ${UTILS_DIR}${UTILS}
	@echo #
	@echo "Building Utilities"
	icc ${UTILS_DIR}${UTILS} -o ${UTILS_DIR}${UTILS_O} ${UTILS_FLAGS}

${UTILS_DIR}${UTILS_PHI_O}: ${UTILS_DIR}${UTILS}
	@echo #
	@echo "Building Utilities for Phi"
	icc ${UTILS_DIR}${UTILS} -o ${UTILS_DIR}${UTILS_PHI_O} ${MIC} ${UTILS_FLAGS}

clean:
	@echo #
	@echo "-- Cleaning Src --"
	rm -f ${OMP_DIR}*.o ${HYBRID_DIR}*.o ${BIN_DIR}*.x ${OMPSS_DIR}/*.o
	rm -f ${OMP_DIR}*.optrpt ${HYBRID_DIR}*.optrpt ${OMPSS_DIR}/*.optrpt
	rm -f ${UTILS_DIR}*.o
