# Makefile for Defmod
export ARCH=linux_ib_intel12_openmpi
export PETSC_DIR=/opt/neu/petsc/${ARCH}

FFLAGS = -cpp -I${PETSC_DIR}/include -O3
FPPFLAGS = 

-include ${PETSC_DIR}/conf/variables
-include ${PETSC_DIR}/conf/rules

OBJS = m_seplaw.o m_global.o m_local.o m_elems.o m_io.o m_utils.o main.o
#.PHONY : print_vars
#print_vars : 
#	echo ${FLINKER}
m_utils.o : m_utils.f90
m_elems.o : m_elems.f90 m_utils.o
m_io.o    : m_io.f90 m_elems.o
m_local.o : m_local.f90 m_elems.o
m_seplaw.o: m_seplaw.f90 
m_global.o: m_global.f90 m_local.o m_seplaw.o
main.o    : main.f90 m_global.o

all: ${OBJS}
	-${FLINKER} ${FFLAGS} ${OBJS} -o ../defmod ${PETSC_LIB}
