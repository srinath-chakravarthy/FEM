# Makefile for Defmod

FFLAGS =
FPPFLAGS =

-include ${PETSC_DIR}/conf/variables
-include ${PETSC_DIR}/conf/rules

OBJS = m_seplaw.o m_global.o m_local.o m_elems.o m_io.o m_utils.o main.o  

m_utils.o : m_utils.F90
m_elems.o : m_elems.F90 m_utils.o
m_io.o    : m_io.F90 m_elems.o
m_local.o : m_local.F90 m_elems.o
m_seplaw.o: m_seplaw.F90
m_global.o: m_global.F90 m_local.o
main.o    : main.F90 m_global.o

all: ${OBJS}
	-${FLINKER} ${OBJS} -o ../defmod ${PETSC_LIB}