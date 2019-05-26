# Makefile for
# @author Soheil Hajian

# The compiler
FC	= gfortran
#FC 	= f95
# Flags
FCFLAGS = # intialize
#FCFLAGS+= -O2	# optimize
FCFLAGS += -fbounds-check
FCFLAGS += -fstack-check
FCFLAGS += -fimplicit-none # checks the implicit-none
## double precision computation
FCFLAGS	+= -fdefault-real-8
## OpenMP for parallel computing
FCFLAGS += -fopenmp 
#
FCFLAGS += -Wall
FCFLAGS += -llapack 
#FCFLAGS += -pg 
#FCFLAGS += -lblas

PETSC_DIR = /home/soheil/Downloads/petsc-3.5.0

# check if we use PETSC or not
# $(PETSC) is a given argument in the terminal
ifeq ($(PETSC),1)
	CFLAGS = -I${PETSC_DIR}/include	
	FFLAGS = -I${PETSC_DIR}/include/finclude
	SOURCESC =
	SOURCESF = main_petsc.F90
	OBJ = $(SOURCESF:.F90=.o)
	CLEANFILES = ${OBJ} 

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

spmat.o: 
fd.o: geo2d.o spmat.o
main_petsc.o: fd.o
main: ${OBJ} fd.o geo2d.o spmat.o
	-${FLINKER} -o test $^ ${PETSC_SYS_LIB}	


#

else
# List of executables
PROGRAMS = main

all: $(PROGRAMS)

# Inclusions
geo2d.o: 
spmat.o:
fd.o: geo2d.o spmat.o
main.o: fd.o
main: fd.o geo2d.o spmat.o

# Building the PROGRAMS; I don't know for what '$@' stands.
# '$^' name of all prerequists 
%: %.o
	$(FC) $(FCFLAGS) -o $@ $^
# Building libraries and modules
%.o: %.f90
	$(FC) $(FCFLAGS) -c $<

clean:
	rm -f *.o *.mod

veryclean:
	rm -f *.o *.mod *.out *.mesh

backup:
	mkdir -p backup
	cp makefile *.f90 *.sh backup

endif

