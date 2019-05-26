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


# List of executables
PROGRAMS = scratch

all: $(PROGRAMS)

# Inclusions
geo2d.o: 
spmat.o:
fd.o: geo2d.o spmat.o
coarsetools.o: geo2d.o spmat.o
scratch.o: fd.o
scratch: fd.o geo2d.o spmat.o

main.o: fd.o
main: fd.o geo2d.o spmat.o

test1.o: fd.o coarsetools.o
test1: fd.o geo2d.o spmat.o coarsetools.o

# Building the PROGRAMS; I don't know for what '$@' stands.
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
