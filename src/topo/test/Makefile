#----------------------------------------------------------------------
# PATHS (need to be adjusted)
#----------------------------------------------------------------------
# PPM Library
PPMDIR     = ../../..
# METIS Library
METISDIR   = ../../../../metis4
# MPI Library
MPIDIR     = /usr/local

CPP = cpp
FC = mpif90 --free-form -fbounds-check
#----------------------------------------------------------------------

TESTINC=-I$(PPMDIR)/include -I$(MPIDIR)/include -L$(PPMDIR)/lib -lppm -L$(METISDIR)/lib -lmetis
sources=$(wildcard *.f)
targets=$(patsubst %.f,%,$(sources))


new : clean all

all : $(targets)

$(targets) : $(sources)
	$(CPP) -P $@.f > __$@.f
	$(FC) $(OPTS) __$@.f $(TESTINC) $(LDFLAGS) $(LIBS) -o $@
	rm -f __$@.f

clean :
	rm -f __*
	rm -f $(targets)
	rm -f ppm_dbg*
	rm -f ppm_log*
	
