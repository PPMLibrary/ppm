
TESTINC=-I../../../include -L../../../lib -lppm
sources=$(wildcard *.f)
targets=$(patsubst %.f,%,$(sources))

all : $(targets)

$(targets) : $(sources)
	$(FC) $(OPTS) $< $(TESTINC) $(LDFLAGS) $(LIBS) -o $@ 
