TESTINC=-I../../../include -I/usr/local/include -L../../../lib -lppm
sources=$(wildcard *.f)
targets=$(patsubst %.f,%,$(sources))
CPP = cpp
FC = mpif90

new : clean all

all : $(targets)

$(targets) : $(sources)
	rm -f __$@.f
	$(CPP) -P $@.f > __$@.f
	$(FC) $(OPTS) __$@.f $(TESTINC) $(LDFLAGS) $(LIBS) -o $@
	rm -f __$@.f
	./$@

clean :
	rm -f __*
	rm -f $(targets)
