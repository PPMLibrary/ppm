# Makefile to generate source files from the templates

.SUFFIXES: .def .f
.PHONY: clean all

TYPES= ppm_t_container_check.f ppm_t_main_abstr_check.f ppm_t_mapping_check.f ppm_t_neighlist_check.f ppm_t_operator_discr_check.f particles_stats_check.f

all: $(TYPES) type.f.in

#.def.f:
%.f:%.def type.f.in typefc.pl
	./typefc.pl $< $@

remake: clean all
#install: all
#	cp -v *.f ../checkp

clean:
	rm -vf *.f
