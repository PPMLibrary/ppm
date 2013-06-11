# Makefile to generate source files from the templates

.SUFFIXES: .def .f
.PHONY: clean

TYPES= ppm_t_container_check.f ppm_t_main_abstr_check.f ppm_t_mapping_check.f ppm_t_neighlist_check.f ppm_t_operator_discr_check.f particles_stats_check.f

all: $(TYPES)

.def.f:
	./typefc.pl $< $@

#install: all
#	cp -v *.f ../checkp

clean:
	rm *.f
