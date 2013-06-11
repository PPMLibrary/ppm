# Makefile to generate source files from the templates

.SUFFIXES: .def .f
.PHONY: clean

TYPES=ppm_t_neighlist_check.f ppm_c_neighlist_check.f ppm_t_container_check.f ppm_c_part_prop_check.f ppm_v_main_abstr_check.f ppm_t_mapping_check.f ppm_t_part_mapping_check.f

all: $(TYPES)

.def.f:
	./mkh5type $< $@

install: all
	cp -v *.f ../checkp

clean:
	rm *.f
