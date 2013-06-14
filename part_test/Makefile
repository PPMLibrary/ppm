.PHONY: all

all:
	$(MAKE) -C gen

test: all
	LD_LIBRARY_PATH=/home/kevin/hdf5-1.8.11/hdf5/lib mpirun -np 1 bin/part_test
