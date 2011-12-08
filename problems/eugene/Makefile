PROBLEMDIR=problems/$(shell basename `pwd`)
export OPENGL=0
export OPT=-O3
export QUADRUPOLE=0
export OPENMP=0
export MPI=1
export CC=mpicc

all:
	# Setup link to different modules
	ln -fs ../$(PROBLEMDIR)/gravity_direct_direct.c ../../src/gravity.c
	ln -fs ../$(PROBLEMDIR)/boundaries_shear_inflow.c ../../src/boundaries.c
	ln -fs integrator_sei.c ../../src/integrator.c
	ln -fs collisions_tree.c ../../src/collisions.c
	# Setup link to problem file
	ln -fs ../$(PROBLEMDIR)/problem.c ../../src/problem.c
	# Compile
	$(MAKE) -C ../../src/
	# Copy result
	cp ../../src/nbody .

doc: all
	$(MAKE) -C ../../src/ doc

clean:
	rm -v 25/m*
	rm -v 50/m*
	rm -v 100/m*
	rm -v 200/m*
	rm -v 25/posi*
	rm -v 50/posi*
	rm -v 100/posi*
	rm -v 200/posi*
