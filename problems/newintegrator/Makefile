PROBLEMDIR=problems/$(shell basename `pwd`)
export OPT=-O3
export OPENGL=1
#export QUADRUPOLE=1
#export MPI=1
#export OPENMP=1
#export LIBPNG=1

all:
	# Setup link to different modules
	ln -fs ../$(PROBLEMDIR)/gravity_direct_testparticle.c ../../src/gravity.c
	#ln -fs gravity_direct.c ../../src/gravity.c
	ln -fs boundaries_none.c ../../src/boundaries.c
	ln -fs ../$(PROBLEMDIR)/integrator_fivestep.c ../../src/integrator.c
	#ln -fs integrator_leapfrog.c ../../src/integrator.c
	ln -fs collisions_none.c ../../src/collisions.c
	# Setup link to problem file
	ln -fs ../$(PROBLEMDIR)/problem.c ../../src/problem.c
	# Compile
	$(MAKE) -C ../../src/
	# Copy result
	cp ../../src/nbody .

doc: all
	$(MAKE) -C ../../src/ doc

clean:
	$(MAKE) -C ../../src/ clean
	rm -vf nbody

test: all
	./nbody 2e-4
	./nbody 1e-4
	./nbody 2e-5
	./nbody 1e-5
	./nbody 2e-6
	./nbody 1e-6
	./nbody 2e-7
	./nbody 1e-7
