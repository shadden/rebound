PROBLEMDIR=$(shell basename `dirname \`pwd\``)"/"$(shell basename `pwd`)
export OPENGL=0
export OPT=-O3
export QUADRUPOLE=0
export MPI=0

all:
	# Setup link to different modules
	ln -fs gravity_none.c ../../src/gravity.c
	ln -fs boundaries_open.c ../../src/boundaries.c
	ln -fs ../$(PROBLEMDIR)/integrator_impulse2ndorder.c ../../src/integrator.c
	ln -fs collisions_none.c ../../src/collisions.c
	# Setup link to problem file
	ln -fs ../$(PROBLEMDIR)/problem.c ../../src/problem.c
	# Compile
	$(MAKE) -C ../../src/
	# Copy result
	cp ../../src/nbody .

leapfrog:
	# Setup link to different modules
	ln -fs gravity_direct.c ../../src/gravity.c
	ln -fs boundaries_open.c ../../src/boundaries.c
	ln -fs integrator_leapfrog.c ../../src/integrator.c
	ln -fs collisions_none.c ../../src/collisions.c
	# Setup link to problem file
	ln -fs ../$(PROBLEMDIR)/problem.c ../../src/problem.c
	# Compile
	$(MAKE) -C ../../src/
	# Copy result
	cp ../../src/nbody .

doc: all
	$(MAKE) -C ../../src/ doc

test:
	rm -v error.txt
	./nbody --dt=1
	./nbody --dt=0.2
	./nbody --dt=0.1
	./nbody --dt=0.02
	./nbody --dt=0.01
	./nbody --dt=0.002
	./nbody --dt=0.001
	./nbody --dt=0.0002
	./nbody --dt=0.0001

clean:
	$(MAKE) -C ../../src/ clean
	rm -vf nbody
