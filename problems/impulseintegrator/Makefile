PROBLEMDIR=$(shell basename `dirname \`pwd\``)"/"$(shell basename `pwd`)
export OPENGL=1
export OPENMP=1
export OPT=-O3
export CC=gcc 
export QUADRUPOLE=0
export MPI=0

impulse:
	# Setup link to different modules
	ln -fs gravity_none.c ../../src/gravity.c
	ln -fs boundaries_open.c ../../src/boundaries.c
	ln -fs collisions_none.c ../../src/collisions.c
	ln -fs ../$(PROBLEMDIR)/integrator_impulse.c ../../src/integrator.c
	ln -fs ../$(PROBLEMDIR)/problem.c ../../src/problem.c
	# Compile
	$(MAKE) -C ../../src/
	# Copy result
	cp ../../src/nbody .

impulse2ndorder:
	# Setup link to different modules
	ln -fs gravity_none.c ../../src/gravity.c
	ln -fs boundaries_open.c ../../src/boundaries.c
	ln -fs collisions_none.c ../../src/collisions.c
	ln -fs ../$(PROBLEMDIR)/integrator_impulse2ndorder.c ../../src/integrator.c
	ln -fs ../$(PROBLEMDIR)/problem.c ../../src/problem.c
	# Compile
	$(MAKE) -C ../../src/
	# Copy result
	cp ../../src/nbody .

leapfrog:
	# Setup link to different modules
	#ln -fs gravity_tree.c ../../src/gravity.c
	ln -fs gravity_direct.c ../../src/gravity.c
	ln -fs boundaries_open.c ../../src/boundaries.c
	ln -fs collisions_none.c ../../src/collisions.c
	ln -fs integrator_leapfrog.c ../../src/integrator.c
	ln -fs ../$(PROBLEMDIR)/problem.c ../../src/problem.c
	# Compile
	$(MAKE) -C ../../src/
	# Copy result
	cp ../../src/nbody .

euler:
	# Setup link to different modules
	ln -fs gravity_direct.c ../../src/gravity.c
	ln -fs boundaries_open.c ../../src/boundaries.c
	ln -fs collisions_none.c ../../src/collisions.c
	ln -fs integrator_euler.c ../../src/integrator.c
	ln -fs ../$(PROBLEMDIR)/problem.c ../../src/problem.c
	# Compile
	$(MAKE) -C ../../src/
	# Copy result
	cp ../../src/nbody .

modifiedeuler:
	# Setup link to different modules
	ln -fs gravity_direct.c ../../src/gravity.c
	ln -fs boundaries_open.c ../../src/boundaries.c
	ln -fs collisions_none.c ../../src/collisions.c
	ln -fs ../$(PROBLEMDIR)/integrator_modifiedeuler.c ../../src/integrator.c
	ln -fs ../$(PROBLEMDIR)/problem.c ../../src/problem.c
	# Compile
	$(MAKE) -C ../../src/
	# Copy result
	cp ../../src/nbody .

radau15:
	# Setup link to different modules
	ln -fs gravity_direct.c ../../src/gravity.c
	ln -fs boundaries_open.c ../../src/boundaries.c
	ln -fs collisions_none.c ../../src/collisions.c
	ln -fs integrator_radau15.c ../../src/integrator.c
	ln -fs ../$(PROBLEMDIR)/problem.c ../../src/problem.c
	# Compile
	$(MAKE) -C ../../src/
	# Copy result
	cp ../../src/nbody .

doc: all
	$(MAKE) -C ../../src/ doc

testall:
	touch error_leapfrog.txt
	rm -v error_*.txt
	make radau15
	./nbody --dt=0.00001 --write=1	
	make leapfrog
	make test | grep Diff
	make euler
	make test | grep Diff
	make modifiedeuler
	make test | grep Diff
	make impulse
	make test | grep Diff
	make impulse2ndorder
	make test | grep Diff
	make radau15
	make test | grep Diff

test:
	./nbody --dt=1
	./nbody --dt=0.5
	./nbody --dt=0.25
	./nbody --dt=0.125
	./nbody --dt=0.1
	./nbody --dt=0.0625
	./nbody --dt=0.05
	./nbody --dt=0.03125
	./nbody --dt=0.025
	./nbody --dt=0.0125
	./nbody --dt=0.01
	./nbody --dt=0.00625
	./nbody --dt=0.005
	./nbody --dt=0.003125
	./nbody --dt=0.0025
	./nbody --dt=0.00125
	./nbody --dt=0.001
	./nbody --dt=0.000625
	./nbody --dt=0.0005
	./nbody --dt=0.0003125
	./nbody --dt=0.00025
	./nbody --dt=0.000125
	./nbody --dt=0.0001
	./nbody --dt=0.0000625
	./nbody --dt=0.00005
	./nbody --dt=0.00003125
	./nbody --dt=0.000025
	./nbody --dt=0.0000125
	./nbody --dt=0.00001

clean:
	$(MAKE) -C ../../src/ clean
	rm -vf nbody
