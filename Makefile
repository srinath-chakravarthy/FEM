all:
	cd src; $(MAKE) all
clean:
	cd src; $(MAKE) clean
test:
	mpiexec.mpich -n 2 ./defmod -f examples/two_quads_qs.inp 
test-lin:
	mpiexec.mpich -n 2 ./defmod -f examples/lin_two_quads_qs.inp 
test-generated:
	mpiexec.mpich -n 2 ./defmod -f examples/generated_example.inp
