N ?= 100										# Make default variable
Np ?= 2

program: main.c
	mpic++ -std=c++11 -fopenmp -o program main.c

test: program							# Have program as prerequisit
	mpirun -np $(Np) ./program $(N)

utest: program
	./program -1

vtest: program
	./program -2

convtest: program
	mpirun -np 2 ./program -3 1;
	mpirun -np 4 ./program -3 1;
	mpirun -np 8 ./program -3 1;
	mpirun -np 16 ./program -3 1;
	mpirun -np 32 ./program -3 1;
	mpirun -np 64 ./program -3 1;

	mpirun -np 2 ./program -3 2;
	mpirun -np 4 ./program -3 2;
	mpirun -np 8 ./program -3 2;
	mpirun -np 16 ./program -3 2;
	mpirun -np 32 ./program -3 2;
	mpirun -np 64 ./program -3 2;

	mpirun -np 2 ./program -3 4;
	mpirun -np 4 ./program -3 4;
	mpirun -np 8 ./program -3 4;
	mpirun -np 16 ./program -3 4;
	mpirun -np 32 ./program -3 4;
	mpirun -np 64 ./program -3 4;

	mpirun -np 2 ./program -3 8;
	mpirun -np 4 ./program -3 8;
	mpirun -np 8 ./program -3 8;
	mpirun -np 16 ./program -3 8;
	mpirun -np 32 ./program -3 8;
	mpirun -np 64 ./program -3 8;
