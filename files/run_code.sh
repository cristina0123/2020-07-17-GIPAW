#!/bin/bash
#Setting number of processors to be used.
nproc=4
#Setting calculation type. 
#calc="Python"
calc="Fortran"
if [ $calc == "Fortran" ] ; then
	echo "Running Fortran test code."
	mpif90 MPI_Tutorial.f90
	mpirun -np $nproc ./a.out
elif [ $calc == "Python" ] ; then
	echo "Running Python test code."
	mpirun -np $nproc python MPI_Tutorial.py
fi
