import numpy as np
#This import statment initializes the MPI process
from mpi4py import MPI

#This subroutine obtains the nodes and weights for the 
#Gauss-Laguerre Quadrature.
def weights_abscissas( N ):
    return np.polynomial.laguerre.laggauss(N)

#triple_inte is the integral calculation routine. The integral is calculated 
#using Gauss-Laguerre Quadrature. 
def triple_inte( nn, N, x=None, w=None ):
    #Checking if weights and abscissas are calculated
    if x is None and w is None:
        x, w = weights_abscissas(N)
    #Calculation of triple integral
    sol = 0
    for i in range(N):
        for j in range(N):
            for k in range(N):
                sol += w[i]*w[j]*w[k]*np.exp(-nn*(x[i]**2+x[j]**2+x[k]**2))
    return sol

#Subroutine for finding out how many calculations each processor
#needs to compute. The calcs array holds this information.
#calcs[:,1] is the integral to be computed
#calcs[:,2] holds the processor number that will calculate this integral
def split_tasks( nproc, calc_num ):
    calcs = np.zeros( ( calc_num, 2 ), dtype=int )
    procnum = 1
    for i in range( calc_num ):
        calcs[ i, 0 ] = i + 1
        calcs[ i, 1 ] = procnum
        if procnum == nproc - 1:
            procnum = 1
        else:
            procnum += 1
    return calcs

#Subroutine for the task runners. They receive the tasks to be
#completed from the main processor.
def proc( iproc, icomm, calc_num, integrals, N=80 ):
    x, w = weights_abscissas( N )
    calcs = np.zeros( ( calc_num, 2 ), dtype=int )
    calcs = icomm.bcast( calcs, root = 0 )
    #Calculating the integrals
    for nn in range( calc_num ):
        if calcs[ nn, 1 ] == iproc:
            integrals[ nn ] = triple_inte( nn+1, N, x=x, w=w )
    return integrals

#Starting MPI
t1 = MPI.Wtime()                #Gets the initial time
icomm = MPI.COMM_WORLD
nproc = icomm.size              #Gets the number of processors
iproc = icomm.Get_rank()        #Gets the processor number
root_process = 0                #Sets the main processor as 0. All others are task runners.

#Setting up arrays and variables
calc_num = 1000
integrals = np.zeros( ( calc_num ), dtype=float )
res = np.zeros( ( calc_num ), dtype=float )

#Splitting program into main and task runners
#Main task splits the tasks for the task runners and send out the data.
#task runners receive data and calculate the integrals.
if iproc == root_process:
    calcs = split_tasks( nproc, calc_num )
    calcs = icomm.bcast( calcs, root = 0 )
else:
    integrals = proc( iproc, icomm, calc_num, integrals )

#Merge all the integral values into the res array. The integrals array holds zero
#unless the integral was calculated by the processor. To merge them we sum each
#array element from all processors. Since we're adding zero to the integral value,
#res will hold the integral value.
res = icomm.reduce( integrals, op=MPI.SUM, root=0 )

#Print integrals to check output
if iproc == root_process:
    #Used to check the integral values
    #for i in range( calc_num ):
    #    print( i+1, res[i] )
    t2 = MPI.Wtime()
    print('Code with', nproc, 'processors took', t2-t1, 'seconds.')