PROGRAM MPI_Tutorial
    IMPLICIT NONE
	INCLUDE "mpif.h"
    INTEGER :: icomm, ierr, nproc, iproc, root_process
    INTEGER :: calc_num, i
    REAL*8, DIMENSION(:), ALLOCATABLE :: integrals, res
    REAL :: t1, t2

    !Initialize the MPI process
    CALL MPI_INIT( ierr )
    !Initial time check
    t1 = MPI_Wtime()
    icomm = MPI_COMM_WORLD
    CALL MPI_COMM_SIZE( icomm, nproc, ierr )
    CALL MPI_COMM_RANK( icomm, iproc, ierr )

    !root_process is the variable that determines the main processor.
    !all other processors are task runners.
    root_process = 0

    !Setting up arrays
    calc_num = 1000
    ALLOCATE(integrals(calc_num))
    ALLOCATE(res(calc_num))
    integrals(:) = 0.D0
    res(:) = 0.D0

    !This IF statement splits the code by the main processor
    !and the task runners.
    IF ( iproc == root_process ) THEN
        IF ( nproc == 1 ) THEN
            CALL single( calc_num, integrals )
        ELSE
            CALL mproc( nproc, icomm, ierr, calc_num )
        END IF
    END IF
    CALL proc( iproc, icomm, ierr, calc_num, integrals )

    !Wait here for all to finish before ending MPI process
    CALL MPI_BARRIER( icomm, ierr )

    !Merge all the integral values into the res array. The integrals array holds zero
    !unless the integral was calculated by the processor. To merge them we sum each
    !array element from all processors. Since we're adding zero to the integral value,
    !res will hold the integral value.
    CALL MPI_REDUCE( integrals, res, calc_num*2, MPI_REAL, MPI_SUM, 0, icomm, ierr )

    !This do loop prints out the integral values. Used to make sure code runs correctly.
    DO i=1,calc_num
        IF ( iproc == root_process ) print *, i, res(i)
    END DO

    !End time check
    t2 = MPI_Wtime()
    IF ( iproc == root_process ) &
        & print *, 'Code with', nproc, 'processors took', t2-t1, 'seconds.'
    CALL MPI_FINALIZE( ierr )
END PROGRAM

SUBROUTINE single( calc_num, integrals )
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: calc_num
    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: calcs
    REAL*8, ALLOCATABLE, DIMENSION(:) :: x, w
    REAL*8, DIMENSION(calc_num), INTENT(INOUT) :: integrals
    INTEGER :: i, n
    ALLOCATE(calcs(calc_num,2))
    DO i=1,calc_num
        calcs(i,1) = i
        calcs(i,2) = 0
    END DO

    !n is the number of quadrature points to use
    !Change here if needed as well as in the proc subroutine
    n = 150
    ALLOCATE(x(n))
    ALLOCATE(w(n))

    !These two routines calculate the integral.
    CALL gaulag( x, w, n )
    CALL calc_inte( calcs, 0, x, w, n, calc_num, integrals )
END SUBROUTINE

!Subroutine for the main processor. This is the one that sends out
!tasks for the other processor to complete.
SUBROUTINE mproc( nproc, icomm, ierr, calc_num )
    IMPLICIT NONE
	INCLUDE "mpif.h"
    INTEGER, INTENT(IN) :: nproc, icomm, ierr, calc_num
    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: calcs
    INTEGER :: i
    ALLOCATE(calcs(calc_num,2))

    !split_tasks tells the each task runner what integrals they will calculate.
    CALL split_tasks( nproc, calc_num, calcs )

    !This MPI_BCAST call sends the calcs array to each processor.
    CALL MPI_BCAST( calcs, 2*calc_num, MPI_INTEGER, 0, icomm, ierr )
END SUBROUTINE mproc

!Subroutine for the task runners. They receive the tasks to be
!completed from the main processor.
SUBROUTINE proc( iproc, icomm, ierr, calc_num, integrals )
    IMPLICIT NONE
	INCLUDE "mpif.h"
    INTEGER, INTENT(IN) :: iproc, icomm, ierr, calc_num
    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: calcs
    REAL*8, ALLOCATABLE, DIMENSION(:) :: x, w
    REAL*8, DIMENSION(calc_num), INTENT(INOUT) :: integrals
    INTEGER :: i, n

    !n is the number of quadrature points to use
    !Change here if needed as well as in the single subroutine
    n = 150
    ALLOCATE(calcs(calc_num,2))
    ALLOCATE(x(n))
    ALLOCATE(w(n))
    
    !This MPI_BCAST call accepts the calcs array from the root processor.
    IF ( .NOT. iproc == 0 )CALL MPI_BCAST( calcs, 2*calc_num, MPI_INTEGER, 0, icomm, ierr )

    !These two routines calculate the integral.
    CALL gaulag( x, w, n )
    CALL calc_inte( calcs, iproc, x, w, n, calc_num, integrals )
END SUBROUTINE proc

!Subroutine for finding out how many calculations each processor
!needs to compute. The calcs array holds this information.
!calcs(:,1) is the integral to be computed
!calcs(:,2) holds the processor number that will calculate this integral
SUBROUTINE split_tasks( nproc, calc_num, calcs )
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nproc, calc_num
    INTEGER, DIMENSION(calc_num,2), INTENT(INOUT) :: calcs
    INTEGER :: i, procnum
    !Setting the processor number for each calculation
    procnum = 0
    DO i=1,calc_num
        calcs(i,1) = i
        calcs(i,2) = procnum
        IF ( procnum == nproc - 1 ) THEN
            procnum = 0
        ELSE
            procnum = procnum + 1
        END IF
    END DO
    !This prints each integral to be calculated and what processor will
    !calculate it. Uncomment for debugging.
    !DO i=1,calc_num
    !    print *, calcs(i,1), calcs(i,2)
    !END DO
END SUBROUTINE split_tasks

!calc_inte is the subroutine that calculates the integrals that each processor
!is supposed to based on the calcs array.
SUBROUTINE calc_inte( calcs, iproc, x, w, n, calc_num, integrals )
    IMPLICIT NONE
    INTEGER, DIMENSION(calc_num,2), INTENT(IN) :: calcs
    INTEGER, INTENT(IN) :: iproc, n, calc_num
    REAL*8, DIMENSION(n), INTENT(IN) :: x, w
    REAL*8, DIMENSION(calc_num), INTENT(INOUT) :: integrals
    REAL*8 :: integra
    INTEGER :: i

    DO i=1,calc_num
        IF ( calcs( i, 2 ) == iproc) integrals(i) = integra( x, w, n, i )
    END DO
END SUBROUTINE

!integra is the integral calculation routine. The integral is calculated 
!using Gauss-Laguerre Quadrature. 
FUNCTION integra( x, w, n, nn )
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n, nn
    REAL*8, DIMENSION(n), INTENT(IN) :: x, w
    INTEGER :: i, j, k
    REAL*8 :: integra, inte_func
    integra = 0.D0
    
    DO i=1,n
        DO j=1,n
            DO k=1,n
                integra = integra + w(i)*w(j)*w(k)*inte_func( x(i), x(j), x(k), nn )
            END DO
        END DO
    END DO
    RETURN
END FUNCTION

!This is the function that is being integrated without e^(-x-y-z)
!according to the Gauss-Laguerre Quadrature.
FUNCTION inte_func( x, y, z, n )
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: x, y, z
    INTEGER, INTENT(IN) :: n
    REAL*8 :: inte_func

    inte_func = EXP(-n*(x**2.D0+z**2.D0+y**2.D0))
END FUNCTION

!This subroutine calculates the nodes and weights for the 
!Gauss-Laguerre Quadrature. Obtained from Numerical Recipes in FORTRAN 90.
SUBROUTINE gaulag(x,w,n)
    IMPLICIT NONE
	INTEGER :: nfi,n
	REAL*8, DIMENSION(n) :: x, w
	REAL*8, DIMENSION(n) :: anu
	REAL*8, PARAMETER :: eps=3.0d-13
	REAL*8, PARAMETER :: pi=4.d0*atan(1.d0)
	INTEGER :: j,its,i
	INTEGER, PARAMETER :: maxit=10
	REAL*8, PARAMETER :: c1=9.084064d-1,c2=5.214976d-2
	REAL*8, PARAMETER :: c3=2.579930d-3,c4=3.986126d-3
	REAL*8, DIMENSION(n) :: rhs,r2,r3,theta,p1,p2,p3
	REAL*8, DIMENSION(n) :: pp(n),z(n),z1(n)
    LOGICAL :: unfinished(n)
	x=0.d0
	w=0.d0
	anu = 4.d0*n+2.d0
	rhs(1) = (4.d0*n-1.d0)*pi/anu(1)
	do i=2,n
		rhs(i)=(4.d0*n-4.d0*i+3.d0)*pi/anu(i)
	end do
	r3 = rhs**(1.d0/3.d0)
	r2 = r3**2.d0
	theta = r3*(c1+r2*(c2+r2*(c3+r2*c4)))
	z = anu*(Cos(theta)**2.d0)
	unfinished = .true.
	do its=1,maxit
		where (unfinished)
			p1=1.d0
			p2=0.d0
		end where
		do j=1,n
			where (unfinished)
				p3=p2
				p2=p1
				p1=((2.d0*real(j,8)-1.d0-z)*p2-(real(j,8)-1.d0)*p3)/real(j,8)
			end where
		end do
		where (unfinished)
			pp=(n*p1-n*p2)/z
			z1=z
			z=z1-p1/pp
			unfinished = (abs(z-z1) > eps*z)
		end where
		if (.not. any(unfinished)) exit
	end do
	if (its == maxit+1) print *, 'too many iterations in gaulag'
	x=z
	w=-1.d0/(pp*real(n,8)*p2)
End Subroutine