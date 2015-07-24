SUBROUTINE FREE_WORKERS(nworkers)

  ! This subroutine sends dummy position arrays to each slave,
  ! and uses a tag value that is larger then the total number of
  ! walkers.  The slaves should interpret this tag as a signal
  ! to break out of their event loops.
  !
  ! Inputs
  ! ------
  !
  ! ndim [integer]:
  !   The dimension of the parameter space
  !
  ! nwalkers [integer]:
  !   The number of walkers
  !
  ! nworkers [integer]:
  !   The number of worker processes that need to be closed
  !
        
  USE MPI
  IMPLICIT NONE
  
  INTEGER, INTENT(in) :: nworkers
  INTEGER :: k, ierr, FREE=99, dummy=0
  
  DO k=1,nworkers
     CALL MPI_SEND(dummy, 1, MPI_INTEGER, &
          k, FREE, MPI_COMM_WORLD, ierr)
  ENDDO
  
END SUBROUTINE FREE_WORKERS
