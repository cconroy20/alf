subroutine free_workers(nworkers)
        !
        ! This subroutine sends dummy position arrays to each slave,
        ! and uses a tag value that is larger then the total number of
        ! walkers.  The slaves should interpret this tag as a signal
        ! to break out of their event loops.
        !
        ! Inputs
        ! ------
        !
        ! ndim [integer]:
        !   The dimension of the parameter space.
        !
        ! nwalkers [integer]:
        !   The number of walkers.
        !
        ! nworkers [integer]:
        !   The number of worker processes that need to be closed
        !
        
  use mpi
  implicit none
  
  integer, intent(in) :: nworkers
  integer :: k, ierr, FREE=99, dummy=0
  
  do k=1,nworkers
     call MPI_SEND(dummy, 1, MPI_INTEGER, &
          k, FREE, MPI_COMM_WORLD, ierr)
  enddo
  
end subroutine free_workers
