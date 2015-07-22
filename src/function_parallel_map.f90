      subroutine function_parallel_map(ndim, nk, nworkers, pos, lnpout)
        !
        ! This subroutine sends rows of the pos array to whatever
        ! function is set up to receive them in a different process,
        ! and collects the results.  This could probably be done with
        ! scatter/gather, but I don't know how to use those yet.
        !
        ! Inputs
        ! ------
        !
        ! ndim [integer]:
        !   The dimension of the parameter space.
        !
        ! nk [integer]:
        !   The number of walkers.
        !        
        ! nworkers [integer]:
        !   The number of worker processes.
        !
        ! pos [double precision (ndim, nk)]:
        !   The positions of the walkers in parameter space.
        !
        ! Outputs
        ! ------
        !
        ! lnpout [double precision (nwalkers)]:
        !   The value of the log-probability function at positions `pos`.
        !

        use mpi
        implicit none

        integer, intent(in) :: ndim, nk, nworkers
        double precision, intent(in), dimension(ndim,nk) :: pos
        double precision, intent(out), dimension(nk) :: lnpout

        integer :: k, ierr, status(MPI_STATUS_SIZE), BEGIN=0
        integer :: npos, walk_per_work, extra, offset
        !integer, dimension(nk) :: rqst, rstat

        !rqst = MPI_REQUEST_NULL
        !Compute useful numbers for worker to walker ratio
        walk_per_work = nk/nworkers
        ! number of extra jobs or positions that need to be spread
        ! amongs the first set of workers
        extra = mod(nk,nworkers)
        
        ! Send chunks of new positions to workers
        offset = 1
        do k=1,nworkers
           if (k .le. extra) then
              !add an extra position
              npos = walk_per_work + 1
           else
              npos = walk_per_work
           endif
           ! Tell the worker how many positions to expect
           call MPI_SEND(npos, 1, MPI_INTEGER, &
                k, BEGIN, MPI_COMM_WORLD, ierr)
           ! Dispatch proposals to worker to figure out lnp
           call MPI_SEND(pos(1,offset), ndim*npos, MPI_DOUBLE_PRECISION, &
                k, BEGIN, MPI_COMM_WORLD, ierr)
           
           ! now increment offset
           offset = offset + npos
        enddo

        ! Loop over the workers to get the proposal lnp
        offset=1
        do k=1,nworkers
           if (k .le. extra) then
              !add an extra position
              npos = walk_per_work + 1
           else
              npos = walk_per_work
           endif

           ! get the lnps from the workers and store
           call MPI_RECV(lnpout(offset), npos, MPI_DOUBLE_PRECISION, &
                k, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
           
           offset = offset + npos
        enddo
        
      end subroutine function_parallel_map
      
