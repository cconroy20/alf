subroutine emcee_advance_mpi (ndim, nwalkers, a, pin, lpin, &
     pout, lpout, accept, nworkers)

  ! This subroutine advances an ensemble of walkers using the
  ! Goodman & Weare stretch move.
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
  ! a [double precision]:
  !   The proposal scale (a tuning parameter). Using `a=2` is almost
  !   always the right move.
  !
  ! pin [double precision (ndim, nwalkers)]:
  !   The starting positions of the walkers in parameter space.
  !
  ! lpin [double precision (nwalkers)]:
  !   The value of the log-probability function at positions `pin`.
  !
  ! nworkers [integer]:
  !   The number of available workers, not including the master
  !   process, which is assumed to be workerid=0
  
  ! Outputs
  ! -------
  !
  ! pout [double precision (ndim, nwalkers)]:
  !   The final positions of the walkers in parameter space.
  !
  ! lpout [double precision (nwalkers)]:
  !   The value of the log-probability function at positions `pout`.
  !
  ! accept [integer (nwalkers)]:
  !   A binary list indicating whether or not each proposal was
  !   accepted.
  
  USE mpi
  USE alf_vars; USE alf_utils, ONLY : func, myran, function_parallel_map
  IMPLICIT NONE

  integer, intent(in) :: ndim, nwalkers
  double precision, intent(in) :: a
  double precision, intent(in), dimension(ndim,nwalkers) :: pin
  double precision, intent(in), dimension(nwalkers) :: lpin

  double precision, intent(out), dimension(ndim,nwalkers) :: pout
  double precision, intent(out), dimension(nwalkers) :: lpout
  integer, intent(out), dimension(nwalkers) :: accept

  integer :: k, ri
  double precision :: r, z, diff
        
  INTEGER, intent(in) :: nworkers
  DOUBLE PRECISION, dimension(nwalkers) :: zarr, lpnew
  DOUBLE PRECISION, dimension(ndim,nwalkers) :: qarr
                
  !rqst = MPI_REQUEST_NULL

  ! Loop over the walkers to propose new positions
  do k=1,nwalkers
     
     ! Compute a random stretch factor and store it
     z = (a - 1.d0) * myran() + 1.d0
     z = z * z / a
     zarr(k) = z
     
     ! Select the helper walker.
     ri = ceiling((nwalkers-1) * myran())
     if (ri .ge. k) then
        ri = ri + 2
     endif
     
     ! Compute the proposal position and store it
     qarr(:,k) = (1.d0 - z) * pin(:, ri) + z * pin(:, k)
     
  enddo
  
  call function_parallel_map(ndim, nwalkers, nworkers, qarr, lpnew)
  
  ! Now loop over walkers to accept/reject, and update
  do k=1,nwalkers
     
     diff = (ndim - 1.d0) * log(zarr(k)) + lpnew(k) - lpin(k)
     
     ! Accept or reject.
     if (diff .ge. 0.d0) then
        accept(k) = 1
     else
        if (diff .ge. log(myran())) then
           accept(k) = 1
        else
           accept(k) = 0
        endif
     endif
     
     ! Do the update.
     if (accept(k) .eq. 1) then
        pout(:, k) = qarr(:, k)
        lpout(k) = lpnew(k)
     else
        pout(:, k) = pin(:, k)
        lpout(k) = lpin(k)
     endif
     
  enddo
  
end subroutine emcee_advance_mpi
      
