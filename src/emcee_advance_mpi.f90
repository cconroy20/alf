SUBROUTINE EMCEE_ADVANCE_MPI (ndim, nwalkers, a, pin, lpin, &
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
  
  USE mpi; USE alf_vars
  USE alf_utils, ONLY : func, myran, function_parallel_map
  IMPLICIT NONE

  INTEGER, INTENT(in)  :: ndim, nwalkers
  REAL(DP), INTENT(in) :: a
  REAL(DP), INTENT(in), DIMENSION(ndim,nwalkers) :: pin
  REAL(DP), INTENT(in), DIMENSION(nwalkers) :: lpin
  INTEGER, INTENT(in) :: nworkers
  REAL(DP), INTENT(inout), DIMENSION(ndim,nwalkers) :: pout
  REAL(DP), INTENT(inout), DIMENSION(nwalkers) :: lpout
  INTEGER, INTENT(out), DIMENSION(nwalkers) :: accept

  INTEGER  :: k, ri
  REAL(DP) :: z, diff
  REAL(DP), DIMENSION(nwalkers) :: zarr, lpnew
  REAL(DP), DIMENSION(ndim,nwalkers) :: qarr
                
  !---------------------------------------------------------------!
  !---------------------------------------------------------------!

  !Loop over the walkers to propose new positions
  DO k=1,nwalkers
     
     !Compute a random stretch factor and store it
     z = (a - 1.d0) * myran() + 1.d0
     z = z * z / a
     zarr(k) = z
     
     !Select the helper walker
     ri = CEILING((nwalkers-1) * myran())
     IF (ri.GE.k) THEN
        ri = ri + 2
     ENDIF
     
     !Compute the proposal position and store it
     qarr(:,k) = (1.d0 - z) * pin(:, ri) + z * pin(:, k)
     
  ENDDO
  
  CALL FUNCTION_PARALLEL_MAP(ndim, nwalkers, nworkers, qarr, lpnew)

  !Now loop over walkers to accept/reject, and update
  DO k=1,nwalkers
     
     diff = (ndim - 1.d0) * LOG(zarr(k)) + lpnew(k) - lpin(k)
     
     !Accept or reject
     IF (diff.GE.0.d0) THEN
        accept(k) = 1
     ELSE
        IF (diff.GE.LOG(myran())) THEN
           accept(k) = 1
        ELSE
           accept(k) = 0
        ENDIF
     ENDIF
     
     !Do the update
     IF (accept(k).EQ.1) then
        pout(:, k) = qarr(:, k)
        lpout(k)   = lpnew(k)
     ELSE
        pout(:, k) = pin(:, k)
        lpout(k)   = lpin(k)
     ENDIF
     
  ENDDO
  
END SUBROUTINE EMCEE_ADVANCE_MPI

