SUBROUTINE EMCEE_ADVANCE(ndim,nwalkers,a,pin,lpin,pout,lpout,accept)

  ! This routine was written by Dan Foreman-Mackey, based on his
  ! python code emcee.

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
  !---------------------------------------------------------------!

  USE alf_vars; USE alf_utils, ONLY : func, myran
  IMPLICIT NONE

  INTEGER, INTENT(in) :: ndim, nwalkers
  REAL(DP), INTENT(in) :: a
  REAL(DP), INTENT(in), DIMENSION(ndim,nwalkers) :: pin
  REAL(DP), INTENT(in), DIMENSION(nwalkers) :: lpin
  
  REAL(DP), INTENT(out), DIMENSION(ndim,nwalkers) :: pout
  REAL(DP), INTENT(out), DIMENSION(nwalkers) :: lpout
  INTEGER, INTENT(out),  DIMENSION(nwalkers) :: accept
  
  INTEGER  :: k, ri
  REAL(DP) :: z, lp, diff
  REAL(DP), DIMENSION(ndim) :: q
 
  !---------------------------------------------------------------!
  !---------------------------------------------------------------!

  DO k=1,nwalkers
           
     ! Compute a random stretch factor
     z = (a-1.d0) * myran() + 1.d0
     z = z*z/a
     
     ! Select the helper walker
     ri = CEILING((nwalkers-1) * myran())
     IF (ri.GE.k) THEN
        ri = ri+1
        q  = pin(:,ri+1)
     ELSE
        q = pout(:,ri)
     ENDIF
     
     ! Compute the proposal position
     q = (1.d0-z) * q + z*pin(:,k)

     ! Compute the new ln-probability
     lp   = -0.5*func(q) !note: func returns chi^2
     diff = (ndim-1.d0) * LOG(z) + lp - lpin(k)
     
     ! Accept or reject
     IF (diff.GE.0.d0) THEN
        accept(k) = 1
     ELSE
        IF (diff.GE.LOG(myran())) THEN
           accept(k) = 1
        ELSE
           accept(k) = 0
        ENDIF
     ENDIF
     
     ! Do the update
     IF (accept(k).EQ.1) THEN
        pout(:,k) = q
        lpout(k)  = lp
     ELSE
        pout(:,k) = pin(:,k)
        lpout(k)  = lpin(k)
     ENDIF
     
  ENDDO

  
END SUBROUTINE EMCEE_ADVANCE
