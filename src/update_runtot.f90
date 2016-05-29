SUBROUTINE UPDATE_RUNTOT(runtot,inarr,m2l,m2lmw)

  !routine to update the array that holds the running totals
  !of the parameters, parameters^2, etc. in order to compute
  !averages and errors.

  USE alf_vars
  IMPLICIT NONE

  REAL(DP), INTENT(inout), DIMENSION(3,npar+2*nfil) :: runtot
  REAL(DP), INTENT(in), DIMENSION(nfil) :: m2l,m2lmw
  REAL(DP), INTENT(in), DIMENSION(npar) :: inarr

  !------------------------------------------------------------!

  runtot(1,:)      = runtot(1,:)+1.

  runtot(2,1:npar) = runtot(2,1:npar) + inarr
  runtot(3,1:npar) = runtot(3,1:npar) + inarr**2

  runtot(2,npar+1:npar+nfil) = runtot(2,npar+1:npar+nfil)+m2l
  runtot(3,npar+1:npar+nfil) = runtot(3,npar+1:npar+nfil)+m2l**2

  runtot(2,npar+nfil+1:npar+2*nfil) = &
       runtot(2,npar+nfil+1:npar+2*nfil)+m2lmw
  runtot(3,npar+nfil+1:npar+2*nfil) = &
       runtot(3,npar+nfil+1:npar+2*nfil)+m2lmw**2
  

END SUBROUTINE UPDATE_RUNTOT
