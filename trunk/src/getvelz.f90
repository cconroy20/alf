FUNCTION GETVELZ()

  !function to estimate the recession velocity
  !this routine is used to get a first-guess at the velocity
  !so that the subsequent Powell minimization (in specfit)
  !coverges faster.  uses 4100<lambda<5000A

  USE sfvars; USE nr, ONLY : locate
  USE sfutils, ONLY : linterp,contnormspec
  IMPLICIT NONE

  REAL(DP) :: getvelz, chi2,lo,hi
  INTEGER, PARAMETER :: nv=1000
  INTEGER :: i,i1,i2
  REAL(DP), DIMENSION(nl) :: mflx,dflx
  REAL(DP), DIMENSION(nv) :: tvz,tchi2,tvza
  TYPE(TDATA), DIMENSION(nl) :: idata
  REAL(DP), DIMENSION(ndat) :: tlam

  !------------------------------------------------------!
 
  chi2 = huge_number
  getvelz = 0.0

  DO i=1,nv

     tvz(i) = REAL(i)*11-1E3

     !de-redshift the data and interpolate to model wave array
     tlam      = data%lam / (1+tvz(i)/clight*1E5)
     idata(1:nl_fit)%flx = linterp(tlam(1:datmax),&
          data(1:datmax)%flx,sspgrid%lam(1:nl_fit))
     idata(1:nl_fit)%err = linterp(tlam(1:datmax),&
          data(1:datmax)%err,sspgrid%lam(1:nl_fit))
     idata(1:nl_fit)%wgt = linterp(tlam(1:datmax),&
          data(1:datmax)%wgt,sspgrid%lam(1:nl_fit))
     
     lo = MAX(l1(1),tlam(1))+50
     hi = MIN(l2(1),tlam(datmax))-50
     IF (lo.GE.hi) CYCLE

     CALL CONTNORMSPEC(sspgrid%lam,idata%flx,idata%err,lo,hi,dflx)
     !use a 5 Gyr Zsol SSP
     CALL CONTNORMSPEC(sspgrid%lam,10**sspgrid%logfkrpa(4,:),&
          idata%wgt*SQRT(10**sspgrid%logfkrpa(3,:)),lo,hi,mflx)

     i1 = MIN(MAX(locate(sspgrid%lam,lo),1),nl_fit-1)
     i2 = MIN(MAX(locate(sspgrid%lam,hi),2),nl_fit)
     tchi2(i) = SUM(idata(i1:i2)%flx**2/idata(i1:i2)%err**2*&
          (dflx(i1:i2)-mflx(i1:i2))**2) / (i2-i1)
     
     IF (tchi2(i).LT.chi2) THEN
        chi2    = tchi2(i)
        getvelz = tvz(i)
     ENDIF

  ENDDO

  !test to see if the solution is good
  !we take all points with delta(chi2/dof)<1 and
  !ask how large is the range in velocities
  tchi2 = tchi2 - MINVAL(tchi2)
  tvza  = getvelz
  DO i=1,nv
     IF (tchi2(i).LT.1.0) tvza(i)=tvz(i)
  ENDDO

  IF ((MAXVAL(tvza)-MINVAL(tvza)).GT.1E3) THEN
     WRITE(*,'("   Failed to find a redshift solution, setting velz=0.0")')
     getvelz = 0.0
  ENDIF



END FUNCTION GETVELZ
