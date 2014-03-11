FUNCTION GETVELZ()

  !function to estimate the recession velocity
  !this routine is used to get a first-guess at the velocity
  !so that the subsequent Powell minimization (in specfit)
  !coverges faster.  uses 4100<lambda<5000A

  USE sfvars; USE nr, ONLY : locate
  USE sfutils, ONLY : linterp,contnormspec
  IMPLICIT NONE

  REAL(DP) :: getvelz, chi2,lo,hi,tchi2,tvz
  INTEGER, PARAMETER :: nv=500
  INTEGER :: i,i1,i2
  REAL(DP), DIMENSION(nl) :: mflx,dflx
  TYPE(TDATA), DIMENSION(nl) :: idata
  REAL(DP), DIMENSION(ndat) :: tlam

  !------------------------------------------------------!

  lo = l1(1)+100
  hi = l2(1)-100

  chi2 = huge_number

  DO i=1,nv

     tvz = REAL(i)*20-1E3

     !de-redshift the data and interpolate to model wave array
     tlam      = data%lam / (1+tvz/clight*1E5)
     idata%flx = linterp(tlam,data%flx,sspgrid%lam)
     idata%err = linterp(tlam,data%err,sspgrid%lam)
     idata%wgt = linterp(tlam,data%wgt,sspgrid%lam)
     
     CALL CONTNORMSPEC(sspgrid%lam,idata%flx,idata%wgt,lo,hi,dflx)
     CALL CONTNORMSPEC(sspgrid%lam,10**sspgrid%logfkrpa(nage,:),&
          idata%wgt,lo,hi,mflx)

     i1 = MIN(MAX(locate(sspgrid%lam,lo),1),nl-1)
     i2 = MIN(MAX(locate(sspgrid%lam,hi),2),nl)
     tchi2 = SUM(idata(i1:i2)%flx**2/idata(i1:i2)%err**2*&
          (dflx(i1:i2)-mflx(i1:i2))**2)

     IF (tchi2.LT.chi2) THEN
        chi2    = tchi2
        getvelz = tvz
     ENDIF

  ENDDO


END FUNCTION GETVELZ
