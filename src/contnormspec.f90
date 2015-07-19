SUBROUTINE NPOLY(x,arr)
  USE nrtype
  IMPLICIT NONE
  REAL(SP),INTENT(IN) :: x
  REAL(SP), DIMENSION(:), INTENT(OUT) :: arr
  INTEGER :: i

  DO i=1,SIZE(arr)
     arr(i) = x**(i-1)
  ENDDO

END SUBROUTINE NPOLY

!---------------------------------------------------------------!
!---------------------------------------------------------------!

SUBROUTINE CONTNORMSPEC(lam,flx,err,il1,il2,flxout,coeff)

  !routine to continuum normalize a spectrum by a high-order
  !polynomial.  The order of the polynomial is determined by
  !n=(lam_max-lam_min)/100.  Only normalized over the input
  !min/max wavelength range

  USE alf_vars; USE nr, ONLY : locate, lfit, svdfit
  IMPLICIT NONE

  REAL(DP), DIMENSION(nl), INTENT(in) :: lam,flx,err
  REAL(DP), INTENT(in) :: il1,il2
  REAL(DP), DIMENSION(nl), INTENT(inout) :: flxout
  REAL(DP), DIMENSION(ncoeff), OPTIONAL :: coeff
  REAL(DP), PARAMETER :: buff=0.0
  REAL(DP), DIMENSION(nl) :: poly
  REAL(DP), DIMENSION(20) :: tcoeff=0.0
  LOGICAL, DIMENSION(20)  :: mask=.TRUE.
  REAL(DP), DIMENSION(20,20) :: covar
  REAL(DP) :: ml,chi2sqr
  INTEGER  :: i1,i2,npow ,i
  INTERFACE
     SUBROUTINE NPOLY(x,arr)
       USE nrtype
       IMPLICIT NONE
       REAL(SP),INTENT(IN) :: x
       REAL(SP), DIMENSION(:), INTENT(OUT) :: arr
     END SUBROUTINE NPOLY
  END INTERFACE

  !---------------------------------------------------------------!
  !---------------------------------------------------------------!

  !divide by a power-law of degree npow. one degree per 100A.
  !don't let things get out of hand (force Npow<=npolymax)
  npow = MIN(NINT((il2-il1)/100.0),npolymax)
  
  i1 = MIN(MAX(locate(lam,il1-buff),1),nl-1)
  i2 = MIN(MAX(locate(lam,il2+buff),2),nl)
  ml = (il1+il2)/2.0

  !simple linear least squares polynomial fit
  CALL LFIT(lam(i1:i2)-ml,flx(i1:i2),err(i1:i2),tcoeff(1:npow+1),&
       mask(1:npow+1),covar(1:npow+1,1:npow+1),chi2sqr,npoly)

  IF (PRESENT(coeff)) THEN
     coeff(1:npow+1) = tcoeff(1:npow+1)
  ENDIF

  poly = 0.0
  DO i=1,npow+1 
     poly = poly + tcoeff(i)*(lam-ml)**(i-1)
  ENDDO
  flxout = flx / poly



END SUBROUTINE CONTNORMSPEC
