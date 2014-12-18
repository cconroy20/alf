SUBROUTINE GETM2L(msto,lam,spec,pos,m2l,mw)

  !compute mass-to-light ratio in several filters (AB mags)

  USE sfvars; USE sfutils, ONLY : tsum, getmass
  IMPLICIT NONE

  REAL(DP), DIMENSION(nl), INTENT(in) :: lam,spec
  REAL(DP), INTENT(in)      :: msto
  TYPE(PARAMS), INTENT(in)  :: pos
  REAL(DP), DIMENSION(nfil), INTENT(out) :: m2l
  INTEGER, OPTIONAL :: mw
  REAL(DP), DIMENSION(nfil) :: mag
  REAL(DP), DIMENSION(nl)   :: aspec
  REAL(DP) :: mass
  INTEGER  :: i

  !---------------------------------------------------------------!
  !---------------------------------------------------------------!

  !convert to the proper units
  aspec  = spec*lsun/1E6*lam**2/clight/1E8/4/mypi/pc2cm**2

  IF (PRESENT(mw)) THEN
     mass = getmass(msto,krpa_imf1,krpa_imf2,krpa_imf3)
  ELSE
     mass = getmass(msto,pos%imf1,pos%imf2,krpa_imf3)
  ENDIF

  !loop over the three filters
  DO i=1,3
     mag(i) = tsum(lam,aspec*filters(:,i)/lam)
     IF (mag(i).LE.0.0) THEN
        m2l(i) = 0.0
     ELSE 
        mag(i) = -2.5*LOG10(mag(i))-48.60
        m2l(i) = mass/10**(2./5*(magsun(i)-mag(i)))
     ENDIF
  ENDDO

END SUBROUTINE GETM2L
