SUBROUTINE GETM2L(msto,lam,spec,pos,m2l,mw)

  !compute mass-to-light ratio in several filters (AB mags)

  USE sfvars; USE sfutils, ONLY : tsum, getmass
  IMPLICIT NONE

  REAL(DP), DIMENSION(nl), INTENT(in) :: lam,spec
  REAL(DP), INTENT(in) :: msto
  TYPE(PARAMS), INTENT(in)   :: pos
  REAL(DP), DIMENSION(nfil), INTENT(out) :: m2l
  REAL(DP), DIMENSION(nfil) :: mag
  REAL(DP), DIMENSION(nl)   :: aspec
  REAL(DP) :: mass
  INTEGER, OPTIONAL :: mw

  !---------------------------------------------------------------!
  !---------------------------------------------------------------!

  !convert to the proper units
  aspec  = spec*lsun/1E6*lam**2/clight/1E8/4/mypi/pc2cm**2

  mag(1) = -2.5*LOG10(tsum(lam,aspec*fil(1,:)/lam))-48.6
  mag(2) = -2.5*LOG10(tsum(lam,aspec*fil(2,:)/lam))-48.6
  mag(3) = 0.0 !-2.5*LOG10(tsum(lam,aspec*fil(3,:)/lam))-48.6

  IF (PRESENT(mw)) THEN
     mass = getmass(msto,krpa_imf1,krpa_imf2,krpa_imf3)
  ELSE
     mass = getmass(msto,pos%imf1,pos%imf2,krpa_imf3)
  ENDIF

  m2l  = mass / 10**(2./5*(magsun-mag))
  m2l(3) = 0.0

END SUBROUTINE GETM2L
