SUBROUTINE GETM2L(lam,spec,pos,m2l,mw)

  !compute mass-to-light ratio in several filters (AB mags)

  USE sfvars; USE sfutils, ONLY : tsum, getmass
  IMPLICIT NONE

  REAL(DP), DIMENSION(nl), INTENT(in) :: lam,spec
  TYPE(PARAMS), INTENT(in)   :: pos
  REAL(DP), DIMENSION(nfil), INTENT(out) :: m2l
  REAL(DP), DIMENSION(nfil) :: mag
  REAL(DP), DIMENSION(nl)   :: aspec
  REAL(DP) :: mass, mto=1.10, imf1=1.3,imf2=2.3,imf3=2.3
  INTEGER, OPTIONAL :: mw

  !---------------------------------------------------------------!
  !---------------------------------------------------------------!

  !convert to the proper units
  aspec  = spec*lsun/1E6*lam**2/clight/1E8/4/mypi/pc2cm**2

  mag(1) = -2.5*LOG10(tsum(lam,aspec*fil(1,:)/lam))-48.6
  mag(2) = -2.5*LOG10(tsum(lam,aspec*fil(2,:)/lam))-48.6
  mag(3) = 0.0 !-2.5*LOG10(tsum(lam,aspec*fil(3,:)/lam))-48.6

  IF (PRESENT(mw)) THEN
     mass = getmass(mto,imf1,imf2,imf3)
  ELSE
     mass = getmass(mto,pos%imf1,pos%imf2,imf3)
  ENDIF

  m2l  = mass / 10**(2./5*(magsun-mag))
  m2l(3) = 0.0

END SUBROUTINE GETM2L
