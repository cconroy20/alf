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
  INTEGER :: i

  !---------------------------------------------------------------!
  !---------------------------------------------------------------!

  !convert to the proper units
  aspec  = spec*lsun/1E6*lam**2/clight/1E8/4/mypi/pc2cm**2

  DO i=1,3
     mag(i) = tsum(lam,aspec*fil(i,:)/lam)
     IF (mag(i).LE.0.0) THEN
        mag(i)=99.
     ELSE 
        mag(i) = -2.5*LOG10(mag(i))-48.60
     ENDIF
  ENDDO

  IF (PRESENT(mw)) THEN
     mass = getmass(msto,krpa_imf1,krpa_imf2,krpa_imf3)
  ELSE
     mass = getmass(msto,pos%imf1,pos%imf2,krpa_imf3)
  ENDIF

  DO i=1,3 
     IF (mag(i).EQ.99) THEN
        m2l(i) = 0.0
     ELSE
        m2l(i)  = mass/10**(2./5*(magsun(i)-mag(i)))
        IF (m2l(i).GT.50.0) m2l(i)=0.0
     ENDIF
  ENDDO

END SUBROUTINE GETM2L
