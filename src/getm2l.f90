SUBROUTINE GETM2L(msto,lam,spec,pos,m2l,mw)

  !compute mass-to-light ratio in several filters (AB mags)

  USE alf_vars; USE alf_utils, ONLY : tsum, getmass
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
     mass = getmass(imflo,msto,krpa_imf1,krpa_imf2,krpa_imf3)
  ELSE
     IF (imf_type.EQ.0) THEN
        !single power-law IMF with a fixed lower-mass cutoff
        mass = getmass(imflo,msto,pos%imf1,pos%imf1,krpa_imf3)
     ELSE IF (imf_type.EQ.1) THEN
        !double power-law IMF with a fixed lower-mass cutoff
        mass = getmass(imflo,msto,pos%imf1,pos%imf2,krpa_imf3)
     ELSE IF (imf_type.EQ.2) THEN
        !single powerlaw index with variable low-mass cutoff
        mass = getmass(pos%imf2,msto,pos%imf1,pos%imf1,krpa_imf3)
     ELSE IF (imf_type.EQ.3) THEN
        !double powerlaw index with variable low-mass cutoff
        mass = getmass(pos%imf3,msto,pos%imf1,pos%imf2,krpa_imf3)
     ENDIF
  ENDIF

  !loop over the three filters
  DO i=1,3
     mag(i) = tsum(lam,aspec*filters(:,i)/lam)
     IF (mag(i).LE.0.0) THEN
        m2l(i) = 0.0
     ELSE 
        mag(i) = -2.5*LOG10(mag(i))-48.60
        m2l(i) = mass/10**(2./5*(magsun(i)-mag(i)))
        IF (m2l(i).GT.100.) m2l(i)=0.0
     ENDIF
  ENDDO

END SUBROUTINE GETM2L
