SUBROUTINE GETM2L(lam,spec,pos,m2l,mw)

  !compute mass-to-light ratio in several filters (AB mags)

  USE alf_vars; USE alf_utils, ONLY : tsum, getmass
  IMPLICIT NONE

  REAL(DP), DIMENSION(nl), INTENT(in) :: lam,spec
  TYPE(PARAMS), INTENT(in)  :: pos
  REAL(DP), DIMENSION(nfil), INTENT(out) :: m2l
  INTEGER, OPTIONAL :: mw
  REAL(DP), DIMENSION(nfil) :: mag
  REAL(DP), DIMENSION(nl)   :: aspec
  REAL(DP) :: mass,msto
  INTEGER  :: i

  !---------------------------------------------------------------!
  !---------------------------------------------------------------!

  !convert to the proper units
  aspec  = spec*lsun/1E6*lam**2/clight/1E8/4/mypi/pc2cm**2

  msto = 10**(msto_t0+msto_t1*pos%logage) * &
       ( msto_z0 + msto_z1*pos%zh + msto_z2*pos%zh**2 )
  
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
        mass = getmass(pos%imf3,msto,pos%imf1,pos%imf1,krpa_imf3)
     ELSE IF (imf_type.EQ.3) THEN
        !double powerlaw index with variable low-mass cutoff
        mass = getmass(pos%imf3,msto,pos%imf1,pos%imf2,krpa_imf3)
     ELSE IF (imf_type.EQ.4) THEN
        !non-parametric IMF for 0.08-1.0 Msun; Salpeter slope at >1 Msun
        mass = getmass(imflo,msto,pos%imf1,pos%imf2,krpa_imf3,&
             pos%imf3,pos%imf4)
     ENDIF
  ENDIF

  !loop over the filters
  DO i=1,nfil
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
