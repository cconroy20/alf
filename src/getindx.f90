FUNCTION INTIND(lam,spec,lo,hi)

  !perform integral over spectrum for index computation

  USE alf_vars
  USE alf_utils, ONLY : tsum
  USE nr, ONLY : locate
  IMPLICIT NONE

  INTEGER :: ll1,ll2,nn
  REAL(DP), INTENT(in), DIMENSION(:) :: lam, spec
  REAL(DP), INTENT(in) :: lo,hi
  REAL(DP) :: f1,f2,intind

  !---------------------------------------------------------------!
  !---------------------------------------------------------------!

  nn = SIZE(lam)

  !take care of the ends
  ll1 = MAX(MIN(locate(lam(1:nn),lo),nn-1),1)
  ll2 = MAX(MIN(locate(lam(1:nn),hi),nn-1),2)
  f1 = (spec(ll1+1)-spec(ll1))/(lam(ll1+1)-lam(ll1))*&
       (lo-lam(ll1))+spec(ll1)
  f2 = (spec(ll2+1)-spec(ll2))/(lam(ll2+1)-lam(ll2))*&
       (hi-lam(ll2))+spec(ll2)

  IF (ll1.EQ.ll2) THEN
     intind = (f2+f1)/2.*(hi-lo)
  ELSE
     intind = TSUM(lam(ll1+1:ll2),spec(ll1+1:ll2))
     intind = intind + (lam(ll1+1)-lo)*(f1+spec(ll1+1))/2.
     intind = intind + (hi-lam(ll2))*(f2+spec(ll2))/2.
  ENDIF

END FUNCTION INTIND

!------------------------------------------------------------!
!------------------------------------------------------------!
!------------------------------------------------------------!

SUBROUTINE GETINDX(lambda,spec,indices)

  !routine to calculate indices from an input spectrum
  !indices are defined in fsps/data/allindices.dat

  USE alf_vars
  USE alf_utils, ONLY : intind
  IMPLICIT NONE

  INTEGER :: j,nn
  REAL(DP), INTENT(in), DIMENSION(:) :: spec,lambda
  REAL(DP), INTENT(inout), DIMENSION(nindx) :: indices
  REAL(DP) :: intfifc,cb,cr,lr,lb

  !---------------------------------------------------------------!
  !---------------------------------------------------------------!

  indices = 999.
  nn = SIZE(lambda)

  DO j=1,nindx

     !blue continuum
     cb = intind(lambda,spec,indxdef(3,j),indxdef(4,j))
     cb = cb / (indxdef(4,j)-indxdef(3,j))     
     lb = (indxdef(3,j)+indxdef(4,j))/2.

     IF (indxdef(7,j).NE.4) THEN 

        !red continuum 
        cr = intind(lambda,spec,indxdef(5,j),indxdef(6,j))
        cr = cr / (indxdef(6,j)-indxdef(5,j))
        lr = (indxdef(5,j)+indxdef(6,j))/2.
        
        !compute integral(fi/fc)
        !NB: fc here is a linear interpolation between the red and blue.
        intfifc = intind(lambda,spec/((cr-cb)/(lr-lb)*(lambda-lb)+cb),&
             indxdef(1,j),indxdef(2,j))

     ELSE
        intfifc = intind(lambda,spec,indxdef(1,j),indxdef(2,j))
        intfifc = intfifc/(indxdef(2,j)-indxdef(1,j))
     ENDIF

     IF (indxdef(7,j).EQ.1.) THEN
        !compute magnitude
        indices(j) = -2.5*LOG10(intfifc/(indxdef(2,j)-indxdef(1,j)))
     ELSE IF (indxdef(7,j).EQ.2.) THEN
        !compute EW (in Ang)
        indices(j) = (indxdef(2,j)-indxdef(1,j)) - intfifc
     ELSE IF (indxdef(7,j).EQ.3.) THEN
        !compute Dn4000
        !NB: this only works with cr and cb computed above b/c 
        !the wavelength intervals for blue and red are the 
        !same for these indices, otherwise a slightly 
        !different cr and cb would have to be computed
        indices(j) = cr/cb
     ELSE IF (indxdef(7,j).EQ.4.) THEN
        !compute magnitude from a flux ratio
        indices(j) = -2.5*LOG10(intfifc/cb)
     ENDIF

     !set dummy values for indices defined off of the wavelength grid
     IF (indxdef(6,j).GT.lambda(nn)) indices(j) = 999.0
     IF (indxdef(3,j).LT.lambda(1))  indices(j) = 999.0

  ENDDO

END SUBROUTINE GETINDX
