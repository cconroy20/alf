FUNCTION LINTERP(xin,yin,xout)

  !routine to linearly interpolate a function yin(xin) at xout

  USE sfvars; USE nr, ONLY: locate
  IMPLICIT NONE
  REAL(DP), DIMENSION(:), INTENT(in) :: xin,yin,xout
  REAL(DP), DIMENSION(SIZE(xout)) :: linterp
  INTEGER :: klo,n,n2,i

  !---------------------------------------------------------------!
  !---------------------------------------------------------------!

  n   = SIZE(xin)
  n2  = SIZE(xout)

  DO i=1,n2
     klo = MAX(MIN(locate(xin(1:n),xout(i)),n-1),1)
     linterp(i) = yin(klo) + (yin(klo+1)-yin(klo))*&
          (xout(i)-xin(klo))/(xin(klo+1)-xin(klo))
  ENDDO

END FUNCTION LINTERP
