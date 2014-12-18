SUBROUTINE LINTERP3(xin,yin1,yin2,yin3,xout,yout1,yout2,yout3)

  !routine to linearly interpolate three functions with the 
  !same input x-axis to the same output x-axis.
  !this routine was created because the locate function is 
  !slow so we only have to call it once, rather than thrice.

  USE sfvars; USE nr, ONLY: locate
  IMPLICIT NONE
  REAL(DP), DIMENSION(:), INTENT(in)    :: xin,yin1,yin2,yin3,xout
  REAL(DP), DIMENSION(:), INTENT(inout) :: yout1,yout2,yout3
  INTEGER :: klo,n,n2,i
  REAL(DP) :: dx

  !---------------------------------------------------------------!
  !---------------------------------------------------------------!

  n   = SIZE(xin)
  n2  = SIZE(xout)

  DO i=1,n2
     klo = MAX(MIN(locate(xin(1:n),xout(i)),n-1),1)
     dx  = (xout(i)-xin(klo))/(xin(klo+1)-xin(klo))
     yout1(i) = yin1(klo)+dx*(yin1(klo+1)-yin1(klo))
     yout2(i) = yin2(klo)+dx*(yin2(klo+1)-yin2(klo))
     yout3(i) = yin3(klo)+dx*(yin3(klo+1)-yin3(klo))
  ENDDO

END SUBROUTINE LINTERP3
