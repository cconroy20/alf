FUNCTION TSUM(xin,yin)

  !simple trapezoidal integration of tabulated function (xin,yin)

  USE sfvars
  IMPLICIT NONE

  REAL(DP), DIMENSION(:), INTENT(in) :: xin,yin
  REAL(DP) :: tsum
  INTEGER  :: nn

  !---------------------------------------------------------------!
  !---------------------------------------------------------------!

  nn = SIZE(xin)

  tsum = SUM( (xin(2:nn)-xin(1:nn-1)) * (yin(2:nn)+yin(1:nn-1))/2. )


END FUNCTION TSUM
