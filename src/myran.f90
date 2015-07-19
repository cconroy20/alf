FUNCTION MYRAN()

  !turns ran1 into a function, rather than a subroutine

  USE alf_vars
  USE nr, ONLY : ran1
  IMPLICIT NONE

  REAL(DP) :: myran

  !---------------------------------------------------------------!
  !---------------------------------------------------------------!

  CALL ran1(myran)

END FUNCTION MYRAN
