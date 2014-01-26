SUBROUTINE ADD_RESPONSE(spec,pos,range,dr,vr,solar,plus,minus)

  USE sfvars
  IMPLICIT NONE

  INTEGER, INTENT(in)  :: vr
  REAL(SP), INTENT(in) :: range
  REAL(DP), INTENT(in) :: pos,dr
  REAL(DP), DIMENSION(nl), INTENT(inout) :: spec
  REAL(DP), DIMENSION(nl) :: tmp,tmpr
  REAL(DP), DIMENSION(nage_rfcn,nl), INTENT(in) :: plus,solar
  REAL(DP), DIMENSION(nage_rfcn,nl), INTENT(in), OPTIONAL :: minus

  !---------------------------------------------------------------!
  !---------------------------------------------------------------!
 
  IF (PRESENT(minus)) THEN 

     IF (pos.GT.0.0) THEN
        tmpr = dr*plus(vr+1,:)/solar(vr+1,:) + (1-dr)*plus(vr,:)/solar(vr,:)
        tmp  = 1 + (tmpr-1)*pos/range
        spec = spec * tmp
     ELSE IF (pos.LT.0.0) THEN     
        tmpr = dr*minus(vr+1,:)/solar(vr+1,:) + (1-dr)*minus(vr,:)/solar(vr,:)
        tmp  = 1 + (tmpr-1)*ABS(pos)/range
        spec = spec * tmp
     ENDIF

  ELSE
 
     tmpr = dr*plus(vr+1,:)/solar(vr+1,:) + (1-dr)*plus(vr,:)/solar(vr,:)
     tmp  = 1 + (tmpr-1)*pos/range
     spec = spec * tmp
 
  ENDIF


END SUBROUTINE ADD_RESPONSE
