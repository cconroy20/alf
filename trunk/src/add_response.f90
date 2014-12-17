SUBROUTINE ADD_RESPONSE(spec,pos,range,dr,vr,solar,plus,minus)

  USE sfvars
  IMPLICIT NONE

  INTEGER, INTENT(in)  :: vr
  REAL(SP), INTENT(in) :: range
  REAL(DP), INTENT(in) :: pos,dr
  REAL(DP), DIMENSION(nl), INTENT(inout) :: spec
  REAL(DP), DIMENSION(nl) :: tmpr
  REAL(DP), DIMENSION(nl,nage_rfcn), INTENT(in) :: plus,solar
  REAL(DP), DIMENSION(nl,nage_rfcn), INTENT(in), OPTIONAL :: minus

  !---------------------------------------------------------------!
  !---------------------------------------------------------------!
 
  IF (PRESENT(minus)) THEN 

     IF (pos.GT.0.0) THEN
        tmpr(1:nl_fit) = dr*plus(1:nl_fit,vr+1)/solar(1:nl_fit,vr+1) + &
             (1-dr)*plus(1:nl_fit,vr)/solar(1:nl_fit,vr)
        spec(1:nl_fit) = spec(1:nl_fit) * (1 + (tmpr(1:nl_fit)-1)*pos/range)
     ELSE 
        tmpr(1:nl_fit) = dr*minus(1:nl_fit,vr+1)/solar(1:nl_fit,vr+1) + &
             (1-dr)*minus(1:nl_fit,vr)/solar(1:nl_fit,vr)
        spec(1:nl_fit) = spec(1:nl_fit) * (1 + (tmpr-1)*ABS(pos)/range)
     ENDIF

  ELSE
 
     tmpr(1:nl_fit) = dr*plus(1:nl_fit,vr+1)/solar(1:nl_fit,vr+1) + &
          (1-dr)*plus(1:nl_fit,vr)/solar(1:nl_fit,vr)
     spec(1:nl_fit) = spec(1:nl_fit) * (1 + (tmpr(1:nl_fit)-1)*pos/range)
 
  ENDIF

END SUBROUTINE ADD_RESPONSE
