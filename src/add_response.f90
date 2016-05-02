SUBROUTINE ADD_RESPONSE(spec,pos,range,dr,vr,dm,vm,solar,plus,minus)

  USE alf_vars
  IMPLICIT NONE

  INTEGER, INTENT(in)  :: vr,vm
  REAL(SP), INTENT(in) :: range
  REAL(DP), INTENT(in) :: pos,dr,dm
  REAL(DP), DIMENSION(nl), INTENT(inout) :: spec
  REAL(DP), DIMENSION(nl) :: tmpr
  REAL(DP), DIMENSION(nl,nage_rfcn,nzmet), INTENT(in) :: plus,solar
  REAL(DP), DIMENSION(nl,nage_rfcn,nzmet), INTENT(in), OPTIONAL :: minus

  !---------------------------------------------------------------!
  !---------------------------------------------------------------!
 
  !in this case we have both + and - response functions
  IF (PRESENT(minus)) THEN 

     IF (pos.GT.0.0) THEN

        tmpr(1:nl_fit) = &
             dr*dm*plus(1:nl_fit,vr+1,vm+1)/solar(1:nl_fit,vr+1,vm+1) + &
             (1-dr)*dm*plus(1:nl_fit,vr,vm+1)/solar(1:nl_fit,vr,vm+1) + &
             dr*(1-dm)*plus(1:nl_fit,vr+1,vm)/solar(1:nl_fit,vr+1,vm) + &
             (1-dr)*(1-dm)*plus(1:nl_fit,vr,vm)/solar(1:nl_fit,vr,vm)

        spec(1:nl_fit) = spec(1:nl_fit) * (1+(tmpr(1:nl_fit)-1)*pos/range)

     ELSE 

        tmpr(1:nl_fit) = &
             dr*dm*minus(1:nl_fit,vr+1,vm+1)/solar(1:nl_fit,vr+1,vm+1) + &
             (1-dr)*dm*minus(1:nl_fit,vr,vm+1)/solar(1:nl_fit,vr,vm+1) + &
             dr*(1-dm)*minus(1:nl_fit,vr+1,vm)/solar(1:nl_fit,vr+1,vm) + &
             (1-dr)*(1-dm)*minus(1:nl_fit,vr,vm)/solar(1:nl_fit,vr,vm)

        spec(1:nl_fit) = spec(1:nl_fit) * (1+(tmpr(1:nl_fit)-1)*ABS(pos)/range)

     ENDIF

  ELSE
 
     tmpr(1:nl_fit) = &
          dr*dm*plus(1:nl_fit,vr+1,vm+1)/solar(1:nl_fit,vr+1,vm+1) + &
          (1-dr)*dm*plus(1:nl_fit,vr,vm+1)/solar(1:nl_fit,vr,vm+1) + &
          dr*(1-dm)*plus(1:nl_fit,vr+1,vm)/solar(1:nl_fit,vr+1,vm) + &
          (1-dr)*(1-dm)*plus(1:nl_fit,vr,vm)/solar(1:nl_fit,vr,vm)

     spec(1:nl_fit) = spec(1:nl_fit) * (1+(tmpr(1:nl_fit)-1)*pos/range)
 
  ENDIF

END SUBROUTINE ADD_RESPONSE
