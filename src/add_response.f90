SUBROUTINE ADD_RESPONSE(spec,pos,range,dr,vr,dm,vm,solar,plus,minus)

  !perform bilinear interpolation over the age and metallicity-
  !dependent response functions

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

        tmpr = dr*dm*plus(:,vr+1,vm+1)/solar(:,vr+1,vm+1) + &
             (1-dr)*dm*plus(:,vr,vm+1)/solar(:,vr,vm+1) + &
             dr*(1-dm)*plus(:,vr+1,vm)/solar(:,vr+1,vm) + &
             (1-dr)*(1-dm)*plus(:,vr,vm)/solar(:,vr,vm)

        spec = spec * (1+(tmpr-1)*pos/range)

     ELSE 

        tmpr = dr*dm*minus(:,vr+1,vm+1)/solar(:,vr+1,vm+1) + &
             (1-dr)*dm*minus(:,vr,vm+1)/solar(:,vr,vm+1) + &
             dr*(1-dm)*minus(:,vr+1,vm)/solar(:,vr+1,vm) + &
             (1-dr)*(1-dm)*minus(:,vr,vm)/solar(:,vr,vm)

        spec = spec * (1+(tmpr-1)*ABS(pos)/range)

     ENDIF

  ELSE
 
     tmpr = dr*dm*plus(:,vr+1,vm+1)/solar(:,vr+1,vm+1) + &
          (1-dr)*dm*plus(:,vr,vm+1)/solar(:,vr,vm+1) + &
          dr*(1-dm)*plus(:,vr+1,vm)/solar(:,vr+1,vm) + &
          (1-dr)*(1-dm)*plus(:,vr,vm)/solar(:,vr,vm)

     spec = spec * (1+(tmpr-1)*pos/range)
 
  ENDIF

END SUBROUTINE ADD_RESPONSE
