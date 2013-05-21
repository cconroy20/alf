SUBROUTINE SET_WAVE_LIMITS(file)

  !define the wavelength boundaries used in the fit.  Each 
  !interval is fit separately.  If the red or blue limit
  !is beyond the data limit, then that interval is exluded from
  !the fit

  USE sfvars
  IMPLICIT NONE
  
  CHARACTER(50), INTENT(in) :: file

  !---------------------------------------------------------------!
  !---------------------------------------------------------------!

  IF (file(1:4).EQ.'sdss') THEN
     l1 = (/0.40,0.48,0.58,0.80/) * 1E4
     l2 = (/0.48,0.58,0.64,0.88/) * 1E4
  ELSE IF (file(1:1).EQ.'n'.OR.file(1:3).EQ.'m32'.OR.&
       file(1:1).EQ.'b') THEN
     l1   = (/0.40,0.46,0.800,0.963/) * 1E4
     l2   = (/0.46,0.55,0.892,1.015/) * 1E4
  ELSE IF (file(1:4).EQ.'ages'.OR.file(1:2).EQ.'SN') THEN
     l1   = (/0.40,0.46,0.6,0.7/) * 1E4
     l2   = (/0.46,0.55,0.7,0.8/) * 1E4
  ELSE
     l1 = (/0.40,0.46,0.6,0.7/) * 1E4
     l2 = (/0.46,0.535,0.7,0.8/) * 1E4
  ENDIF

END SUBROUTINE SET_WAVE_LIMITS
