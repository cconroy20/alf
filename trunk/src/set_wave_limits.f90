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

  IF (file(1:4).EQ.'sdss'.OR.file(1:3).EQ.'ucd') THEN
     l1 = (/0.40,0.47,0.58,0.80/) * 1E4
     l2 = (/0.47,0.58,0.64,0.88/) * 1E4
  ELSE IF (file(1:1).EQ.'n'.OR.file(1:3).EQ.'m32'.OR.&
       file(1:1).EQ.'b'.OR.file(1:5).EQ.'manga'.OR.&
       file(1:5).EQ.'deep2') THEN
     l1  = (/0.40,0.46,0.800,0.963/) * 1E4
     l2  = (/0.46,0.55,0.892,1.015/) * 1E4
  ELSE IF (file(1:4).EQ.'ages'.OR.file(1:2).EQ.'SN'&
       .OR.file(1:7).EQ.'mosfire') THEN
     l1  = (/0.40,0.47,0.6,0.7/) * 1E4
     l2   = (/0.47,0.55,0.7,0.8/) * 1E4
  ELSE IF (file(1:3).EQ.'eso') THEN
     l1  = (/0.4025,0.45,0.55,0.8/) * 1E4
     l2  = (/0.45,0.55,0.6,0.885/) * 1E4
  ELSE IF (file(1:5).EQ.'usher') THEN
     l1  = (/0.4025,0.45,0.6530,0.8160/) * 1E4
     l2  = (/0.45,0.55,0.6720,0.8870/) * 1E4
  ELSE IF (file(1:4).EQ.'keck') THEN
     l1  = (/0.40,0.47,0.544,0.7/) * 1E4
     l2  = (/0.47,0.544,0.7,0.8/) * 1E4	
  ELSE
     l1 = (/0.40,0.46,0.6,0.7/) * 1E4
     l2 = (/0.46,0.54,0.7,0.8/) * 1E4
  ENDIF

END SUBROUTINE SET_WAVE_LIMITS
