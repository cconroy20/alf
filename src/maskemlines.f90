SUBROUTINE MASKEMLINES(zred,sigma)

  !routine to mask emission lines

  USE alf_vars
  IMPLICIT NONE

  REAL(DP), INTENT(in)      :: zred,sigma
  REAL(DP), DIMENSION(neml) :: wave0,wavel,waveh
  INTEGER :: i,j, m=2

  !---------------------------------------------------------------!

  !central wavelengths of (potentially) strong emission lines
  !convert to observed frame
  wave0 = emlines * (1+zred/clight*1E5)
  !mask within +/-m sigma of the central wavelength
  wavel = wave0 - m*wave0*sigma/clight*1E5
  waveh = wave0 + m*wave0*sigma/clight*1E5

  DO i=1,datmax
     DO j=1,neml
        IF (data(i)%lam.GE.wavel(j).AND.data(i)%lam.LE.waveh(j)) THEN
           data(i)%wgt = huge_number
           data(i)%err = huge_number
        ENDIF
     ENDDO
  ENDDO
  
END SUBROUTINE MASKEMLINES
