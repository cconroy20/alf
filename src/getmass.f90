FUNCTION GETMASS(mlo,mto,imf1,imf2,imf3)

  !compute mass in stars and remnants (normalized to 1 Msun at t=0)
  !assume an IMF that runs from 0.08 to 100 Msun.

  USE alf_vars
  IMPLICIT NONE

  !turnoff mass
  REAL(DP), INTENT(in) :: mlo,mto,imf1,imf2,imf3
  REAL(DP) :: imfnorm, getmass
  REAL(DP), PARAMETER :: bhlim=40.0,nslim=8.5
  REAL(DP) :: m2=0.5,m3=1.0

  !---------------------------------------------------------------!
  !---------------------------------------------------------------!

  IF (mlo.GT.m2) THEN
     WRITE(*,*) 'GETMASS ERROR: mlo>m2'
     STOP
  ENDIF

  getmass  = 0.d0

  !normalize the weights so that 1 Msun formed at t=0
  imfnorm = (m2**(-imf1+2)-mlo**(-imf1+2))/(-imf1+2) + &
       m2**(-imf1+imf2)*(m3**(-imf2+2)-m2**(-imf2+2))/(-imf2+2) + &
       m2**(-imf1+imf2)*(imfhi**(-imf3+2)-m3**(-imf3+2))/(-imf3+2)

  !stars still alive
  getmass = getmass + (m2**(-imf1+2)-mlo**(-imf1+2))/(-imf1+2)
  IF (mto.LT.m3) THEN
     getmass = getmass + m2**(-imf1+imf2)*(mto**(-imf2+2)-m2**(-imf2+2))/(-imf2+2)
  ELSE
     getmass = getmass + m2**(-imf1+imf2)*(m3**(-imf2+2)-m2**(-imf2+2))/(-imf2+2) + &
          m2**(-imf1+imf2)*(mto**(-imf3+2)-m3**(-imf3+2))/(-imf3+2)
  ENDIF
  getmass = getmass/imfnorm


  !BH remnants
  !40<M<imf_up leave behind a 0.5*M BH
  getmass = getmass + &
      0.5*m2**(-imf1+imf2)*(imfhi**(-imf3+2)-bhlim**(-imf3+2))/(-imf3+2)/imfnorm

  !NS remnants
  !8.5<M<40 leave behind 1.4 Msun NS
  getmass = getmass + &
        1.4*m2**(-imf1+imf2)*(bhlim**(-imf3+1)-nslim**(-imf3+1))/(-imf3+1)/imfnorm

  !WD remnants
  !M<8.5 leave behind 0.077*M+0.48 WD
  IF (mto.LT.m3) THEN
     getmass = getmass + &
          0.48*m2**(-imf1+imf2)*(nslim**(-imf3+1)-m3**(-imf3+1))/(-imf3+1)/imfnorm
     getmass = getmass + &
          0.48*m2**(-imf1+imf2)*(m3**(-imf2+1)-mto**(-imf2+1))/(-imf2+1)/imfnorm
     getmass = getmass + &
          0.077*m2**(-imf1+imf2)*(nslim**(-imf3+2)-m3**(-imf3+2))/(-imf3+2)/imfnorm
     getmass = getmass + &
          0.077*m2**(-imf1+imf2)*(m3**(-imf2+2)-mto**(-imf2+2))/(-imf2+2)/imfnorm
  ELSE
     getmass = getmass + &
          0.48*m2**(-imf1+imf2)*(nslim**(-imf3+1)-mto**(-imf3+1))/(-imf3+1)/imfnorm
     getmass = getmass + &
          0.077*m2**(-imf1+imf2)*(nslim**(-imf3+2)-mto**(-imf3+2))/(-imf3+2)/imfnorm
  ENDIF

  RETURN

END FUNCTION GETMASS
