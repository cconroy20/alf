FUNCTION GETMASS(mto,imf1,imf2,imf3)

  !compute mass in stars and remnants (normalized to 1 Msun at t=0)
  !assume an IMF that runs from 0.1 to 100 Msun.

  USE sfvars
  USE nr, ONLY : qromb; USE sfutils, ONLY : imf
  IMPLICIT NONE

  !turnoff mass
  REAL(DP), INTENT(in) :: mto
  REAL(DP), INTENT(in) :: imf1,imf2,imf3
  REAL(DP) :: imfnorm, getmass
  REAL(DP), PARAMETER :: bhlim=40.0,nslim=8.5

  !---------------------------------------------------------------!
  !---------------------------------------------------------------!

  imf_alpha(1) = imf1
  imf_alpha(2) = imf2
  imf_alpha(3) = imf3

  getmass  = 0.d0
  imf_flag = 0

  !normalize the weights so that 1 Msun formed at t=0
  imf_flag = imf_flag+10
  imfnorm  = qromb(imf,imflo,imfhi)
  imf_flag = imf_flag-10

  !stars still alive
  imf_flag = imf_flag+10
  getmass  = getmass + qromb(imf,imflo,mto)/imfnorm
  imf_flag = imf_flag-10

  !BH remnants
  !40<M<imf_up leave behind a 0.5*M BH
  imf_flag = imf_flag+10
  getmass  = getmass + 0.5*qromb(imf,bhlim,imfhi)/imfnorm
  imf_flag = imf_flag-10

  !NS remnants
  !8.5<M<40 leave behind 1.4 Msun NS
  getmass = getmass + 1.4*qromb(imf,nslim,bhlim)/imfnorm

  !WD remnants
  !M<8.5 leave behind 0.077*M+0.48 WD
  getmass  = getmass + 0.48*qromb(imf,mto,nslim)/imfnorm
  imf_flag = imf_flag+10
  getmass  = getmass + 0.077*qromb(imf,mto,nslim)/imfnorm
  imf_flag = imf_flag-10

  RETURN

END FUNCTION GETMASS
