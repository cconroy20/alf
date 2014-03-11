SUBROUTINE SETUP_PARAMS(pos,prlo,prhi,velz)

  !define the first position (pos), the step size (dstep),
  !and the lower and upper bounds on the priors (prlo, prhi)

  USE sfvars; USE sfutils, ONLY : myran
  IMPLICIT NONE

  TYPE(PARAMS), INTENT(inout) :: pos,prlo,prhi
  REAL(DP), OPTIONAL :: velz
  INTEGER :: i
  
  !---------------------------------------------------------------!
  !---------------------------------------------------------------!

  !setup the first position
  pos%logage    = myran()*0.2+0.8
  pos%feh       = myran()*1.0-0.5
  pos%ah        = myran()*1.0-0.5
  pos%nhe       = myran()*1.0-0.5
  pos%ch        = myran()*1.0-0.5
  pos%nh        = myran()*1.0-0.5
  pos%nah       = myran()*1.0-0.5
  pos%mgh       = myran()*1.0-0.5
  pos%sih       = myran()*1.0-0.5
  pos%kh        = myran()*1.0-0.5
  pos%cah       = myran()*1.0-0.5
  pos%tih       = myran()*1.0-0.5
  pos%vh        = myran()*1.0-0.5
  pos%crh       = myran()*1.0-0.5
  pos%mnh       = myran()*1.0-0.5
  pos%coh       = myran()*1.0-0.5
  pos%nih       = myran()*1.0-0.5
  pos%cuh       = myran()*1.0-0.5
  pos%rbh       = myran()*1.0-0.5
  pos%srh       = myran()*1.0-0.5
  pos%yh        = myran()*1.0-0.5
  pos%zrh       = myran()*1.0-0.5
  pos%bah       = myran()*1.0-0.5
  pos%euh       = myran()*1.0-0.5
  pos%teff      = myran()*80-40
  pos%imf1      = myran()*0.8-0.4 + 1.3
  pos%imf2      = myran()*0.8-0.4 + 2.3
  pos%logfy     = myran()*2-3
  pos%logm7g    = myran()*2-3
  pos%hotteff   = myran()*5+15
  pos%loghot    = myran()*2-3
  pos%chi2      = huge_number
  pos%sigma     = myran()*100+100.0
  pos%sigma2    = myran()*100+100.0
  pos%velz2     = myran()*10-5
  DO i=1,neml
     pos%logemnorm(i) = myran()*2-3
  ENDDO

  IF (PRESENT(velz)) THEN
     pos%velz  = velz + (myran()*10-5)
  ELSE
     pos%velz  = myran()*10-5
  ENDIF


  !priors (low)
  prlo%logage    = LOG10(0.5)
  prlo%feh       = -1.0
  prlo%ah        = -1.0
  prlo%nhe       = -1.0
  prlo%ch        = -1.0
  prlo%nh        = -1.0
  prlo%nah       = -1.0
  prlo%mgh       = -1.0
  prlo%sih       = -1.0
  prlo%kh        = -1.0
  prlo%cah       = -1.0
  prlo%tih       = -1.0
  prlo%vh        = -1.0
  prlo%crh       = -1.0
  prlo%mnh       = -1.0
  prlo%coh       = -1.0
  prlo%nih       = -1.0
  prlo%cuh       = -1.0
  prlo%rbh       = -1.0
  prlo%srh       = -1.0
  prlo%yh        = -1.0
  prlo%zrh       = -1.0
  prlo%bah       = -1.0
  prlo%euh       = -1.0
  prlo%teff      = -100.0
  prlo%imf1      = 0.5
  prlo%imf2      = 0.5
  prlo%logfy     = -5.0
  prlo%logm7g    = -5.0
  prlo%hotteff   = 10.0
  prlo%loghot    = -5.0
  prlo%sigma     = 20.0
  prlo%sigma2    = 20.0
  prlo%velz      = -1E4
  prlo%velz2     = -1E4
  prlo%logemnorm = -5.0
  prlo%logcoeff  = -30.

  !priors (high)
  prhi%logage    = LOG10(20.0)
  prhi%feh       = 1.0
  prhi%ah        = 1.0
  prhi%nhe       = 1.0
  prhi%ch        = 1.0
  prhi%nh        = 1.0
  prhi%nah       = 1.0
  prhi%mgh       = 1.0
  prhi%sih       = 1.0
  prhi%kh        = 1.0
  prhi%cah       = 1.0
  prhi%tih       = 1.0
  prhi%vh        = 1.0
  prhi%crh       = 1.0
  prhi%mnh       = 1.0
  prhi%coh       = 1.0
  prhi%nih       = 1.0
  prhi%cuh       = 1.0
  prhi%rbh       = 1.0
  prhi%srh       = 1.0
  prhi%yh        = 1.0
  prhi%zrh       = 1.0
  prhi%bah       = 1.0
  prhi%euh       = 1.0
  prhi%teff      = 100.0
  prhi%imf1      = 3.5
  prhi%imf2      = 3.5
  prhi%logfy     = -0.1
  prhi%logm7g    = -0.1
  prhi%hotteff   = 30.0
  prhi%loghot    = -0.1
  prhi%sigma     = 1E4
  prhi%sigma2    = 1E4
  prhi%velz      = 1E4
  prhi%velz2     = 1E4
  prhi%logemnorm = 1.0
  prhi%logcoeff  = 10.

END SUBROUTINE SETUP_PARAMS
