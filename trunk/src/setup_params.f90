SUBROUTINE SETUP_PARAMS(pos,dstep,prlo,prhi,velz)

  !define the first position (pos), the step size (dstep),
  !and the lower and upper bounds on the priors (prlo, prhi)

  USE sfvars; USE sfutils, ONLY : myran
  IMPLICIT NONE

  TYPE(PARAMS), INTENT(inout) :: pos,prlo,prhi,dstep
  REAL(DP), OPTIONAL :: velz

  !---------------------------------------------------------------!
  !---------------------------------------------------------------!

  !setup the first position
  pos%logage    = myran()*0.2+0.8
  pos%feh       = myran()*0.2-0.1
  pos%ah        = myran()*0.2-0.1
  pos%nhe       = myran()*0.2-0.1
  pos%ch        = myran()*0.2-0.1
  pos%nh        = myran()*0.2-0.1
  pos%nah       = myran()*0.2-0.1
  pos%mgh       = myran()*0.2-0.1
  pos%sih       = myran()*0.2-0.1
  pos%kh        = myran()*0.2-0.1
  pos%cah       = myran()*0.2-0.1
  pos%tih       = myran()*0.2-0.1
  pos%vh        = myran()*0.2-0.1
  pos%crh       = myran()*0.2-0.1
  pos%mnh       = myran()*0.2-0.1
  pos%coh       = myran()*0.2-0.1
  pos%nih       = myran()*0.2-0.1
  pos%cuh       = myran()*0.2-0.1
  pos%rbh       = myran()*0.2-0.1
  pos%srh       = myran()*0.2-0.1
  pos%yh        = myran()*0.2-0.1
  pos%zrh       = myran()*0.2-0.1
  pos%bah       = myran()*0.2-0.1
  pos%euh       = myran()*0.2-0.1
  pos%teff      = myran()*20-10
  pos%imf1      = myran()*0.8-0.4 + 1.3
  pos%imf2      = myran()*0.8-0.4 + 2.3
  pos%logfy     = myran()*0.2-3
  pos%logm7g    = myran()*0.2-3
  pos%hotteff   = myran()*5+15
  pos%loghot    = myran()*0.2-3
  pos%chi2      = huge_number
  pos%sigma     = 200.0
  pos%sigma2    = 200.0
  pos%velz2     = 0.0 
  pos%logemnorm = -2.0
  pos%logcoeff  = -20.0
  pos%logcoeff(1) = -1.0

  IF (PRESENT(velz)) THEN
     pos%velz  = velz
  ELSE
     pos%velz  = 0.0
  ENDIF

  !setup the step size
  dstep%logage    = mcstep * 20 ! * 500
  dstep%feh       = mcstep
  dstep%ah        = mcstep * 10
  dstep%nhe       = mcstep * 200
  dstep%ch        = mcstep
  dstep%nh        = mcstep
  dstep%nah       = mcstep * 10
  dstep%mgh       = mcstep
  dstep%sih       = mcstep
  dstep%kh        = mcstep * 200
  dstep%cah       = mcstep
  dstep%tih       = mcstep
  dstep%vh        = mcstep * 20
  dstep%crh       = mcstep * 20
  dstep%mnh       = mcstep * 20
  dstep%coh       = mcstep * 20
  dstep%nih       = mcstep * 20
  dstep%cuh       = mcstep * 100
  dstep%rbh       = mcstep * 200
  dstep%srh       = mcstep * 200
  dstep%yh        = mcstep * 200
  dstep%zrh       = mcstep * 200
  dstep%bah       = mcstep * 200
  dstep%euh       = mcstep * 200
  dstep%teff      = mcstep * 800
  dstep%imf1      = mcstep * 100
  dstep%imf2      = mcstep * 100
  dstep%logfy     = mcstep * 100
  dstep%logm7g    = mcstep * 100
  dstep%hotteff   = mcstep * 100
  dstep%loghot    = mcstep * 100
  dstep%sigma     = mcstep * 100
  dstep%sigma2    = mcstep * 100
  dstep%velz      = mcstep * 50
  dstep%velz2     = mcstep * 50
  dstep%logemnorm = mcstep * 100
  dstep%logcoeff  = mcstep * 100

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
