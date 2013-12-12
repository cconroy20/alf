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
  pos%age  = myran()*4+6
  pos%feh  = myran()*0.2-0.1
  pos%afe  = myran()*0.2-0.1
  pos%nhe  = myran()*0.2-0.1
  pos%cfe  = myran()*0.2-0.1
  pos%nfe  = myran()*0.2-0.1
  pos%nafe = myran()*0.2-0.1
  pos%mgfe = myran()*0.2-0.1
  pos%sife = myran()*0.2-0.1
  pos%kfe  = myran()*0.2-0.1
  pos%cafe = myran()*0.2-0.1
  pos%tife = myran()*0.2-0.1
  pos%vfe  = myran()*0.2-0.1
  pos%crfe = myran()*0.2-0.1
  pos%mnfe = myran()*0.2-0.1
  pos%cofe = myran()*0.2-0.1
  pos%nife = myran()*0.2-0.1
  pos%cufe = myran()*0.2-0.1
  pos%rbfe = myran()*0.2-0.1
  pos%srfe = myran()*0.2-0.1
  pos%yfe  = myran()*0.2-0.1
  pos%zrfe = myran()*0.2-0.1
  pos%bafe = myran()*0.2-0.1
  pos%eufe = myran()*0.2-0.1
  pos%teff = myran()*20-10
  pos%imf1 = myran()*0.8-0.4 + 1.3
  pos%imf2 = myran()*0.8-0.4 + 2.3
  pos%logfy   = myran()*0.2-3
  pos%logm7g  = myran()*0.2-3
  pos%hotteff = myran()*5+15
  pos%loghot  = myran()*0.2-3
  pos%chi2    = huge_number
  pos%sigma   = 200.0
  pos%sigma2  = 200.0
  pos%velz2   = 0.0 
  pos%logemnorm = -2.0

  IF (PRESENT(velz)) THEN
     pos%velz  = velz
  ELSE
     pos%velz  = 0.0
  ENDIF

  !setup the step size
  dstep%age  = mcstep * 50 * 10
  dstep%feh  = mcstep
  dstep%afe  = mcstep * 10
  dstep%nhe  = mcstep * 200
  dstep%cfe  = mcstep
  dstep%nfe  = mcstep
  dstep%nafe = mcstep * 10
  dstep%mgfe = mcstep
  dstep%sife = mcstep
  dstep%kfe  = mcstep * 200
  dstep%cafe = mcstep
  dstep%tife = mcstep
  dstep%vfe  = mcstep * 20
  dstep%crfe = mcstep * 20
  dstep%mnfe = mcstep * 20
  dstep%cofe = mcstep * 20
  dstep%nife = mcstep * 20
  dstep%cufe = mcstep * 100
  dstep%rbfe = mcstep * 200
  dstep%srfe = mcstep * 200
  dstep%yfe  = mcstep * 200
  dstep%zrfe = mcstep * 200
  dstep%bafe = mcstep * 200
  dstep%eufe = mcstep * 200
  dstep%teff = mcstep * 800
  dstep%imf1 = mcstep * 100
  dstep%imf2 = mcstep * 100
  dstep%logfy   = mcstep * 100
  dstep%logm7g  = mcstep * 100
  dstep%hotteff = mcstep * 100
  dstep%loghot  = mcstep * 100
  dstep%sigma   = mcstep * 100
  dstep%sigma2  = mcstep * 100
  dstep%velz    = mcstep * 50
  dstep%velz2   = mcstep * 50
  dstep%logemnorm = mcstep * 100

  !priors (low)
  prlo%age  = 0.5
  prlo%feh  = -0.6
  prlo%afe  = -0.6
  prlo%nhe  = -0.6
  prlo%cfe  = -0.6
  prlo%nfe  = -0.6
  prlo%nafe = -0.6
  prlo%mgfe = -0.6
  prlo%sife = -0.6
  prlo%kfe  = -0.6
  prlo%cafe = -0.6
  prlo%tife = -0.6
  prlo%vfe  = -0.6
  prlo%crfe = -0.6
  prlo%mnfe = -0.6
  prlo%cofe = -0.6
  prlo%nife = -0.6
  prlo%cufe = -0.6
  prlo%rbfe = -0.6
  prlo%srfe = -0.6
  prlo%yfe  = -0.6
  prlo%zrfe = -0.6
  prlo%bafe = -0.6
  prlo%eufe = -0.6
  prlo%teff = -100.0
  prlo%imf1 = 0.5
  prlo%imf2 = 0.5
  prlo%logfy   = -5.0
  prlo%logm7g  = -5.0
  prlo%hotteff = 10.0
  prlo%loghot  = -5.0
  prlo%sigma   = 20.0
  prlo%sigma2  = 20.0
  prlo%velz    = -1E4
  prlo%velz2   = -1E4
  prlo%logemnorm = -5.0

  !priors (high)
  prhi%age  = 15.0
  prhi%feh  = 0.6
  prhi%afe  = 0.6
  prhi%nhe  = 0.6
  prhi%cfe  = 0.6
  prhi%nfe  = 1.0
  prhi%nafe = 1.3
  prhi%mgfe = 0.6
  prhi%sife = 0.6
  prhi%kfe  = 0.6
  prhi%cafe = 0.6
  prhi%tife = 0.6
  prhi%vfe  = 0.6
  prhi%crfe = 0.6
  prhi%mnfe = 0.6
  prhi%cofe = 0.6
  prhi%nife = 0.6
  prhi%cufe = 0.6
  prhi%rbfe = 0.6
  prhi%srfe = 0.6
  prhi%yfe  = 0.6
  prhi%zrfe = 0.6
  prhi%bafe = 0.6
  prhi%eufe = 0.6
  prhi%teff = 100.0
  prhi%imf1 = 3.5
  prhi%imf2 = 3.5
  prhi%logfy   = -0.1
  prhi%logm7g  = -0.1
  prhi%hotteff = 30.0
  prhi%loghot  = -0.1
  prhi%sigma   = 1E4
  prhi%sigma2  = 1E4
  prhi%velz    = 1E4
  prhi%velz2   = 1E4
  prhi%logemnorm = 1.0

END SUBROUTINE SETUP_PARAMS
