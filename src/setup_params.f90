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
  prlo%feh  = -1.0
  prlo%afe  = -1.0
  prlo%nhe  = -1.0
  prlo%cfe  = -1.0
  prlo%nfe  = -1.0
  prlo%nafe = -1.0
  prlo%mgfe = -1.0
  prlo%sife = -1.0
  prlo%kfe  = -1.0
  prlo%cafe = -1.0
  prlo%tife = -1.0
  prlo%vfe  = -1.0
  prlo%crfe = -1.0
  prlo%mnfe = -1.0
  prlo%cofe = -1.0
  prlo%nife = -1.0
  prlo%cufe = -1.0
  prlo%rbfe = -1.0
  prlo%srfe = -1.0
  prlo%yfe  = -1.0
  prlo%zrfe = -1.0
  prlo%bafe = -1.0
  prlo%eufe = -1.0
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
  prhi%age  = 20.0
  prhi%feh  = 1.0
  prhi%afe  = 1.0
  prhi%nhe  = 1.0
  prhi%cfe  = 1.0
  prhi%nfe  = 1.0
  prhi%nafe = 1.0
  prhi%mgfe = 1.0
  prhi%sife = 1.0
  prhi%kfe  = 1.0
  prhi%cafe = 1.0
  prhi%tife = 1.0
  prhi%vfe  = 1.0
  prhi%crfe = 1.0
  prhi%mnfe = 1.0
  prhi%cofe = 1.0
  prhi%nife = 1.0
  prhi%cufe = 1.0
  prhi%rbfe = 1.0
  prhi%srfe = 1.0
  prhi%yfe  = 1.0
  prhi%zrfe = 1.0
  prhi%bafe = 1.0
  prhi%eufe = 1.0
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
