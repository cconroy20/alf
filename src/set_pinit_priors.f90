SUBROUTINE SET_PINIT_PRIORS(pos,prlo,prhi,velz)

  !define the first position (pos), and the lower and upper bounds 
  !on the priors (prlo, prhi).  The priors are defined in such a way
  !that if the user defines a prior limit that is **different from
  !the default parameter set**, then that value overrides the defaults below

  USE alf_vars; USE alf_utils, ONLY : myran
  IMPLICIT NONE

  TYPE(PARAMS), INTENT(inout) :: pos,prlo,prhi
  TYPE(PARAMS) :: test
  REAL(DP), OPTIONAL :: velz
  
  !---------------------------------------------------------------!
  !---------------------------------------------------------------!

  !setup the first position
  pos%logage    = myran()*0.4+0.6
  pos%zh        = myran()*1.0-1.0
  pos%feh       = myran()*0.4-0.2
  pos%ah        = myran()*0.4-0.2
  pos%nhe       = myran()*0.4-0.2
  pos%ch        = myran()*0.4-0.2
  pos%nh        = myran()*0.4-0.2
  pos%nah       = myran()*0.4-0.2
  pos%mgh       = myran()*0.4-0.2
  pos%sih       = myran()*0.4-0.2
  pos%kh        = myran()*0.4-0.2
  pos%cah       = myran()*0.4-0.2
  pos%tih       = myran()*0.4-0.2
  pos%vh        = myran()*0.4-0.2
  pos%crh       = myran()*0.4-0.2
  pos%mnh       = myran()*0.4-0.2
  pos%coh       = myran()*0.4-0.2
  pos%nih       = myran()*0.4-0.2
  pos%cuh       = myran()*0.4-0.2
  pos%srh       = myran()*0.4-0.2
  pos%bah       = myran()*0.4-0.2
  pos%euh       = myran()*0.4-0.2
  pos%teff      = myran()*80-40
  pos%imf1      = myran()*0.6-0.3 + 1.3
  IF (fit_2ximf.EQ.1) THEN
     pos%imf2        = myran()*0.6-0.3 + 2.3
  ELSE
     pos%imf2        = myran()*0.1 + 0.1
  ENDIF
  pos%logfy          = myran()*1-4
  pos%fy_logage      = myran()*0.3
  pos%logm7g         = myran()*1-4
  pos%hotteff        = myran()*5+15
  pos%loghot         = myran()*1-4
  pos%chi2           = huge_number
  pos%sigma          = myran()*100+100
  pos%sigma2         = myran()*100+100
  pos%velz2          = myran()*10-5
  pos%logtrans       = myran()*3-4.
  pos%logemline_h    = myran()*1-4
  pos%logemline_oiii = myran()*1-4
  pos%logemline_ni   = myran()*1-4
  pos%logemline_nii  = myran()*1-4
  pos%logemline_sii  = myran()*1-4
  pos%jitter         = myran()*0.5+0.75

  IF (PRESENT(velz)) THEN
     pos%velz  = velz + (myran()*10-5)
  ELSE
     pos%velz  = myran()*10-5
  ENDIF

  !these pr=test statements allow the user to pre-set
  !specific priors at the beginning of alf.f90; those
  !choices are then not overwritten below

  !priors (low)
  IF (fit_type.EQ.0) THEN
     !in this case we're fitting a two component model
     !so dont allow them to overlap in age
     IF (prlo%logage.EQ.test%logage) prlo%logage = LOG10(3.0)
  ELSE
     !in this case we have a single age model, so it needs to
     !cover the full range
     IF (prlo%logage.EQ.test%logage) prlo%logage = LOG10(0.5)
  ENDIF
  IF (prlo%zh.EQ.test%zh) prlo%zh             = -1.8
  IF (prlo%feh.EQ.test%feh) prlo%feh          = -0.3
  IF (prlo%ah.EQ.test%ah) prlo%ah             = -0.3
  IF (prlo%nhe.EQ.test%nhe) prlo%nhe          = -0.3
  IF (prlo%ch.EQ.test%ch) prlo%ch             = -0.3
  IF (prlo%nh.EQ.test%nh) prlo%nh             = -0.3
  IF (prlo%nah.EQ.test%nah) prlo%nah          = -0.3
  IF (prlo%mgh.EQ.test%mgh) prlo%mgh          = -0.3
  IF (prlo%sih.EQ.test%sih) prlo%sih          = -0.3
  IF (prlo%kh.EQ.test%kh) prlo%kh             = -0.3
  IF (prlo%cah.EQ.test%cah) prlo%cah          = -0.3
  IF (prlo%tih.EQ.test%tih) prlo%tih          = -0.3
  IF (prlo%vh.EQ.test%vh) prlo%vh             = -0.3
  IF (prlo%crh.EQ.test%crh) prlo%crh          = -0.3
  IF (prlo%mnh.EQ.test%mnh) prlo%mnh          = -0.3
  IF (prlo%coh.EQ.test%coh) prlo%coh          = -0.3
  IF (prlo%nih.EQ.test%nih) prlo%nih          = -0.3
  IF (prlo%cuh.EQ.test%cuh) prlo%cuh          = -0.3
  IF (prlo%srh.EQ.test%srh) prlo%srh          = -0.3
  IF (prlo%bah.EQ.test%bah) prlo%bah          = -0.6
  IF (prlo%euh.EQ.test%euh) prlo%euh          = -0.6
  IF (prlo%teff.EQ.test%teff) prlo%teff       = -75.0
  IF (prlo%imf1.EQ.test%imf1) prlo%imf1       = 0.5
  IF (fit_2ximf.EQ.1) THEN
     IF (prlo%imf2.EQ.test%imf2) prlo%imf2    = 0.5
  ELSE
     IF (prlo%imf2.EQ.test%imf2) prlo%imf2    = 0.08
  ENDIF
  IF (prlo%logfy.EQ.test%logfy) prlo%logfy    = -6.0
  IF (prlo%fy_logage.EQ.test%fy_logage) prlo%fy_logage = LOG10(0.5)
  IF (prlo%logm7g.EQ.test%logm7g) prlo%logm7g   = -6.0
  IF (prlo%hotteff.EQ.test%hotteff) prlo%hotteff= 8.0
  IF (prlo%loghot.EQ.test%loghot) prlo%loghot   = -6.0
  IF (prlo%sigma.EQ.test%sigma) prlo%sigma      = 10.0
  IF (prlo%sigma2.EQ.test%sigma2) prlo%sigma2   = 10.0
  IF (prlo%velz.EQ.test%velz) prlo%velz         = -1E3
  IF (prlo%velz2.EQ.test%velz2) prlo%velz2      = -1E3
  IF (prlo%logtrans.EQ.test%logtrans) prlo%logtrans = -6.0
  IF (prlo%logemline_h.EQ.test%logemline_h) prlo%logemline_h          = -6.0
  IF (prlo%logemline_oiii.EQ.test%logemline_oiii) prlo%logemline_oiii = -6.0
  IF (prlo%logemline_sii.EQ.test%logemline_sii) prlo%logemline_sii    = -6.0
  IF (prlo%logemline_ni.EQ.test%logemline_ni) prlo%logemline_ni       = -6.0
  IF (prlo%logemline_nii.EQ.test%logemline_nii) prlo%logemline_nii    = -6.0
  IF (prlo%jitter.EQ.test%jitter) prlo%jitter    = 0.1

 
  !priors (high)
  !if you change the prior on the age, also change the max 
  !age allowed in getmodel
  IF (prhi%logage.EQ.test%logage) prhi%logage = LOG10(15.0)
  IF (prhi%zh.EQ.test%zh) prhi%zh             = 0.3
  IF (prhi%feh.EQ.test%feh) prhi%feh          = 0.5
  IF (prhi%ah.EQ.test%ah) prhi%ah             = 0.5
  IF (prhi%nhe.EQ.test%nhe) prhi%nhe          = 0.5
  IF (prhi%ch.EQ.test%ch) prhi%ch             = 0.5
  IF (prhi%nh.EQ.test%nh) prhi%nh             = 1.0
  IF (prhi%nah.EQ.test%nah) prhi%nah          = 1.0
  IF (prhi%mgh.EQ.test%mgh) prhi%mgh          = 0.5
  IF (prhi%sih.EQ.test%sih) prhi%sih          = 0.5
  IF (prhi%kh.EQ.test%kh) prhi%kh             = 0.5
  IF (prhi%cah.EQ.test%cah) prhi%cah          = 0.5
  IF (prhi%tih.EQ.test%tih) prhi%tih          = 0.5
  IF (prhi%vh.EQ.test%vh) prhi%vh             = 0.5
  IF (prhi%crh.EQ.test%crh) prhi%crh          = 0.5
  IF (prhi%mnh.EQ.test%mnh) prhi%mnh          = 0.5
  IF (prhi%coh.EQ.test%coh) prhi%coh          = 0.5
  IF (prhi%nih.EQ.test%nih) prhi%nih          = 0.5
  IF (prhi%cuh.EQ.test%cuh) prhi%cuh          = 0.5
  IF (prhi%srh.EQ.test%srh) prhi%srh          = 0.5
  IF (prhi%bah.EQ.test%bah) prhi%bah          = 0.5
  IF (prhi%euh.EQ.test%euh) prhi%euh          = 0.5
  IF (prhi%teff.EQ.test%teff) prhi%teff       = 75.0
  IF (prhi%imf1.EQ.test%imf1) prhi%imf1       = 3.5
  IF (fit_2ximf.EQ.1) THEN
     IF (prhi%imf2.EQ.test%imf2) prhi%imf2       = 3.5
  ELSE
     IF (prhi%imf2.EQ.test%imf2) prhi%imf2       = 0.5
  ENDIF
  IF (prhi%logfy.EQ.test%logfy) prhi%logfy    = -0.7
  IF (prhi%fy_logage.EQ.test%fy_logage) prhi%fy_logage = LOG10(3.0)
  IF (prhi%logm7g.EQ.test%logm7g) prhi%logm7g   = -1.0
  IF (prhi%hotteff.EQ.test%hotteff) prhi%hotteff= 30.0
  IF (prhi%loghot.EQ.test%loghot) prhi%loghot   = -1.0
  IF (prhi%sigma.EQ.test%sigma) prhi%sigma      = 1E3
  IF (prhi%sigma2.EQ.test%sigma2) prhi%sigma2   = 1E3
  IF (prhi%velz.EQ.test%velz) prhi%velz         = 2E4
  IF (prhi%velz2.EQ.test%velz2) prhi%velz2      = 1E3
  IF (prhi%logtrans.EQ.test%logtrans) prhi%logtrans = 1.0
  IF (prhi%logemline_h.EQ.test%logemline_h) prhi%logemline_h          = 1.0
  IF (prhi%logemline_oiii.EQ.test%logemline_oiii) prhi%logemline_oiii = 1.0
  IF (prhi%logemline_sii.EQ.test%logemline_sii) prhi%logemline_sii    = 1.0
  IF (prhi%logemline_ni.EQ.test%logemline_ni) prhi%logemline_ni       = 1.0
  IF (prhi%logemline_nii.EQ.test%logemline_nii) prhi%logemline_nii    = 1.0
  IF (prhi%jitter.EQ.test%jitter) prhi%jitter    = 10.0


  !reset the initial parameters if the priors have been altered
  !this should be done for every parameter to ensure that the altered parameters
  !do not fall outside of the prior range.  In practice only these params are 
  !frequently altered.
  IF (prhi%logtrans.NE.test%logtrans) &
       pos%logtrans = myran()*(prhi%logtrans-prlo%logtrans)+prlo%logtrans
  IF (prhi%logfy.NE.test%logfy) &
       pos%logfy = myran()*(prhi%logfy-prlo%logfy)+prlo%logfy
  IF (prhi%loghot.NE.test%loghot) &
       pos%loghot = myran()*(prhi%loghot-prlo%loghot)+prlo%loghot
  IF (prhi%logm7g.NE.test%logm7g) &
       pos%logm7g = myran()*(prhi%logm7g-prlo%logm7g)+prlo%logm7g
  IF (prhi%logemline_h.NE.test%logemline_h) &
       pos%logemline_h = myran()*(prhi%logemline_h-prlo%logemline_h)+prlo%logemline_h
  IF (prhi%logemline_oiii.NE.test%logemline_oiii) &
       pos%logemline_oiii = myran()*(prhi%logemline_oiii-prlo%logemline_oiii)+prlo%logemline_oiii
  IF (prhi%logemline_sii.NE.test%logemline_sii) &
       pos%logemline_sii = myran()*(prhi%logemline_sii-prlo%logemline_sii)+prlo%logemline_sii
  IF (prhi%logemline_ni.NE.test%logemline_ni) &
       pos%logemline_ni = myran()*(prhi%logemline_ni-prlo%logemline_ni)+prlo%logemline_ni
  IF (prhi%logemline_nii.NE.test%logemline_nii) &
       pos%logemline_nii = myran()*(prhi%logemline_nii-prlo%logemline_nii)+prlo%logemline_nii

  IF (prhi%zh.NE.test%zh.OR.prlo%zh.NE.test%zh) &
       pos%zh = myran()*(prhi%zh-prlo%zh)+prlo%zh

 IF (prhi%teff.NE.test%teff.OR.prlo%teff.NE.test%teff) &
       pos%teff = myran()*(prhi%teff-prlo%teff)+prlo%teff

END SUBROUTINE SET_PINIT_PRIORS
