SUBROUTINE SET_PINIT_PRIORS(pos,prlo,prhi,velz)

  !define the first position (pos), and the lower and upper bounds 
  !on the priors (prlo, prhi).  The priors are defined in such a way
  !that if the user defines a prior limit that is **different from
  !the default parameter set**, then that value overrides the defaults below

  USE alf_vars; USE alf_utils, ONLY : myran,str2arr
  IMPLICIT NONE

  TYPE(PARAMS), INTENT(inout) :: pos,prlo,prhi
  TYPE(PARAMS) :: test,tprlo,tprhi
  INTEGER :: i
  REAL(DP) :: tmps
  REAL(DP), OPTIONAL :: velz
  REAL(DP), DIMENSION(npar) :: prloarr1=0.,prhiarr1=0.,tprloarr1=0.
  REAL(DP), DIMENSION(npar) :: testarr1=0.,posarr1=0.,tprhiarr1=0.

  !---------------------------------------------------------------!
  !---------------------------------------------------------------!

  tprlo = prlo
  tprhi = prhi
  

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
  pos%teff      = myran()*80.0-40.
  pos%logfy          = myran()*2-4
  pos%fy_logage      = myran()*0.3
  pos%logm7g         = myran()*1-4
  pos%hotteff        = myran()*5+15
  pos%loghot         = myran()*1-4
  pos%chi2           = huge_number
  pos%sigma          = myran()*100+100
  pos%sigma2         = myran()*100+100
  pos%velz2          = myran()*10-5
  pos%logtrans       = myran()*4-4.
  pos%logemline_h    = myran()*2-4
  pos%logemline_oiii = myran()*2-4
  pos%logemline_ni   = myran()*2-4
  pos%logemline_nii  = myran()*2-4
  pos%logemline_sii  = myran()*2-4
  pos%jitter         = myran()*0.5+0.75
  pos%logsky         = myran()*2-4

  IF (imf_type.LE.3) THEN
     pos%imf1      = myran()*1.0-0.3 + 1.3
     pos%imf3      = myran()*0.1 + 0.1
  ELSE
     pos%imf1      = myran()*1+0.5
     pos%imf3      = myran()*1
  ENDIF
  IF (imf_type.EQ.0.OR.imf_type.EQ.1.OR.imf_type.EQ.3) THEN
     pos%imf2        = myran()*1.5-0.75 + 2.0
  ELSE IF (imf_type.EQ.2) THEN
     pos%imf2        = myran()*0.1 + 0.1
  ELSE IF (imf_type.EQ.4) THEN
     pos%imf2      = myran()*0.5+0.5
  ENDIF
  pos%imf4         = myran()*0.5

  IF (PRESENT(velz)) THEN
     IF (ABS(pos%velz).LE.tiny_number) THEN
        pos%velz  = myran()*1E4-1E3
     ELSE
        pos%velz  = velz + (myran()*10-5)
     ENDIF
  ELSE
     pos%velz  = myran()*1E4-1E3
  ENDIF

  IF (nonpimf_alpha.EQ.2.3) THEN 
     pos%imf1 = myran()*1.0-0.5
     pos%imf2 = myran()*1.0-0.5
     pos%imf3 = myran()*1.0-0.5
     pos%imf4 = myran()*1.0-0.5
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
  IF (prlo%teff.EQ.test%teff) prlo%teff       = -50.0
  IF (prlo%logfy.EQ.test%logfy) prlo%logfy      = -6.0
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
  IF (prlo%logsky.EQ.test%logsky) prlo%logsky    = -6.0

  IF (imf_type.LE.3) THEN
     IF (prlo%imf1.EQ.test%imf1) prlo%imf1       = 0.5
     IF (prlo%imf3.EQ.test%imf3) prlo%imf3       = 0.08
  ELSE
     IF (prlo%imf1.EQ.test%imf1) prlo%imf1    = -3.0
     IF (prlo%imf3.EQ.test%imf3) prlo%imf3    = -3.0
  ENDIF
  IF (imf_type.EQ.0.OR.imf_type.EQ.1.OR.imf_type.EQ.3) THEN
     IF (prlo%imf2.EQ.test%imf2) prlo%imf2    = 0.5
  ELSE IF (imf_type.EQ.2) THEN
     IF (prlo%imf2.EQ.test%imf2) prlo%imf2    = 0.08
  ELSE IF (imf_type.EQ.4) THEN
     IF (prlo%imf2.EQ.test%imf2) prlo%imf2    = -3.0
  ENDIF
  IF (prlo%imf4.EQ.test%imf4) prlo%imf4    = -3.0


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
  IF (prhi%teff.EQ.test%teff) prhi%teff       = 50.0
  IF (prhi%logfy.EQ.test%logfy) prhi%logfy      = -0.1
  IF (prhi%fy_logage.EQ.test%fy_logage) prhi%fy_logage = LOG10(3.0)
  IF (prhi%logm7g.EQ.test%logm7g) prhi%logm7g   = -1.0
  IF (prhi%hotteff.EQ.test%hotteff) prhi%hotteff= 30.0
  IF (prhi%loghot.EQ.test%loghot) prhi%loghot   = -1.0
  IF (prhi%sigma.EQ.test%sigma) prhi%sigma      = 1E3
  IF (prhi%sigma2.EQ.test%sigma2) prhi%sigma2   = 1E3
  IF (prhi%velz.EQ.test%velz) prhi%velz         = 1E5
  IF (prhi%velz2.EQ.test%velz2) prhi%velz2      = 1E3
  IF (prhi%logtrans.EQ.test%logtrans) prhi%logtrans = 1.0
  IF (prhi%logemline_h.EQ.test%logemline_h) prhi%logemline_h          = 1.0
  IF (prhi%logemline_oiii.EQ.test%logemline_oiii) prhi%logemline_oiii = 1.0
  IF (prhi%logemline_sii.EQ.test%logemline_sii) prhi%logemline_sii    = 1.0
  IF (prhi%logemline_ni.EQ.test%logemline_ni) prhi%logemline_ni       = 1.0
  IF (prhi%logemline_nii.EQ.test%logemline_nii) prhi%logemline_nii    = 1.0
  IF (prhi%jitter.EQ.test%jitter) prhi%jitter    = 10.0
  IF (prhi%logsky.EQ.test%logsky) prhi%logsky    = 2.0

  IF (imf_type.LE.3) THEN
     IF (prhi%imf1.EQ.test%imf1) prhi%imf1       = 3.5
     IF (prhi%imf3.EQ.test%imf3) prhi%imf3       = 0.4
  ELSE
     IF (prhi%imf1.EQ.test%imf1) prhi%imf1 = 3.0
     IF (prhi%imf3.EQ.test%imf3) prhi%imf3 = 3.0
  ENDIF
  IF (imf_type.EQ.0.OR.imf_type.EQ.1.OR.imf_type.EQ.3) THEN
     IF (prhi%imf2.EQ.test%imf2) prhi%imf2    = 3.5
  ELSE IF (imf_type.EQ.2) THEN
     IF (prhi%imf2.EQ.test%imf2) prhi%imf2    = 0.5
  ELSE IF (imf_type.EQ.4) THEN
     IF (prhi%imf2.EQ.test%imf2) prhi%imf2    = 3.0
  ENDIF
  IF (prhi%imf4.EQ.test%imf4) prhi%imf4 = 3.0

  !--------------------------------------------------------------------------!
  !-------reset the initial parameters if the priors have been altered-------!
  !--------------------------------------------------------------------------!

  !str->arr
  CALL STR2ARR(1,tprlo,tprloarr1)   
  CALL STR2ARR(1,tprhi,tprhiarr1)
  CALL STR2ARR(1,prlo,prloarr1)   
  CALL STR2ARR(1,prhi,prhiarr1)
  CALL STR2ARR(1,test,testarr1)
  CALL STR2ARR(1,pos,posarr1)

  !test the priors and if the priors have been altered then 
  !re-initialize the parameters within the prior range
  !NB: why not simply always initialize the starting position this way?
  DO i=3,npar
     IF (prhiarr1(i).LE.prloarr1(i)) THEN
        WRITE(*,*) 'SET_PINIT_PRIORS ERROR: prhi < prlo!', i
        STOP
     ENDIF
     !reset the initial parameters randomly within the prior range
     !if both the initial priors have changed *and* the position is outside
     !of the new priors
     IF (tprhiarr1(i).NE.testarr1(i).OR.tprloarr1(i).NE.testarr1(i)) THEN
        IF (posarr1(i).GE.prhiarr1(i).OR.posarr1(i).LE.prloarr1(i)) &
             posarr1(i) = myran()*(prhiarr1(i)-prloarr1(i)) + prloarr1(i)
     ENDIF
  ENDDO

  !arr->str
  CALL STR2ARR(2,pos,posarr1)


END SUBROUTINE SET_PINIT_PRIORS
