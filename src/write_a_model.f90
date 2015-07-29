PROGRAM WRITE_A_MODEL

  !write a model to file

  USE alf_vars; USE alf_utils
  USE nr, ONLY : gasdev,locate,powell,ran1
  USE ran_state, ONLY : ran_seed,ran_init

  IMPLICIT NONE

  INTEGER  :: i
  REAL(DP) :: s2n,lmin,lmax,ires=0.
  REAL(DP), DIMENSION(nl) :: mspec,lam,err,gdev
  CHARACTER(100)  :: file=''
  TYPE(PARAMS)   :: pos
  
  !-----------------------------------------------------------!
  !-----------------------------------------------------------!

  !initialize the random number generator
  CALL INIT_RANDOM_SEED()
  !compute an array of gaussian deviates
  CALL GASDEV(gdev)

  file = 'age+10.0_sn+100_sigma200_Ha-6.spec'
  s2n  = 1000.0
  lmin = 3800.
  lmax = 10000.
  pos%sigma  = 200.0
  pos%logage = 1.0

  ires = 10.

  !pos%sigma   = 10.0
  !pos%logage  = LOG10(8.0)
  pos%feh     = 0.0
  pos%ah      = 0.0
  pos%nhe     = 0.0
  pos%ch      = 0.0
  pos%nh      = 0.0
  pos%nah     = 0.0
  pos%mgh     = 0.0
  pos%sih     = 0.0
  pos%kh      = 0.0
  pos%cah     = 0.0
  pos%tih     = 0.0
  pos%vh      = 0.0
  pos%crh     = 0.0
  pos%mnh     = 0.0
  pos%coh     = 0.0
  pos%nih     = 0.0
  pos%cuh     = 0.0
  pos%rbh     = 0.0
  pos%srh     = 0.0
  pos%yh      = 0.0
  pos%zrh     = 0.0
  pos%bah     = 0.0
  pos%euh     = 0.0
  pos%teff    = 0.0
  pos%imf1    = 1.3
  pos%imf2    = 2.3
  pos%logfy   = -5.0
  pos%sigma2  = 300.
  !pos%velz    = 0.0
  pos%velz2   = 0.0
  pos%logm7g  = -5.0
  pos%hotteff = 20.0
  pos%loghot  = -5.0
  !pos%logemnorm = -10.0

  !force a constant instrumental resolution
  datmax=10000
  DO i=1,datmax
     data(i)%lam=i+3500
  ENDDO
  data(1:datmax)%ires = ires

  !read in the SSPs and bandpass filters
  CALL SETUP()
  lam = sspgrid%lam

  !define the log wavelength grid used in velbroad.f90
  nl_fit = MIN(MAX(locate(lam,lmax+500.0),1),nl)
  dlstep = (LOG(sspgrid%lam(nl_fit))-LOG(sspgrid%lam(1)))/nl_fit
  DO i=1,nl_fit
     lnlam(i) = i*dlstep+LOG(sspgrid%lam(1))
  ENDDO
  l1(1) = lmin
  nlint = 2
  l2(nlint) = lmax

  !get a model spectrum
  CALL GETMODEL(pos,mspec)

  IF (1.EQ.1) THEN
     err   = mspec/s2n
     mspec = mspec + err*gdev
  ELSE
     DO i=1,nl
        IF (lam(i).LT.4600) THEN
           err(i) = mspec(i)/10. 
           mspec(i) = mspec(i) + err(i)*gdev(i)
        ELSE
           err(i) = mspec(i)/30.
           mspec(i) = mspec(i) + err(i)*gdev(i)
        ENDIF
     ENDDO
  ENDIF
  
  !write model spectrum to file
  OPEN(12,FILE=TRIM(SPECFIT_HOME)//'models/'//&
       TRIM(file),STATUS='REPLACE')
  DO i=1,nl
     IF (lam(i).GE.lmin.AND.lam(i).LE.lmax) THEN
        WRITE(12,'(F10.3,2ES12.4,2x,F3.1)') lam(i),mspec(i),err(i),1.0
     ENDIF
  ENDDO
  CLOSE(12)


END PROGRAM WRITE_A_MODEL
