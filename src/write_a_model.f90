PROGRAM WRITE_A_MODEL

  !write a model to file

  USE sfvars; USE sfutils
  USE nr, ONLY : gasdev,ran,locate,powell,ran1
  USE ran_state, ONLY : ran_seed,ran_init

  IMPLICIT NONE

  INTEGER  :: i
  REAL(DP) :: s2n,lmin,lmax
  REAL(DP), DIMENSION(nl) :: mspec,mspec2,lam,err,gdev
  CHARACTER(100)  :: file=''
  TYPE(PARAMS)   :: pos
  
  !-----------------------------------------------------------!
  !-----------------------------------------------------------!

  !initialize the random number generator
  CALL init_random_seed()
  !compute an array of gaussian deviates
  CALL GASDEV(gdev)

  file = 'age+13.0_sn+100.spec'
  s2n  = 100.0
  lmin = 3900.
  lmax = 5600.

  !read in the SSPs and bandpass filters
  CALL SFSETUP()
  lam = sspgrid%lam

  pos%sigma   = 10.0
  pos%logage  = LOG10(13.0)
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
  pos%sigma2  = 0.0
  pos%velz    = 0.0
  pos%velz2   = 0.0
  pos%logm7g  = -5.0
  pos%hotteff = 20.0
  pos%loghot  = -5.0
  pos%logemnorm = -5.0

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
