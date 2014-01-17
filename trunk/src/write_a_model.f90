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

  file = 'age+8.0_mgfe+0.2_feh-0.05_nfe+0.2_cfe+0.2_afe+0.2_sigma+350_sn+10_v2.spec'
  s2n  = 10.0
  lmin = 3900.
  lmax = 5600.

  !read in the SSPs and bandpass filters
  CALL SFSETUP()
  lam = sspgrid%lam

  pos%sigma   = 350.0
  pos%age     = 8.0
  pos%feh     = -0.05
  pos%afe     = 0.2
  pos%nhe     = 0.0
  pos%cfe     = 0.2
  pos%nfe     = 0.2
  pos%nafe    = 0.0
  pos%mgfe    = 0.2
  pos%sife    = 0.0
  pos%kfe     = 0.0
  pos%cafe    = 0.0
  pos%tife    = 0.0
  pos%vfe     = 0.0
  pos%crfe    = 0.0
  pos%mnfe    = 0.0
  pos%cofe    = 0.0
  pos%nife    = 0.0
  pos%cufe    = 0.0
  pos%rbfe    = 0.0
  pos%srfe    = 0.0
  pos%yfe     = 0.0
  pos%zrfe    = 0.0
  pos%bafe    = 0.0
  pos%eufe    = 0.0
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
