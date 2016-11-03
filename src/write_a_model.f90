PROGRAM WRITE_A_MODEL

  !write a model to file

  USE alf_vars; USE alf_utils
  USE nr, ONLY : gasdev,locate,powell,ran1
  USE ran_state, ONLY : ran_seed,ran_init

  IMPLICIT NONE

  INTEGER  :: i
  REAL(DP) :: s2n,lmin,lmax,ires=0.,emnorm
  REAL(DP), DIMENSION(nl) :: mspec,lam,err,gdev
  CHARACTER(100) :: file=''
  CHARACTER(50)  :: infile=''
  TYPE(PARAMS)   :: pos
  
  !-----------------------------------------------------------!
  !-----------------------------------------------------------!

  !instrumental resolution (<10 -> no broadening)
  ires = 1. !100.

  imf_type = 1

  !initialize the random number generator
  CALL INIT_RANDOM_SEED()

  !compute an array of gaussian deviates
  CALL GASDEV(gdev)

  file = 'model2_3.1_1.4_0.48.dat'
  s2n  = 1e5
  lmin = 3800.
  lmax = 11000.

  pos%sigma  = 1. !300.
  pos%logage = 1.13 !LOG10(10.0)
  pos%mgh    = 0.3
  pos%nah    = 0.4856
  pos%zh     = 0.07
  emnorm     = -5.0
  !Kroupa
  pos%imf1   = 3.1
  pos%imf2   = 1.4
  pos%imf3   = 0.08

  IF (imf_type.EQ.4) THEN
     !x=3.3
     pos%imf1 = 2.64
     pos%imf2 = 1.68
     pos%imf3 = 0.872
     pos%imf4 = 0.369
     !bottom-heavy
  !   pos%imf1 = 2.2
  !   pos%imf2 = 2.2
  !   pos%imf3 = 0.573
  !   pos%imf4 = 0.256
     !Kroupa
   !  pos%imf1 = 1.121
   !  pos%imf2 = 0.894
   !  pos%imf3 = 0.573
   !  pos%imf4 = 0.256
     !Salpeter x=2.3
  !   pos%imf1 = 1.741
  !   pos%imf2 = 1.155
  !   pos%imf3 = 0.603
  !   pos%imf4 = 0.256
     !flat
   !  pos%imf1 = -0.222
   !  pos%imf2 =  0.000
   !  pos%imf3 =  0.000
   !  pos%imf4 =  0.000
   !  pos%imf5 =  0.000
  ENDIF

  pos%logemline_h    = emnorm
  pos%logemline_oiii = emnorm
  pos%logemline_sii  = emnorm
  pos%logemline_ni   = emnorm
  pos%logemline_nii  = emnorm

  !force a constant instrumental resolution
  !needs to be done this way for setup.f90 to work
  datmax=lmax-lmin
  DO i=1,datmax
     data(i)%lam=i+lmin
     data(i)%ires = ires
  ENDDO

  !use the LRIS instrumental dispersion
  !infile = 'snl2_lris'
  !CALL READ_DATA(infile)

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
  OPEN(12,FILE=TRIM(ALF_HOME)//'models/'//TRIM(file),STATUS='REPLACE')
  WRITE(12,'("# 0.400 0.470")') 
  WRITE(12,'("# 0.470 0.570")')
  WRITE(12,'("# 0.570 0.640")')
  WRITE(12,'("# 0.800 0.892")')
  WRITE(12,'("# 0.963 1.015")') 
  DO i=1,nl
     IF (lam(i).GE.lmin.AND.lam(i).LE.lmax) THEN
        WRITE(12,'(F10.3,2ES12.4,2x,F4.1,2x,F7.2)') &
             lam(i),mspec(i),err(i),1.0,ires
     ENDIF
  ENDDO
  CLOSE(12)


END PROGRAM WRITE_A_MODEL
