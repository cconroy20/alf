PROGRAM WRITE_A_MODEL

  !write a model to file

  USE alf_vars; USE alf_utils
  USE nr, ONLY : gasdev,locate,powell,ran1
  USE ran_state, ONLY : ran_seed,ran_init

  IMPLICIT NONE

  INTEGER  :: i,j
  REAL(DP) :: s2n,s2np,lmin,lmax,ires=0.,emnorm,msto
  REAL(DP), DIMENSION(nl) :: gspec,mspec,lam,err,gdev
  REAL(DP), DIMENSION(nfil) :: m2l=0.0
  CHARACTER(100) :: file=''
  CHARACTER(50)  :: infile=''
  CHARACTER(2), DIMENSION(10)  :: str=''
  CHARACTER(2) :: is
  TYPE(PARAMS) :: pos
  
  !-----------------------------------------------------------!
  !-----------------------------------------------------------!

  !instrumental resolution (<10 -> no broadening)
  ires     = 1. !100.
  imf_type = 1
  fit_type = 1
  
  str = (/'00','01','02','03','04','05','06','07','08','09'/)

  !initialize the random number generator
  CALL INIT_RANDOM_SEED()

  lmin = 3700.
  lmax = 11000.
  !force a constant instrumental resolution
  !needs to be done this way for setup.f90 to work
  datmax=lmax-lmin
  DO i=1,datmax
     data(i)%lam  = i+lmin
     data(i)%ires = ires
  ENDDO

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

  !loop to generate multiple mock datasets
  DO j=0,0

     !compute an array of gaussian deviates
     CALL GASDEV(gdev)

     !string for indexing of filenames
     IF (j.LE.10) THEN
        is = str(j+1)
     ELSE
        WRITE(is,'(I2)') j
     ENDIF

     file = 'hmodel_t10.0_zp0.0.dat'
     pos%sigma  = 1.

     s2n  = 10000.  !S/N per A
 
     pos%logage = LOG10(10.0)
     pos%zh     = 0.0
     emnorm     = -5.0

     pos%imf1   = 1.3
     pos%imf2   = 2.3
     pos%imf3   = 0.08

     pos%logemline_h    = emnorm
     pos%logemline_oiii = emnorm
     pos%logemline_sii  = emnorm
     pos%logemline_ni   = emnorm
     pos%logemline_nii  = emnorm

     !get a model spectrum
     gspec = 0.0
     mspec = 0.0
     CALL GETMODEL(pos,mspec)

     !compute M/L
     msto = MAX(MIN(10**(msto_t0+msto_t1*pos%logage) * &
          (msto_z0+msto_z1*pos%zh+msto_z2*pos%zh**2),3.0),0.75)
     CALL GETM2L(msto,lam,mspec,pos,m2l)
     !print to screen
     !write(*,*) m2l
     
     DO i=1,nl
        !this is to convert between S/N per A and S/N per pix
        IF (lam(i).LT.7500.) THEN
           s2np = s2n*SQRT(0.9)
        ELSE
           s2np = s2n*SQRT(2.5)
        ENDIF
        err(i)   = mspec(i)/s2np
        gspec(i) = mspec(i) + err(i)*gdev(i)
     ENDDO
 
     !write model spectrum to file
     OPEN(12,FILE=TRIM(ALF_HOME)//'models/'//TRIM(file),STATUS='REPLACE')
     WRITE(12,'("# 0.400 0.470")') 
     WRITE(12,'("# 0.470 0.560")')
     !WRITE(12,'("# 0.570 0.640")')
     !WRITE(12,'("# 0.640 0.800")')
     !WRITE(12,'("# 0.800 0.892")')
     !WRITE(12,'("# 0.963 1.015")') 
     DO i=1,nl
        IF (lam(i).GE.lmin.AND.lam(i).LE.lmax) THEN
           WRITE(12,'(F10.3,2ES12.4,2x,F4.1,2x,F7.2)') &
                lam(i),gspec(i),err(i),1.0,ires
        ENDIF
     ENDDO
     CLOSE(12)

     
  ENDDO


END PROGRAM WRITE_A_MODEL
