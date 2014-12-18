PROGRAM ALF

  !  Master program to fit the absorption line spectrum
  !  of an old (>~1 Gyr), metal-rich ([Fe/H]>~-0.3) stellar
  !  population to the CvD model SSPs.

  ! Some things to keep in mind:
  ! 1. The prior bounds on the parameters are specified in setup_params. 
  !    Always make sure that the output parameters are not hitting a prior.
  ! 2. Make sure that the chain is converged in all relevant parameters
  !    by plotting the chain trace (parameter vs. chain step).
  ! 3. Never ever ever use this code blindly.  Fitting spectra is a 
  !    subtle art and the code can easily fool you if you don't know
  !    what you're doing.  Make sure you understand *why* the code is 
  !    settling on a particular parameter value.  
  
  !---------------------------------------------------------------!
  !---------------------------------------------------------------!

  USE sfvars; USE nr, ONLY : ran,locate,powell,ran1
  USE ran_state, ONLY : ran_seed,ran_init; USE sfutils
  IMPLICIT NONE

  !number of chain steps to print to file
  INTEGER, PARAMETER :: nmcmc=10000
  !length of chain burn-in
  INTEGER, PARAMETER :: nburn=1000000
  !start w/ powell minimization?
  INTEGER, PARAMETER :: dopowell=1
  !number of walkers for emcee
  INTEGER, PARAMETER :: nwalkers=100
  !Powell iteration tolerance
  REAL(DP), PARAMETER :: ftol=0.1

  INTEGER  :: i,j,totacc=0,iter=30
  REAL(DP) :: velz,msto,minchi2=huge_number,fret,wdth,bret=huge_number
  REAL(DP), DIMENSION(nl)   :: mspec=0.0,mspecmw=0.0,lam=0.0
  REAL(DP), DIMENSION(nfil) :: m2l=0.0,m2lmw=0.0
  REAL(DP), DIMENSION(npar) :: oposarr=0.,bposarr=0.0
  REAL(DP), DIMENSION(3,npar+2*nfil) :: runtot=0.0
  REAL(DP), DIMENSION(npar,npar)     :: xi=0.0
  CHARACTER(10) :: time=''
  CHARACTER(50) :: file='',tag=''
  TYPE(PARAMS)  :: opos,prlo,prhi,bpos
  !the next three definitions are for emcee
  REAL(DP), DIMENSION(npar,nwalkers) :: pos_emcee
  REAL(DP), DIMENSION(nwalkers)      :: lp_emcee
  INTEGER,  DIMENSION(nwalkers)      :: accept_emcee
  
  !---------------------------------------------------------------!
  !---------------------------Setup-------------------------------!
  !---------------------------------------------------------------!

  !flag determining the level of complexity
  !0=full, 1=simple, 2=super-simple.  See sfvars for details
  fitsimple = 0
  IF (fitsimple.EQ.1.OR.fitsimple.EQ.2) mwimf=1

  !prhi%logm7g = -2.0

  !initialize the random number generator
  CALL INIT_RANDOM_SEED()

  IF (IARGC().LT.1) THEN
     WRITE(*,*) 'ALF ERROR: You need to specify an input file'
     STOP
  ELSE
     CALL GETARG(1,file)
  ENDIF
  IF (IARGC().GT.1) THEN
     tag(1:1)='_'
     CALL GETARG(2,tag(2:))
  ENDIF

  !write some important variables to screen
  WRITE(*,*) 
  WRITE(*,'(" ************************************")') 
  WRITE(*,'("   dopowell  =",I2)') dopowell
  WRITE(*,'("  fitsimple  =",I2)') fitsimple
  WRITE(*,'("      mwimf  =",I2)') mwimf
  WRITE(*,'("  force_nah  =",I2)') force_nah
  WRITE(*,'("  age-dep Rf =",I2)') use_age_dep_resp_fcns
  WRITE(*,'("  Nburn      = ",I8)') nburn
  WRITE(*,'("  Nchain     = ",I8)') nmcmc
  WRITE(*,'("  filename   = ",A)') TRIM(file)//TRIM(tag)
  WRITE(*,'(" ************************************")') 
 
  CALL date_and_time(TIME=time)
  WRITE(*,*) 
  WRITE(*,*) 'Start Time '//time(1:2)//':'//time(3:4)

  !read in the SSPs and bandpass filters
  CALL SFSETUP()
  lam = sspgrid%lam

  !read in the data and wavelength boundaries
  CALL READ_DATA(file)
  WRITE(*,'("  Using ",I1," wavelength intervals")') nlint
  IF (l2(nlint).GT.lam(nl).OR.l1(1).LT.lam(1)) THEN
     WRITE(*,*) 'ERROR: wavelength boundaries exceed model wavelength grid'
     WRITE(*,'(4F8.1)') l2(nlint),lam(nl),l1(1),lam(1)
     STOP
  ENDIF

  !define the log wavelength grid used in velbroad.f90
  nl_fit = MIN(MAX(locate(lam,l2(nlint)+500.0),1),nl)
  !nl_fit = MIN(MAX(locate(lam,15000.d0),1),nl)
  dlstep = (LOG(sspgrid%lam(nl_fit))-LOG(sspgrid%lam(1)))/nl_fit
  DO i=1,nl_fit
     lnlam(i) = i*dlstep+LOG(sspgrid%lam(1))
  ENDDO

  !masked regions have wgt=0.0.  We'll use wgt as a pseudo-error
  !array in contnormspec, so turn these into large numbers
  data%wgt = MIN(1/(data%wgt+tiny_number),huge_number)
  !fold the masked regions into the errors
  data%err = data%err * data%wgt

  !set initial params, step sizes, and prior ranges
  CALL SETUP_PARAMS(opos,prlo,prhi)

  !make an initial estimate of the redshift
  !we do this to help Powell minimization
  WRITE(*,*) ' Finding redshift...'
  IF (file(1:4).EQ.'cdfs') THEN
     velz = 0.0 
  ELSE 
     velz = getvelz()
  ENDIF
  opos%velz = velz
  WRITE(*,'("    best velocity: ",F6.1)') velz
  

  !convert the structures into their equivalent arrays
  CALL STR2ARR(1,prlo,prloarr)   !str->arr
  CALL STR2ARR(1,prhi,prhiarr)   !str->arr
  CALL STR2ARR(1,opos,oposarr)   !str->arr

  !---------------------------------------------------------------!
  !---------------------Powell minimization-----------------------!
  !---------------------------------------------------------------!

  IF (dopowell.EQ.1) THEN 

     WRITE(*,*) ' Running Powell...'
     powell_fitting = 1
     DO j=1,10
        xi=0.0
        DO i=1,npar
           xi(i,i) = 1E-2
        ENDDO
        fret = huge_number
        CALL SETUP_PARAMS(opos,prlo,prhi,velz=velz)
        CALL STR2ARR(1,opos,oposarr) !str->arr
        CALL POWELL(oposarr(1:npowell),xi(1:npowell,1:npowell),&
             ftol,iter,fret)
        CALL STR2ARR(2,opos,oposarr) !arr->str
        IF (fret.LT.bret) THEN
           bposarr = oposarr
           bpos    = opos
           bret    = fret
        ENDIF
     ENDDO
     powell_fitting = 0

     !use the best-fit Powell position for the first MCMC position
     CALL STR2ARR(2,opos,bposarr) !arr->str
 
     WRITE(*,'("    best velocity: ",F6.1)') opos%velz
     WRITE(*,'("    best sigma:    ",F6.1)') opos%sigma
     WRITE(*,'("    best age:      ",F6.1)') 10**opos%logage

  ENDIF

  IF (maskem.EQ.1) THEN
     !now that we have a good guess of the redshift and velocity dispersion, 
     !mask out regions where emission line contamination may be a problem
     !In full mode, the default is to actually *fit* for emissions lines.
     CALL MASKEMLINES(opos%velz,opos%sigma)
  ENDIF


  !---------------------------------------------------------------!
  !-----------------------------MCMC------------------------------!
  !---------------------------------------------------------------!

  WRITE(*,*) ' Running emcee...'
  CALL FLUSH()

  !open output file
  OPEN(12,FILE=TRIM(SPECFIT_HOME)//TRIM(OUTDIR)//&
       TRIM(file)//TRIM(tag)//'.mcmc',STATUS='REPLACE')
     
  !initialize the walkers
  DO j=1,nwalkers
     
     !random initialization of each walker
     CALL SETUP_PARAMS(opos,prlo,prhi,velz=velz)
     CALL STR2ARR(1,opos,pos_emcee(:,j))

     IF (dopowell.EQ.1) THEN
        !use the best-fit position from Powell, with small
        !random offsets to set up all the walkers, but only
        !do this for the params actually fit in Powell!
        !the first two params are velz and sigma so give them
        !larger variation.
        DO i=1,npowell
           IF (i.LE.2) wdth=10.0
           IF (i.GT.2) wdth=0.3
           pos_emcee(i,j) = bposarr(i) + wdth*(2.*myran()-1.0)
        ENDDO
     ENDIF

     !Compute the initial log-probability for each walker
     lp_emcee(j) = -0.5*func(pos_emcee(:, j))
     
  ENDDO
  
  !burn-in
  WRITE(*,*) '   burning in...'
  WRITE(*,'(A)',advance='no') '      Progress:'
  DO i=1,nburn/nwalkers
     CALL EMCEE_ADVANCE(npar,nwalkers,2.d0,pos_emcee,&
          lp_emcee,pos_emcee,lp_emcee,accept_emcee)
     IF (i.EQ.nburn/nwalkers/4.*1) THEN
        WRITE (*,'(A)',advance='no') ' ...25%'
        CALL FLUSH
     ENDIF
     IF (i.EQ.nburn/nwalkers/4.*2) WRITE (*,'(A)',advance='no') '...50%'
     IF (i.EQ.nburn/nwalkers/4.*3) WRITE (*,'(A)',advance='no') '...75%'
     IF (i.EQ.nburn/nwalkers/4.*4) WRITE (*,'(A)') '...100%'
  ENDDO
  
  !Run a production chain
  WRITE(*,*) '   production run...'
  DO i=1,nmcmc/nwalkers
     
     CALL EMCEE_ADVANCE(npar,nwalkers,2.d0,pos_emcee,&
          lp_emcee,pos_emcee,lp_emcee,accept_emcee)
     totacc = totacc + SUM(accept_emcee)
     
     DO j=1,nwalkers

        CALL STR2ARR(2,opos,pos_emcee(:,j)) !arr->str

        !compute the main sequence turn-off mass
        msto = 10**( msto_fit0 + msto_fit1*opos%logage )
        msto = MIN(MAX(msto,0.90),3.)
        CALL GETMODEL(opos,mspecmw,mw=1)     !get spectra for MW IMF
        CALL GETM2L(msto,lam,mspecmw,opos,m2lmw,mw=1) !compute M/L_MW
        IF (fitsimple.EQ.0.AND.mwimf.EQ.0) THEN
           CALL GETMODEL(opos,mspec)
           CALL GETM2L(msto,lam,mspec,opos,m2l) ! compute M/L
        ELSE
           m2l = m2lmw
           !these parameters aren't actually being updated
           pos_emcee(nparsimp+1:,j) = 0.0 
        ENDIF

        !write the chain element to file
        WRITE(12,'(ES12.5,1x,99(F9.4,1x))') -2.0*lp_emcee(j),&
             pos_emcee(:, j),m2l,m2lmw

        !keep the model with the lowest chi2
        IF (-2.0*lp_emcee(j).LT.minchi2) THEN
           bposarr = pos_emcee(:, j)
           minchi2 = -2.0*lp_emcee(j)
        ENDIF

        CALL UPDATE_RUNTOT(runtot,pos_emcee(:,j),m2l,m2lmw)

     ENDDO
     CALL FLUSH(12)
     
  ENDDO
  
  !save the best position to the structure
  CALL STR2ARR(2,bpos,bposarr)
  bpos%chi2 = minchi2
  
  CLOSE(12)

  CALL date_and_time(TIME=time)
  WRITE(*,*) 'End Time   '//time(1:2)//':'//time(3:4)
  WRITE(*,*) 
  WRITE(*,'("  Facc: ",F5.2)') REAL(totacc)/REAL(nmcmc)

  !---------------------------------------------------------------!
  !--------------------Write results to file----------------------!
  !---------------------------------------------------------------!

  OPEN(13,FILE=TRIM(SPECFIT_HOME)//TRIM(OUTDIR)//&
       TRIM(file)//TRIM(tag)//'.bestspec',STATUS='REPLACE')
  CALL STR2ARR(1,bpos,bposarr)
  !NB: the model written to file has the lowest chi^2
  fret = func(bposarr,spec=mspec,funit=13)
  CLOSE(13)
 
  !write best-fit parameters
  !here, "best-fit" is the mean of the posterior distributions
  OPEN(14,FILE=TRIM(SPECFIT_HOME)//TRIM(OUTDIR)//&
       TRIM(file)//TRIM(tag)//'.bestp',STATUS='REPLACE')
  WRITE(14,'(ES12.5,1x,99(F9.4,1x))') bpos%chi2, runtot(2,:)/runtot(1,:)
  CLOSE(14)

  !write one sigma errors on parameters
  OPEN(15,FILE=TRIM(SPECFIT_HOME)//TRIM(OUTDIR)//&
       TRIM(file)//TRIM(tag)//'.errp',STATUS='REPLACE')
  WRITE(15,'(ES12.5,1x,99(F9.4,1x))') 0.0, &
       SQRT( runtot(3,:)/runtot(1,:) - runtot(2,:)**2/runtot(1,:)**2 )
  CLOSE(15)

  WRITE(*,*)

END PROGRAM ALF
