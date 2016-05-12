PROGRAM ALF

  !  Master program to fit the absorption line spectrum
  !  of an old (>~1 Gyr) stellar population

  ! Some things to keep in mind:
  ! 1. The prior bounds on the parameters are specified in set_pinit_priors. 
  !    Always make sure that the output parameters are not hitting a prior.
  ! 2. Make sure that the chain is converged in all relevant parameters
  !    by plotting the chain trace (parameter vs. chain step).
  ! 3. Never ever ever use this code blindly.  Fitting spectra is a 
  !    subtle art and the code can easily fool you if you don't know
  !    what you're doing.  Make sure you understand *why* the code is 
  !    settling on a particular parameter value.  
  ! 4. Wavelength-dependent instrumental broadening is included but
  !    will not be accurate in the limit of modest-large redshift b/c
  !    this is implemented in the model restframe at code setup time

  !To Do: 
  !1. add SFH and MDF options

  !---------------------------------------------------------------!
  !---------------------------------------------------------------!

  USE alf_vars; USE alf_utils; USE mpi
  USE nr, ONLY : locate,powell
  USE ran_state, ONLY : ran_seed,ran_init

  IMPLICIT NONE

  !number of chain steps to print to file
  INTEGER, PARAMETER :: nmcmc=100
  !sampling of the walkers for printing
  INTEGER, PARAMETER :: nsample=1
  !length of chain burn-in
  INTEGER, PARAMETER :: nburn=1000  !1E4 seems to be good enough
  !number of walkers
  INTEGER, PARAMETER :: nwalkers=512 !1024

  !start w/ powell minimization?
  INTEGER, PARAMETER :: dopowell=1
  !Powell iteration tolerance
  REAL(DP), PARAMETER :: ftol=0.1
  INTEGER, PARAMETER :: test_time=0
  INTEGER  :: i,j,k,totacc=0,iter=30,npos
  REAL(DP) :: velz,msto,minchi2=huge_number,fret,wdth,bret=huge_number
  REAL(DP), DIMENSION(nl)   :: mspec=0.0,mspecmw=0.0,lam=0.0
  REAL(DP), DIMENSION(nfil) :: m2l=0.0,m2lmw=0.0
  REAL(DP), DIMENSION(npar) :: oposarr=0.,bposarr=0.0
  REAL(DP), DIMENSION(npar,nwalkers) :: mpiposarr=0.0
  REAL(DP), DIMENSION(3,npar+2*nfil) :: runtot=0.0
  REAL(DP), DIMENSION(npar,npar)     :: xi=0.0
  CHARACTER(10) :: time
  REAL(SP)      :: time2
  REAL(SP), DIMENSION(2) :: dumt
  CHARACTER(50) :: file='',tag=''
  TYPE(PARAMS)  :: opos,prlo,prhi,bpos
  !the next three definitions are for emcee
  REAL(DP), DIMENSION(npar,nwalkers) :: pos_emcee_in,pos_emcee_out
  REAL(DP), DIMENSION(nwalkers)      :: lp_emcee_in,lp_emcee_out,lp_mpi
  INTEGER,  DIMENSION(nwalkers)      :: accept_emcee

  !variables for MPI
  INTEGER :: ierr,taskid,ntasks,received_tag,status(MPI_STATUS_SIZE)
  INTEGER :: KILL=99,BEGIN=0
  LOGICAL :: wait=.TRUE.
  INTEGER, PARAMETER :: masterid=0
 
  !---------------------------------------------------------------!
  !---------------------------Setup-------------------------------!
  !---------------------------------------------------------------!

  !flag determining the level of complexity
  !0=full, 1=simple, 2=super-simple.  See sfvars for details
  fit_type   = 1
  !dont fit transmission function in cases where the input
  !spectrum has already been de-redshifted to ~0.0
  fit_trans  = 1
  !fit two-part power-law IMF if fit_oneimf=0
  fit_oneimf = 0

  !set low upper prior limits to kill these parameters
  !prhi%logm7g   = -5.0
  !prhi%loghot   = -5.0
  !prhi%logtrans = -5.0
  !prhi%logfy    = -5.0
  !prhi%logemline_h    = -5.0
  !prhi%logemline_oiii = -5.0
  !prhi%logemline_sii  = -5.0
  !prhi%logemline_nii  = -5.0
  !prhi%logemline_ni   = -5.0

  IF (fit_type.EQ.1.OR.fit_type.EQ.2) mwimf=1

  ! Initialize MPI, and get the total number of processes and
  ! your process number
  CALL MPI_INIT( ierr )
  CALL MPI_COMM_RANK( MPI_COMM_WORLD, taskid, ierr )
  CALL MPI_COMM_SIZE( MPI_COMM_WORLD, ntasks, ierr )

  !initialize the random number generator
  !set each task to sleep for a different length of time
  !so that each task has its own unique random number seed
  CALL SLEEP(taskid)
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

  IF (taskid.EQ.masterid) THEN
     !write some important variables to screen
     WRITE(*,*) 
     WRITE(*,'(" ************************************")') 
     WRITE(*,'("   dopowell  =",I2)') dopowell
     WRITE(*,'("   fit_type  =",I2)') fit_type
     WRITE(*,'("      mwimf  =",I2)') mwimf
     WRITE(*,'("  age-dep Rf =",I2)') use_age_dep_resp_fcns
     WRITE(*,'("  Nwalkers   = ",I5)') nwalkers
     WRITE(*,'("  Nburn      = ",I6)') nburn
     WRITE(*,'("  Nchain     = ",I5)') nmcmc
     WRITE(*,'("  filename   = ",A)') TRIM(file)//TRIM(tag)
     WRITE(*,'(" ************************************")') 
     CALL DATE_AND_TIME(TIME=time)
     CALL DTIME(dumt,time2)
     WRITE(*,*) 
     WRITE(*,*) 'Start Time '//time(1:2)//':'//time(3:4)
     
  ENDIF

  !read in the data and wavelength boundaries
  CALL READ_DATA(file)

  !read in the SSPs and bandpass filters
  CALL SETUP()
  lam = sspgrid%lam

  !we only compute things up to 500A beyond the input fit region
  nl_fit = MIN(MAX(locate(lam,l2(nlint)+500.0),1),nl)
  !nl_fit = MIN(MAX(locate(lam,15000.d0),1),nl)

  !define the log wavelength grid used in velbroad.f90
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
  CALL SET_PINIT_PRIORS(opos,prlo,prhi)
  !convert the structures into their equivalent arrays
  CALL STR2ARR(1,prlo,prloarr)   !str->arr
  CALL STR2ARR(1,prhi,prhiarr)   !str->arr

  ! The worker's only job is to calculate the value of a function
  ! after receiving a parameter vector.
  IF (taskid.NE.masterid) THEN
     
     ! Start event loop
     DO WHILE (wait)

        ! Look for data from the master. This call can accept up
        ! to ``nwalkers`` paramater positions, but it expects
        ! that the actual number of positions is smaller and is
        ! given by the MPI_TAG.  This call does not return until
        ! a set of parameter vectors is received
        CALL MPI_RECV(npos, 1, MPI_INTEGER, &
             masterid, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
        received_tag = status(MPI_TAG)
        IF ((received_tag.EQ.KILL).OR.(npos.EQ.0)) EXIT
        CALL MPI_RECV(mpiposarr(1,1), npos*npar, MPI_DOUBLE_PRECISION, &
             masterid, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
   
        IF (taskid.EQ.1.AND.test_time.EQ.1) THEN
           CALL DATE_AND_TIME(TIME=time)
           WRITE(*,*) '1 Time '//time(1:2)//':'//time(3:4)//':'&
                //time(5:9),npos,taskid
        ENDIF

        !Calculate the probability for these parameter positions
        DO k=1,npos
           lp_mpi(k) = -0.5*func(mpiposarr(:,k))
        ENDDO

         IF (taskid.EQ.1.AND.test_time.EQ.1) THEN
           CALL DATE_AND_TIME(TIME=time)
           WRITE(*,*) '2 Time '//time(1:2)//':'//time(3:4)//':'&
                //time(5:9),npos,taskid
        ENDIF
             
        !Send it back to the master
        CALL MPI_SEND(lp_mpi(1), npos, MPI_DOUBLE_PRECISION, &
             masterid, BEGIN, MPI_COMM_WORLD, ierr)

     ENDDO

  ENDIF
 
  !this is the master process
  IF (taskid.EQ.masterid) THEN
 
     WRITE(*,'("  Fitting ",I1," wavelength intervals")') nlint
     IF (l2(nlint).GT.lam(nl).OR.l1(1).LT.lam(1)) THEN
        WRITE(*,*) 'ERROR: wavelength boundaries exceed model wavelength grid'
        WRITE(*,'(4F8.1)') l2(nlint),lam(nl),l1(1),lam(1)
        STOP
     ENDIF

     !make an initial estimate of the redshift
     !we do this to help Powell minimization
     WRITE(*,*) ' Finding redshift...'
     IF (file(1:4).EQ.'cdfs') THEN
        velz = 0.0 
     ELSE 
        velz = getvelz()
     ENDIF
     opos%velz = velz
     WRITE(*,'("    best velocity: ",F7.1)') velz

     CALL STR2ARR(1,opos,oposarr)   !str->arr

     !initialize the random number generator
     !why is this being done here again?
     CALL INIT_RANDOM_SEED()

     IF (1.EQ.0) THEN
        opos%logemline_h    = -8.0
        opos%logemline_oiii = -8.0
        opos%logemline_nii  = -8.0
        opos%logemline_sii  = -8.0
        opos%logemline_ni   = -8.0
        opos%logage=1.13
        opos%logfy=-5.0
        opos%logm7g=-5.0
        opos%loghot=-5.0
        opos%imf1=1.3
        opos%imf2=2.3
        opos%zh=0.0
        opos%teff=0.0
        msto = MIN(MAX(10**(msto_fit0+msto_fit1*opos%logage),0.8),3.)  
        CALL GETMODEL(opos,mspecmw,mw=1)     !get spectrum for MW IMF
        CALL GETM2L(msto,lam,mspecmw,opos,m2lmw,mw=1) !compute M/L_MW
        !CALL GETM2L(msto,lam,10**sspgrid%logfkrpa(:,nage,nzmet-1),opos,m2lmw,mw=1)
        write(*,'(2F7.2)') m2lmw(1:2)
        CALL GETMODEL(opos,mspec)
        CALL GETM2L(msto,lam,mspec,opos,m2l) ! compute M/L
        write(*,'(2F7.2)') m2l(1:2)
        stop
     ENDIF


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
           CALL SET_PINIT_PRIORS(opos,prlo,prhi,velz=velz)
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
        
        WRITE(*,'("    best velocity: ",F7.1)') opos%velz
        WRITE(*,'("    best sigma:    ",F6.1)') opos%sigma
        WRITE(*,'("    best age:      ",F6.1)') 10**opos%logage
        WRITE(*,'("    best [Z/H]:    ",F6.1)') opos%zh
        
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
        CALL SET_PINIT_PRIORS(opos,prlo,prhi,velz=velz)
        CALL STR2ARR(1,opos,pos_emcee_in(:,j))

        IF (dopowell.EQ.1) THEN
           !use the best-fit position from Powell, with small
           !random offsets to set up all the walkers, but only
           !do this for the params actually fit in Powell!
           !the first two params are velz and sigma so give them
           !larger variation.
           DO i=1,npowell
              IF (i.LE.2) wdth = 10.0
              IF (i.GT.2) wdth = 0.1
              pos_emcee_in(i,j) = bposarr(i) + wdth*(2.*myran()-1.0)
              IF (pos_emcee_in(i,j).LE.prloarr(i)) &
                   pos_emcee_in(i,j)=prloarr(i)+wdth
              IF (pos_emcee_in(i,j).GE.prhiarr(i)) &
                   pos_emcee_in(i,j)=prhiarr(i)-wdth
           ENDDO
        ENDIF

        !Compute the initial log-probability for each walker
        lp_emcee_in(j) = -0.5*func(pos_emcee_in(:, j))
   
        IF (-2.*lp_emcee_in(j).GE.huge_number/2.) THEN
           WRITE(*,*) 'ALF ERROR: initial lnp out of bounds!'
           STOP
        ENDIF

     ENDDO

     !burn-in
     WRITE(*,*) '   burning in...'
     WRITE(*,'(A)',advance='no') '      Progress:'
     DO i=1,nburn
        CALL EMCEE_ADVANCE_MPI(npar,nwalkers,2.d0,pos_emcee_in,&
             lp_emcee_in,pos_emcee_out,lp_emcee_out,accept_emcee,ntasks-1)
        pos_emcee_in = pos_emcee_out
        lp_emcee_in  = lp_emcee_out
        IF (i.EQ.nburn/4.*1) THEN
           WRITE (*,'(A)',advance='no') ' ...25%'
           CALL FLUSH()
        ENDIF
        IF (i.EQ.nburn/4.*2) THEN
           WRITE (*,'(A)',advance='no') '...50%'
           CALL FLUSH()
        ENDIF
        IF (i.EQ.nburn/4.*3) THEN
           WRITE (*,'(A)',advance='no') '...75%'
           CALL FLUSH()
        ENDIF
     ENDDO
     WRITE (*,'(A)') '...100%'
     CALL FLUSH()
     
     !Run a production chain
     WRITE(*,*) '   production run...'
     DO i=1,nmcmc
        
        CALL EMCEE_ADVANCE_MPI(npar,nwalkers,2.d0,pos_emcee_in,&
             lp_emcee_in,pos_emcee_out,lp_emcee_out,accept_emcee,ntasks-1)
        pos_emcee_in = pos_emcee_out
        lp_emcee_in  = lp_emcee_out
        totacc = totacc + SUM(accept_emcee)
        
        DO j=1,nwalkers,nsample
           
           CALL STR2ARR(2,opos,pos_emcee_in(:,j)) !arr->str
           
           !kill the emission lines for computing M/L
           !since unconstrained lines can really mess up R,I bands
           opos%logemline_h    = -8.0
           opos%logemline_oiii = -8.0
           opos%logemline_nii  = -8.0
           opos%logemline_sii  = -8.0
           opos%logemline_ni   = -8.0

           !compute the main sequence turn-off mass
           !NB: Need to update this for other metallicities
           msto = MIN(MAX(10**(msto_fit0+msto_fit1*opos%logage),0.8),3.)           
           CALL GETMODEL(opos,mspecmw,mw=1)     !get spectrum for MW IMF
           CALL GETM2L(msto,lam,mspecmw,opos,m2lmw,mw=1) !compute M/L_MW
           
           IF (mwimf.EQ.0) THEN
              CALL GETMODEL(opos,mspec)
              CALL GETM2L(msto,lam,mspec,opos,m2l) ! compute M/L
           ELSE
              m2l = m2lmw
           ENDIF
           
           IF (fit_type.EQ.1) THEN
              !these parameters aren't actually being updated
              pos_emcee_in(nparsimp+1:,j) = 0.0 
           ELSE IF (fit_type.EQ.2) THEN
              !these parameters aren't actually being updated
              pos_emcee_in(npowell+1:,j) = 0.0 
           ENDIF
           
           !write the chain element to file
           WRITE(12,'(ES12.5,1x,F11.4,99(F9.4,1x))') &
                -2.0*lp_emcee_in(j),pos_emcee_in(:, j),m2l,m2lmw
           
           !keep the model with the lowest chi2
           IF (-2.0*lp_emcee_in(j).LT.minchi2) THEN
              bposarr = pos_emcee_in(:, j)
              minchi2 = -2.0*lp_emcee_in(j)
           ENDIF
           
           CALL UPDATE_RUNTOT(runtot,pos_emcee_in(:,j),m2l,m2lmw)
           
        ENDDO
        
     ENDDO
     
     !save the best position to the structure
     CALL STR2ARR(2,bpos,bposarr)
     bpos%chi2 = minchi2
     
     CLOSE(12)
     
     CALL DATE_AND_TIME(TIME=time)
     CALL DTIME(dumt,time2)
     WRITE(*,*) 'End Time   '//time(1:2)//':'//time(3:4)
     WRITE(*,'(" Elapsed Time: ",F6.2," hr")') time2/3600.
     WRITE(*,*) 
     WRITE(*,'("  Facc: ",F5.2)') REAL(totacc)/REAL(nmcmc*nwalkers)


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
     WRITE(14,'("#  Elapsed Time: ",F6.2," hr")') time2/3600.
     WRITE(14,'("#   dopowell  =",I2)') dopowell
     WRITE(14,'("#   fit_type  =",I2)') fit_type
     WRITE(14,'("#   fit_trans =",I2)') fit_trans
     WRITE(14,'("#   fit_poly  =",I2)') fit_poly
     WRITE(14,'("#      mwimf  =",I2)') mwimf
     WRITE(14,'("#  age-dep Rf =",I2)') use_age_dep_resp_fcns
     WRITE(14,'("#  Nwalkers   = ",I5)') nwalkers
     WRITE(14,'("#  Nburn      = ",I5)') nburn
     WRITE(14,'("#  Nchain     = ",I5)') nmcmc
     WRITE(14,'("#  Ncores     = ",I5)') ntasks
     WRITE(14,'("#  Facc: ",F5.2)') REAL(totacc)/REAL(nmcmc*nwalkers)
     WRITE(14,'(ES12.5,1x,F11.4,99(F9.4,1x))') bpos%chi2,runtot(2,:)/runtot(1,:)
     CLOSE(14)

     !write one sigma errors on parameters
     OPEN(15,FILE=TRIM(SPECFIT_HOME)//TRIM(OUTDIR)//&
          TRIM(file)//TRIM(tag)//'.errp',STATUS='REPLACE')
     WRITE(15,'(ES12.5,1x,F11.4,99(F9.4,1x))') 0.0, &
          SQRT( runtot(3,:)/runtot(1,:) - runtot(2,:)**2/runtot(1,:)**2 )
     CLOSE(15)

     WRITE(*,*)
     WRITE(*,'(" ************************************")') 

     !break the workers out of their event loops so they can close
     CALL FREE_WORKERS(ntasks-1)

  ENDIF

  CALL MPI_FINALIZE(ierr)
 

END PROGRAM ALF
