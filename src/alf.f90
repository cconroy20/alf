PROGRAM ALF

  !  Main program to fit the absorption line spectrum, or indices,
  !  of a quiescent (>1 Gyr) stellar population

  ! Some important points to keep in mind:
  ! 1. The prior bounds on the parameters are specified in set_pinit_priors. 
  !    Always make sure that the output parameters are not hitting a prior.
  ! 2. Make sure that the chain is converged in all relevant parameters
  !    by plotting the chain trace (parameter vs. chain step).
  ! 3. Do not use this code blindly.  Fitting spectra is a 
  !    subtle art and the code can easily fool you if you don't know
  !    what you're doing.  Make sure you understand *why* the code is 
  !    settling on a particular parameter value.  
  ! 4. Wavelength-dependent instrumental broadening can be included but
  !    will not be accurate in the limit of modest-large redshift b/c
  !    this is implemented in the model restframe at code setup time
  ! 5. The code can fit for the atmospheric transmission function but
  !    this will only work if the input data are in the original 
  !    observed frame; i.e., not de-redshifted.
  ! 6. I've found that Nwalkers=1024 and Nburn=~10,000 seems to
  !    generically yield well-converged solutions, but you should test
  !    this yourself by fitting mock data generated with write_a_model

  ! To Do:
  ! 1. Let the Fe-peak elements track Fe in simple mode
  ! 2. Force both young and old components to have the same abundance pattern
  
  !---------------------------------------------------------------------!
  !---------------------------------------------------------------------!

  USE alf_vars; USE alf_utils; USE mpi
  USE nr, ONLY : locate,powell,sort,gasdev
  USE ran_state, ONLY : ran_seed,ran_init

  IMPLICIT NONE

  !number of chain steps to print to file
  INTEGER, PARAMETER :: nmcmc=100
  !inverse sampling of the walkers for printing
  !NB: setting this to >1 currently results in errors in the *sum outputs
  INTEGER, PARAMETER :: nsample=1
  !length of chain burn-in
  INTEGER, PARAMETER :: nburn=1000
  !number of walkers
  INTEGER, PARAMETER :: nwalkers=256 !512 
  !save the chain outputs to file and the model spectra
  INTEGER, PARAMETER :: print_mcmc=1, print_mcmc_spec=0
  !option to re-initialize the parameters around the best-fit
  !solution and execute another round of burn-in
  INTEGER, PARAMETER :: moreburn=0
  
  !start w/ powell minimization?
  INTEGER, PARAMETER  :: dopowell=0
  !Powell iteration tolerance
  REAL(DP), PARAMETER :: ftol=0.1
  !if set, will print to screen timing of likelihood calls
  INTEGER, PARAMETER  :: test_time=0
  !number of Monte Carlo realizations of the noise for index errors
  INTEGER, PARAMETER  :: nmcindx=1000

  INTEGER  :: i,j,k,totacc=0,iter=30,npos,ml
  REAL(DP) :: velz,minchi2=huge_number,fret,wdth,bret=huge_number
  REAL(DP), DIMENSION(nl)   :: mspec=0.0,mspecmw=0.0,lam=0.0
  REAL(DP), DIMENSION(nfil) :: m2l=0.0,m2lmw=0.0
  REAL(DP), DIMENSION(npar) :: oposarr=0.,bposarr=0.0
  REAL(DP), DIMENSION(npar,nwalkers) :: mpiposarr=0.0
  REAL(DP), DIMENSION(3,npar+2*nfil) :: runtot=0.0
  REAL(DP), DIMENSION(npar+2*nfil)   :: cl2p5,cl16,cl50,cl84,cl97p5
  REAL(DP), DIMENSION(npar,npar)     :: xi=0.0
  REAL(DP), DIMENSION(npar+2*nfil,nwalkers*nmcmc/nsample) :: mcmcpar=0.0

  REAL(DP), DIMENSION(nwalkers*nmcmc/nsample) :: sortpos
  CHARACTER(10) :: time
  REAL(SP)      :: time2
  REAL(SP), DIMENSION(2) :: dumt
  CHARACTER(50) :: file='',tag=''
  TYPE(PARAMS)  :: opos,prlo,prhi,bpos,tpos
  REAL(DP)      :: sigma_indx,velz_indx
  REAL(DP), DIMENSION(ndat) :: gdev,tflx
  REAL(DP), DIMENSION(nmcindx,nindx) :: tmpindx=0.
  !REAL(SP), DIMENSION(nmcmc*nwalkers/nsample+1,nl) :: mspec_mcmc=0.0

  !variables for emcee
  REAL(DP), DIMENSION(npar,nwalkers) :: pos_emcee_in,pos_emcee_out
  REAL(DP), DIMENSION(nwalkers)      :: lp_emcee_in,lp_emcee_out,lp_mpi
  INTEGER,  DIMENSION(nwalkers)      :: accept_emcee

  !variables for MPI
  INTEGER :: ierr,taskid,ntasks,received_tag,status(MPI_STATUS_SIZE)
  INTEGER :: KILL=99,BEGIN=0
  LOGICAL :: wait=.TRUE.
  INTEGER, PARAMETER :: parentid=0
 
  !---------------------------------------------------------------!
  !---------------------------Setup-------------------------------!
  !---------------------------------------------------------------!

  !flag specifying if fitting indices or spectra
  fit_indices = 0

  !flag determining the level of complexity
  !0=full, 1=simple, 2=super-simple.  See sfvars for details
  fit_type = 0

  !fit h3 and h4 parameters
  fit_hermite = 0
  
  !type of IMF to fit
  !0=1PL, 1=2PL, 2=1PL+cutoff, 3=2PL+cutoff, 4=non-parametric IMF
  imf_type = 1

  !are the data in the original observed frame?
  observed_frame = 0

  !force a MW (Kroupa) IMF
  mwimf = 0

  !fit two-age SFH or not?  (only considered if fit_type=0)
  fit_two_ages = 0

  !IMF slope within the non-parametric IMF bins
  !0 = flat, 1 = Kroupa, 2 = Salpeter
  nonpimf_alpha = 2

  !turn on/off the use of an external tabulated M/L prior
  extmlpr = 0
  
  !change the prior limits to kill off these parameters
  prhi%logm7g = -5.0
  prhi%teff   =  2.0
  prlo%teff   = -2.0
  prhi%loghot = -4.0
  prlo%loghot = -6.0
  
  !mass of the young component should always be sub-dominant
  prhi%logfy = -0.5

  !---------------------------------------------------------------!
  !--------------Do not change things below this line-------------!
  !---------------unless you know what you are doing--------------!
  !---------------------------------------------------------------!

  !regularize non-parametric IMF (always do this)
  nonpimf_regularize = 1

  !dont fit transmission function in cases where the input
  !spectrum has already been de-redshifted to ~0.0
  IF (observed_frame.EQ.0.OR.fit_indices.EQ.1) THEN
     fit_trans     =  0
     prhi%logtrans = -5.0
     prhi%logsky   = -5.0 
  ELSE
     fit_trans = 1
     !extra smoothing to the transmission spectrum.
     !if the input data has been smoothed by a gaussian
     !in velocity space, set the parameter below to that extra smoothing
     smooth_trans = 0.0
  ENDIF
 
  IF (ssp_type.EQ.'cvd') THEN
     !always limit the [Z/H] range for CvD since
     !these models are actually only at Zsol
     prhi%zh =  0.01
     prlo%zh = -0.01
     IF (imf_type.GT.1) THEN
        WRITE(*,*) 'ALF ERROR, ssp_type=cvd but imf>1'
        STOP
     ENDIF
  ENDIF

  IF (fit_type.EQ.1.OR.fit_type.EQ.2) mwimf=1

  !---------------------------------------------------------------!

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

  IF (taskid.EQ.parentid) THEN
     !write some important variables to screen
     WRITE(*,*) 
     WRITE(*,'(" ************************************")') 
     IF (fit_indices.EQ.1) THEN
        WRITE(*,'(" ***********Index Fitter*************")')
     ELSE
        WRITE(*,'(" **********Spectral Fitter***********")')
     ENDIF
     WRITE(*,'(" ************************************")') 
     WRITE(*,'("   ssp_type  =",A4)') ssp_type
     WRITE(*,'("   fit_type  =",I2)') fit_type
     WRITE(*,'("   imf_type  =",I2)') imf_type
     WRITE(*,'(" fit_hermite =",I2)') fit_hermite
     WRITE(*,'("fit_two_ages =",I2)') fit_two_ages
     IF (imf_type.EQ.4) &
          WRITE(*,'("   nonpimf   =",I2)') nonpimf_alpha
     WRITE(*,'("  obs_frame  =",I2)') observed_frame 
     WRITE(*,'("      mwimf  =",I2)') mwimf
     WRITE(*,'("  age-dep Rf =",I2)') use_age_dep_resp_fcns
     WRITE(*,'("    Z-dep Rf =",I2)') use_z_dep_resp_fcns
     WRITE(*,'("  Nwalkers   = ",I6)') nwalkers
     WRITE(*,'("  Nburn      = ",I6)') nburn
     WRITE(*,'("  Nchain     = ",I6)') nmcmc
     WRITE(*,'("  Ncores     = ",I6)') ntasks
     WRITE(*,'("  filename   = ",A)') TRIM(file)//TRIM(tag)
     WRITE(*,'(" ************************************")') 
     CALL DATE_AND_TIME(TIME=time)
     CALL DTIME(dumt,time2)
     WRITE(*,*) 
     WRITE(*,*) 'Start Time '//time(1:2)//':'//time(3:4)
     
  ENDIF

  !read in the data and wavelength boundaries
  CALL READ_DATA(file,sigma_indx,velz_indx)

  
  IF (fit_indices.EQ.1) THEN

     !fold in the approx data sigma into the "instrumental"
     data%ires = SQRT(data%ires**2+sigma_indx**2)

     !read in the SSPs and bandpass filters
     CALL SETUP()
     lam = sspgrid%lam

     prhi%logemline_h    = -5.0
     prhi%logemline_oii  = -5.0
     prhi%logemline_oiii = -5.0
     prhi%logemline_nii  = -5.0
     prhi%logemline_sii  = -5.0
     prhi%logemline_ni   = -5.0
     prhi%loghot         = -5.0
     prhi%logm7g         = -5.0
     prhi%teff           =  2.0
     prlo%teff           = -2.0
     !we dont use velocities or dispersions here, so this 
     !should be unnecessary, but haven't tested turning them off yet.
     prlo%velz           = -10.
     prhi%velz           =  10.
     prlo%sigma          = sigma_indx-10.
     prhi%sigma          = sigma_indx+10.

     !de-redshift, monte carlo sample the noise, and compute indices
     !NB: need to mask bad pixels!
     DO j=1,nmcindx
        CALL GASDEV(gdev(1:datmax))
        tflx(1:datmax) = linterp(data(1:datmax)%lam/(1+velz_indx),&
             data(1:datmax)%flx+gdev(1:datmax)*data(1:datmax)%err,&
             data(1:datmax)%lam)
        CALL GETINDX(data(1:datmax)%lam,tflx(1:datmax),tmpindx(j,:))
     ENDDO

     !compute mean indices and errors
     DO j=1,nindx
        IF (indx2fit(j).EQ.1) THEN
           data_indx(j)%indx = SUM(tmpindx(:,j))/nmcindx
           data_indx(j)%err  = SQRT( SUM(tmpindx(:,j)**2)/nmcindx - &
                (SUM(tmpindx(:,j))/nmcindx)**2 )
           !write(*,'(I2,2F6.2)') j,data_indx(j)%indx,data_indx(j)%err
        ELSE
           data_indx(j)%indx = 0.0
           data_indx(j)%err  = 999.
        ENDIF
     ENDDO

     nl_fit = nl

  ENDIF

  IF (fit_indices.EQ.0) THEN

     !read in the SSPs and bandpass filters
     CALL SETUP()
     lam = sspgrid%lam

     !interpolate the sky emission model onto the observed wavelength grid
     IF (observed_frame.EQ.1) THEN
        data(1:datmax)%sky = MAX(linterp(lsky,fsky,data(1:datmax)%lam),0.0)
     ELSE
        data%sky = tiny_number
     ENDIF
     data%sky = tiny_number

     !we only compute things up to 500A beyond the input fit region
     nl_fit = MIN(MAX(locate(lam,l2(nlint)+500.0),1),nl)

     !define the log wavelength grid used in velbroad.f90
     dlstep = (LOG(sspgrid%lam(nl_fit))-LOG(sspgrid%lam(1)))/nl_fit
     DO i=1,nl_fit
        lnlam(i) = i*dlstep+LOG(sspgrid%lam(1))
     ENDDO

     !masked regions have wgt=0.0.  We'll use wgt as a pseudo-error
     !array in contnormspec, so turn these into large numbers
     data%wgt = MIN(1/(data%wgt+tiny_number),huge_number)
     !fold the masked regions into the errors
     data%err = MIN(data%err*data%wgt, huge_number)

  ENDIF

  !set initial params, step sizes, and prior ranges
  CALL SET_PINIT_PRIORS(opos,prlo,prhi)
  !convert the structures into their equivalent arrays
  CALL STR2ARR(1,prlo,prloarr)   !str->arr
  CALL STR2ARR(1,prhi,prhiarr)   !str->arr

  
  ! The worker's only job is to calculate the value of a function
  ! after receiving a parameter vector.
  IF (taskid.NE.parentid) THEN
     
     ! Start event loop
     DO WHILE (wait)

        ! Look for data from the parent. This call can accept up
        ! to ``nwalkers`` paramater positions, but it expects
        ! that the actual number of positions is smaller and is
        ! given by the MPI_TAG.  This call does not return until
        ! a set of parameter vectors is received
        CALL MPI_RECV(npos, 1, MPI_INTEGER, &
             parentid, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
        received_tag = status(MPI_TAG)
        IF ((received_tag.EQ.KILL).OR.(npos.EQ.0)) EXIT
        CALL MPI_RECV(mpiposarr(1,1), npos*npar, MPI_DOUBLE_PRECISION, &
             parentid, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
   
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
             
        !Send it back to the parent
        CALL MPI_SEND(lp_mpi(1), npos, MPI_DOUBLE_PRECISION, &
             parentid, BEGIN, MPI_COMM_WORLD, ierr)

     ENDDO

  ENDIF
 
  !this is the parent process
  IF (taskid.EQ.parentid) THEN
 
     !for testing
     IF (1.EQ.0) THEN
        tpos%logage = 1.143
        tpos%imf1   = 3.32
        tpos%imf2   = 2.76
        tpos%imf3   = 0.08
        
        CALL GETMODEL(tpos,mspecmw,mw=1)     !get spectrum for MW IMF
        CALL GETM2L(lam,mspecmw,tpos,m2lmw,mw=1) !compute M/L_MW
        write(*,'(A10,2F7.2)') 'M/L(MW)=', m2lmw(1:2)
        CALL GETMODEL(tpos,mspec)
        CALL GETM2L(lam,mspec,tpos,m2l)
        write(*,'(A10,2F7.2)') 'M/L=', m2l(1:2)
        CALL FREE_WORKERS(ntasks-1)
        CALL MPI_FINALIZE(ierr)
        STOP
     ENDIF

     IF (fit_indices.EQ.0) THEN

        WRITE(*,'("  Fitting ",I1," wavelength intervals")') nlint
        IF (l2(nlint).GT.lam(nl).OR.l1(1).LT.lam(1)) THEN
           WRITE(*,*) 'ERROR: wavelength boundaries exceed model wavelength grid'
           WRITE(*,'(4F8.1)') l2(nlint),lam(nl),l1(1),lam(1)
           STOP
        ENDIF

        !make an initial estimate of the redshift
        IF (file(1:4).EQ.'cdfs'.OR.file(1:5).EQ.'legac') THEN
           WRITE(*,*) 'Setting initial cz to 0.0'
           velz = 0.0 
        ELSE IF (file(1:4).EQ.'df44') THEN
           velz = 6280.00
        ELSE 
           WRITE(*,*) ' Fitting cz...'
           velz = getvelz()
           IF (velz.LT.prlo%velz.OR.velz.GT.prhi%velz) THEN
              WRITE(*,*) 'cz out of prior bounds, setting to 0.0'
              velz = 0.0
           ENDIF
        ENDIF
        opos%velz = velz
        WRITE(*,'("    cz= ",F7.1," (z=",F6.3,")")') &
             opos%velz, opos%velz/3E5
        
     ENDIF

     CALL STR2ARR(1,opos,oposarr)   !str->arr

     !initialize the random number generator
     !why is this being done here again?
     CALL INIT_RANDOM_SEED()

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
   
        !check for initialization errors
        IF (-2.*lp_emcee_in(j).GE.huge_number/2.) THEN
           WRITE(*,*) 'ALF ERROR: initial lnp out of bounds!', j
           DO i=1,npar
              IF (pos_emcee_in(i,j).GT.prhiarr(i).OR.&
                   pos_emcee_in(i,j).LT.prloarr(i)) THEN
                 WRITE(*,*) i, pos_emcee_in(i,j), prloarr(i), prhiarr(i)
              ENDIF
           ENDDO
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


     !burn-in V2
     IF (moreburn.EQ.1) THEN 

        ml    = MAXLOC(lp_emcee_in,1)    
        bposarr  = pos_emcee_in(:,ml) 

        DO j=1,nwalkers
           DO i=1,npar
              IF (i.LE.2) wdth = 10.0
              IF (i.GT.2) wdth = 0.05
              pos_emcee_in(i,j) = bposarr(i) + wdth*(2.*myran()-1.0)
              IF (pos_emcee_in(i,j).LE.prloarr(i)) &
                   pos_emcee_in(i,j)=prloarr(i)+wdth
              IF (pos_emcee_in(i,j).GE.prhiarr(i)) &
                   pos_emcee_in(i,j)=prhiarr(i)-wdth
           ENDDO
           
           !Compute the initial log-probability for each walker
           lp_emcee_in(j) = -0.5*func(pos_emcee_in(:, j))
        ENDDO
        

        WRITE(*,*) '   burning in (V2)...'
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
        
     END IF


     
     !Run a production chain
     WRITE(*,*) '   production run...'
     
     IF (print_mcmc.EQ.1) THEN
        !open output file
        OPEN(12,FILE=TRIM(ALF_HOME)//TRIM(OUTDIR)//&
             TRIM(file)//TRIM(tag)//'.mcmc',STATUS='REPLACE')
     ENDIF

     DO i=1,nmcmc
 
        CALL EMCEE_ADVANCE_MPI(npar,nwalkers,2.d0,pos_emcee_in,&
             lp_emcee_in,pos_emcee_out,lp_emcee_out,accept_emcee,ntasks-1)
        pos_emcee_in = pos_emcee_out
        lp_emcee_in  = lp_emcee_out
        totacc       = totacc + SUM(accept_emcee)
        
        DO j=1,nwalkers,nsample
           
           CALL STR2ARR(2,opos,pos_emcee_in(:,j)) !arr->str
           
           !turn off various parameters for computing M/L
           opos%logemline_h    = -8.0
           opos%logemline_oii  = -8.0
           opos%logemline_oiii = -8.0
           opos%logemline_nii  = -8.0
           opos%logemline_sii  = -8.0
           opos%logemline_ni   = -8.0
           opos%logtrans       = -8.0

           !compute the main sequence turn-off mass vs. t and Z
           CALL GETMODEL(opos,mspecmw,mw=1)     !get spectrum for MW IMF
           CALL GETM2L(lam,mspecmw,opos,m2lmw,mw=1) !compute M/L_MW

           IF (mwimf.EQ.0) THEN
              CALL GETMODEL(opos,mspec)
              CALL GETM2L(lam,mspec,opos,m2l) ! compute M/L
           ELSE
              m2l   = m2lmw
              mspec = mspecmw
           ENDIF
           
           !save each model spectrum
           !mspec_mcmc(1+j+(i-1)*nwalkers/nsample,:) = mspec

           !these parameters aren't actually being updated
           IF (fit_indices.EQ.1) THEN
              pos_emcee_in(1,j) = 0.0
              pos_emcee_in(2,j) = sigma_indx
           ENDIF
           IF (fit_type.EQ.1) THEN
              pos_emcee_in(nparsimp+1:,j) = 0.0
           ELSE IF (fit_type.EQ.2) THEN
              pos_emcee_in(npowell+1:,j) = 0.0 
           ENDIF
           
           IF (print_mcmc.EQ.1) THEN
              !write the chain element to file
              WRITE(12,'(ES12.5,1x,99(F11.4,1x))') &
                   -2.0*lp_emcee_in(j),pos_emcee_in(:,j),m2l,m2lmw
           ENDIF
           
           !keep the model with the lowest chi2
           IF (-2.0*lp_emcee_in(j).LT.minchi2) THEN
              bposarr = pos_emcee_in(:,j)
              minchi2 = -2.0*lp_emcee_in(j)
           ENDIF
           
           CALL UPDATE_RUNTOT(runtot,pos_emcee_in(:,j),m2l,m2lmw)
           
           !save each chain element
           mcmcpar(1:npar,j+(i-1)*nwalkers/nsample) = pos_emcee_in(:,j)
           mcmcpar(npar+1:npar+nfil,j+(i-1)*nwalkers/nsample)        = m2l
           mcmcpar(npar+nfil+1:npar+2*nfil,j+(i-1)*nwalkers/nsample) = m2lmw

        ENDDO
        
     ENDDO
     
     IF (print_mcmc.EQ.1) CLOSE(12)

     !save the best position to the structure
     CALL STR2ARR(2,bpos,bposarr)  !arr->str
     bpos%chi2 = minchi2

     !compute CLs
     DO i=1,npar+2*nfil
        sortpos = mcmcpar(i,:)
        CALL SORT(sortpos)
        cl2p5(i)  = sortpos(INT(0.025*nwalkers*nmcmc/nsample))
        cl16(i)   = sortpos(INT(0.160*nwalkers*nmcmc/nsample))
        cl50(i)   = sortpos(INT(0.500*nwalkers*nmcmc/nsample))
        cl84(i)   = sortpos(INT(0.840*nwalkers*nmcmc/nsample))
        cl97p5(i) = sortpos(INT(0.975*nwalkers*nmcmc/nsample))
     ENDDO
          
     CALL DATE_AND_TIME(TIME=time)
     CALL DTIME(dumt,time2)
     WRITE(*,*) 'End Time   '//time(1:2)//':'//time(3:4)
     WRITE(*,'(" Elapsed Time: ",F5.2," hr")') time2/3600.
     WRITE(*,*) 
     WRITE(*,'("  facc: ",F6.3)') REAL(totacc)/REAL(nmcmc*nwalkers)


     !---------------------------------------------------------------!
     !--------------------Write results to file----------------------!
     !---------------------------------------------------------------!
     
     !write a binary file of the production chain spectra
     IF (print_mcmc_spec.EQ.1) THEN
      !  mspec_mcmc(1,:) = lam
      !  OPEN(11,FILE=TRIM(ALF_HOME)//TRIM(OUTDIR)//&
      !       TRIM(file)//TRIM(tag)//'.spec',FORM='UNFORMATTED',&
      !       STATUS='REPLACE',access='DIRECT',&
      !       recl=(1+nmcmc*nwalkers/nsample)*nl*4)
      !  WRITE(11,rec=1) mspec_mcmc
      !  CLOSE(11)
     ENDIF

     OPEN(13,FILE=TRIM(ALF_HOME)//TRIM(OUTDIR)//&
          TRIM(file)//TRIM(tag)//'.bestspec',STATUS='REPLACE')
     CALL STR2ARR(1,bpos,bposarr)
     !NB: the model written to file has the lowest chi^2
     fret = func(bposarr,spec=mspec,funit=13)
     CLOSE(13)
 
     !write mean of the posterior distributions
     OPEN(14,FILE=TRIM(ALF_HOME)//TRIM(OUTDIR)//&
          TRIM(file)//TRIM(tag)//'.sum',STATUS='REPLACE')
     WRITE(14,'("#   Elapsed Time: ",F6.2," hr")') time2/3600.
     WRITE(14,'("#    ssp_type  =",A4)') ssp_type
     WRITE(14,'("#    fit_type  =",I2)') fit_type
     WRITE(14,'("#    imf_type  =",I2)') imf_type
     WRITE(14,'("#  fit_hermite =",I2)') fit_hermite
     WRITE(14,'("# fit_two_ages =",I2)') fit_two_ages
     WRITE(14,'("#     nonpimf  =",I2)') nonpimf_alpha
     WRITE(14,'("#   obs_frame  =",I2)') observed_frame 
     WRITE(14,'("#    fit_poly  =",I2)') fit_poly
     WRITE(14,'("#       mwimf  =",I2)') mwimf
     WRITE(14,'("#   age-dep Rf =",I2)') use_age_dep_resp_fcns
     WRITE(14,'("#     Z-dep Rf =",I2)') use_z_dep_resp_fcns
     WRITE(14,'("#   Nwalkers   = ",I6)') nwalkers
     WRITE(14,'("#   Nburn      = ",I6)') nburn
     WRITE(14,'("#   Nchain     = ",I6)') nmcmc
     WRITE(14,'("#   Nsample    = ",I6)') nsample
     WRITE(14,'("#   Nwave      = ",I6)') nl
     WRITE(14,'("#   Ncores     = ",I6)') ntasks
     WRITE(14,'("#   facc: ",F6.3)') REAL(totacc)/REAL(nmcmc*nwalkers)
     WRITE(14,'("#   rows: mean posterior, pos(chi^2_min), 1 sigma errors, '//&
          '2.5%, 16%, 50%, 84%, 97.5% CL, lower priors, upper priors ")') 

     !write mean of posteriors
     WRITE(14,'(ES12.5,1x,99(F11.4,1x))') bpos%chi2,runtot(2,:)/runtot(1,:)

     !write position where chi^2=min
     WRITE(14,'(ES12.5,1x,99(F11.4,1x))') bpos%chi2,bposarr,m2l*0.0,m2lmw*0.0

     !write 1 sigma errors
     WRITE(14,'(ES12.5,1x,99(F11.4,1x))') 0.0, &
          SQRT( runtot(3,:)/runtot(1,:) - runtot(2,:)**2/runtot(1,:)**2 )

     !write 2.5%, 16%, 50%, 84%, 97.5% CL
     WRITE(14,'(ES12.5,1x,99(F11.4,1x))') 0.0, cl2p5
     WRITE(14,'(ES12.5,1x,99(F11.4,1x))') 0.0, cl16
     WRITE(14,'(ES12.5,1x,99(F11.4,1x))') 0.0, cl50
     WRITE(14,'(ES12.5,1x,99(F11.4,1x))') 0.0, cl84
     WRITE(14,'(ES12.5,1x,99(F11.4,1x))') 0.0, cl97p5

     !write lower/upper priors
     WRITE(14,'(ES12.5,1x,99(F11.4,1x))') 0.0,prloarr,m2l*0.0,m2lmw*0.0
     WRITE(14,'(ES12.5,1x,99(F11.4,1x))') 0.0,prhiarr,m2l*0.0,m2lmw*0.0

     CLOSE(14)

     WRITE(*,*)
     WRITE(*,'(" ************************************")') 

     !break the workers out of their event loops so they can close
     CALL FREE_WORKERS(ntasks-1)

  ENDIF

  CALL MPI_FINALIZE(ierr)
 

END PROGRAM ALF
