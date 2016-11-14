PROGRAM ALF

  !  Master program to fit the absorption line spectrum
  !  of a quiescent (>1 Gyr) stellar population

  ! Some important points to keep in mind:
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
  ! 5. The code can fit for the atmospheric transmission function but
  !    this will only work if the input data are in the original 
  !    observed frame; i.e., not de-redshifted.
  ! 6. I've found that Nwalkers=1024 and Nburn=~10,000 seems to
  !    generically yield well-converged solutions, but you should test
  !    this yourself by fitting mock data generated with write_a_model

  !---------------------------------------------------------------------!
  !---------------------------------------------------------------------!

  USE alf_vars; USE alf_utils; USE mpi
  USE nr, ONLY : locate,powell,sort
  USE ran_state, ONLY : ran_seed,ran_init

  IMPLICIT NONE

  !number of chain steps to print to file
  INTEGER, PARAMETER :: nmcmc=100
  !inverse sampling of the walkers for printing
  INTEGER, PARAMETER :: nsample=1
  !length of chain burn-in
  INTEGER, PARAMETER :: nburn=20000
  !number of walkers
  INTEGER, PARAMETER :: nwalkers=1020 !1020 (15), 1085 (31)
  !save the chain outputs to file
  INTEGER, PARAMETER :: print_mcmc=1

  !start w/ powell minimization?
  INTEGER, PARAMETER  :: dopowell=0
  !Powell iteration tolerance
  REAL(DP), PARAMETER :: ftol=0.1
  !if set, will print to screen timing of likelihood calls
  INTEGER, PARAMETER  :: test_time=0

  INTEGER  :: i,j,k,totacc=0,iter=30,npos
  REAL(DP) :: velz,msto,minchi2=huge_number,fret,wdth,bret=huge_number
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

  !variables for emcee
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
  fit_type = 0
  !type of IMF to fit
  !0=1PL, 1=2PL, 2=1PL+cutoff, 3=2PL+cutoff, 4=non-parametric IMF
  imf_type = 1
  !are the data in the original observed frame?
  observed_frame = 1
  !IMF slope within the non-parametric IMF bins
  nonpimf_alpha = 2.3
  !force MW IMF
  mwimf    = 0
  !fit two-age SFH or not?
  fit_two_ages = 1
  !regularize non-parametric IMF
  nonpimf_regularize = 1

  !set low upper prior limits to kill off these parameters
  prhi%logm7g = -5.0
  prhi%teff   =  2.0
  prlo%teff   = -2.0

  !correction factor between Salpeter and flat intrabin weights
  !for non-parametric IMF
  IF (nonpimf_alpha.EQ.0.0) THEN
     corr_salp_flat = 0.0
  ENDIF

  !dont fit transmission function in cases where the input
  !spectrum has already been de-redshifted to ~0.0
  IF (observed_frame.EQ.0) THEN
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

  IF (taskid.EQ.masterid) THEN
     !write some important variables to screen
     WRITE(*,*) 
     WRITE(*,'(" ************************************")') 
     WRITE(*,'("   ssp_type  =",A4)') ssp_type
     WRITE(*,'("   fit_type  =",I2)') fit_type
     WRITE(*,'("   imf_type  =",I2)') imf_type
     IF (imf_type.EQ.4) &
          WRITE(*,'("   nonpimf   =",F4.1)') nonpimf_alpha
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
  CALL READ_DATA(file)

  !read in the SSPs and bandpass filters
  CALL SETUP()
  lam = sspgrid%lam

  !interpolate the sky emission model onto the observed wavelength grid
  IF (observed_frame.EQ.1) THEN
     data(1:datmax)%sky = MAX(linterp(lsky,fsky,data(1:datmax)%lam),0.0)
  ELSE
     data%sky = tiny_number
  ENDIF

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

     !for testing
     IF (1.EQ.0) THEN
        tpos%logage = 1.0
        tpos%logfy  = -5.0
        tpos%logm7g = -5.0
        tpos%loghot = -5.0
        tpos%imf1 = 2.3
        tpos%imf2 = 2.3
        tpos%imf3 = 0.08
        tpos%imf4 = 0.0
        tpos%zh   = 0.0
        tpos%teff = 0.0
        msto = 10**(msto_t0+msto_t1*tpos%logage) * &
             ( msto_z0 + msto_z1*tpos%zh + msto_z2*tpos%zh**2 )
        CALL GETMODEL(tpos,mspecmw,mw=1)     !get spectrum for MW IMF
        CALL GETM2L(msto,lam,mspecmw,tpos,m2lmw,mw=1) !compute M/L_MW
        write(*,'("M/L=",2F7.2)') m2lmw(1:2)
        CALL GETMODEL(tpos,mspec)
        CALL GETM2L(msto,lam,mspec,tpos,m2l)
        write(*,'("M/L=",2F7.2)') m2l(1:2)
        STOP
     ENDIF


     !make an initial estimate of the redshift
     IF (file(1:4).EQ.'cdfs'.OR.file(1:5).EQ.'legac') THEN
        velz = 0.0 
     ELSE 
        WRITE(*,*) ' Finding redshift...'
        velz = getvelz()
     ENDIF
     opos%velz = velz
     WRITE(*,'("    cz= ",F7.1," (z=",F6.3,")")') &
          velz, velz/3E5

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
     
     !CALL STR2ARR(1,prlo,prloarr)   !str->arr
     !CALL STR2ARR(1,prhi,prhiarr)   !str->arr

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
           !now increase the prior range on the IMF
         !  prhi%imf1 =  4.0
         !  prlo%imf1 = -4.0
         !  prhi%imf2 =  4.0
         !  prlo%imf2 = -4.0
         !  prhi%imf3 =  4.0
         !  prlo%imf3 = -4.0
         !  prhi%imf4 =  4.0
         !  prlo%imf4 = -4.0
         !  !convert the structures into their equivalent arrays
         !  CALL STR2ARR(1,prlo,prloarr)   !str->arr
         !  CALL STR2ARR(1,prhi,prhiarr)   !str->arr
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
           opos%logemline_oiii = -8.0
           opos%logemline_nii  = -8.0
           opos%logemline_sii  = -8.0
           opos%logemline_ni   = -8.0
           opos%logtrans       = -8.0

           !compute the main sequence turn-off mass vs. t and Z
           msto = MAX(MIN(10**(msto_t0+msto_t1*opos%logage) * &
                (msto_z0+msto_z1*opos%zh+msto_z2*opos%zh**2),3.0),0.75)
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
     WRITE(*,'("  Facc: ",F5.2)') REAL(totacc)/REAL(nmcmc*nwalkers)


     !---------------------------------------------------------------!
     !--------------------Write results to file----------------------!
     !---------------------------------------------------------------!
     
     OPEN(13,FILE=TRIM(ALF_HOME)//TRIM(OUTDIR)//&
          TRIM(file)//TRIM(tag)//'.bestspec',STATUS='REPLACE')
     CALL STR2ARR(1,bpos,bposarr)
     !NB: the model written to file has the lowest chi^2
     fret = func(bposarr,spec=mspec,funit=13)
     CLOSE(13)
 
     !write mean of the posterior distributions
     OPEN(14,FILE=TRIM(ALF_HOME)//TRIM(OUTDIR)//&
          TRIM(file)//TRIM(tag)//'.sum',STATUS='REPLACE')
     WRITE(14,'("#  Elapsed Time: ",F6.2," hr")') time2/3600.
     WRITE(14,'("#   ssp_type  =",A4)') ssp_type
     WRITE(14,'("#   fit_type  =",I2)') fit_type
     WRITE(14,'("#   imf_type  =",I2)') imf_type
     WRITE(14,'("#    nonpimf  =",F4.1)') nonpimf_alpha
     WRITE(14,'("#  obs_frame  =",I2)') observed_frame 
     WRITE(14,'("#   fit_poly  =",I2)') fit_poly
     WRITE(14,'("#      mwimf  =",I2)') mwimf
     WRITE(14,'("#  age-dep Rf =",I2)') use_age_dep_resp_fcns
     WRITE(14,'("#    Z-dep Rf =",I2)') use_z_dep_resp_fcns
     WRITE(14,'("#  Nwalkers   = ",I6)') nwalkers
     WRITE(14,'("#  Nburn      = ",I6)') nburn
     WRITE(14,'("#  Nchain     = ",I6)') nmcmc
     WRITE(14,'("#  Ncores     = ",I6)') ntasks
     WRITE(14,'("#  Facc: ",F5.2)') REAL(totacc)/REAL(nmcmc*nwalkers)
     WRITE(14,'("#  rows: mean posterior, pos(chi^2_min), 1 sigma errors, '//&
          '2.5%, 16%, 50%, 84%, 97.5% CL, lower priors, upper priors ")') 
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
