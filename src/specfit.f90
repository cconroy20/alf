PROGRAM SPECFIT

  ! Master program to fit the continuum-normalized spectrum
  ! of an old (>1 Gyr), metal-rich ([Fe/H]>-0.3) stellar
  ! population to the CvD model.

  USE sfvars; USE nr, ONLY : gasdev,ran,locate,powell,ran1
  USE ran_state, ONLY : ran_seed,ran_init; USE sfutils
  IMPLICIT NONE

  !number of chain steps to run
  INTEGER, PARAMETER :: nmcmc=1E5
  !estimated burn-in length
  INTEGER, PARAMETER :: nburn=1E5
  !start w/ powell minimization?
  INTEGER, PARAMETER :: dopowell=1
  !total length of output mcmc file
  INTEGER, PARAMETER :: nmax=1E4

  !down-sample the output chains by this factor
  INTEGER, PARAMETER :: sample=nmcmc/nmax
  !Powell iteration tolerance
  REAL(DP), PARAMETER :: ftol=0.1
  INTEGER  :: i,j,k,totacc=0,stat,i1,i2,iter=30,vv
  REAL(DP) :: mass,mwmass,fret,bret=huge_number,deltachi2,velz,s2n,msto
  REAL(DP), DIMENSION(nl)   :: mspec,aspec,mspecmw,aspecmw
  REAL(DP), DIMENSION(nl)   :: dflx,mflx,lam, err1=1.0
  REAL(DP), DIMENSION(nfil) :: mag,m2l,m2lmw
  REAL(DP), DIMENSION(ndat) :: tlam
  REAL(DP), DIMENSION(npar) :: oposarr=0.,nposarr=0.,bposarr=0.0
  REAL(DP), DIMENSION(npar) :: granarr=0.,dsteparr=0.
  REAL(DP), DIMENSION(3,npar+2*nfil) :: runtot=0.0
  REAL(DP), DIMENSION(npar,npar) :: xi=0.0
  CHARACTER(10) :: time=''
  CHARACTER(50) :: file='',tag=''
  TYPE(PARAMS)  :: npos,opos,prlo,prhi,dstep,bpos
  TYPE(TDATA), DIMENSION(nl) :: idata
  
  !---------------------------------------------------------------!
  !---------------------------Setup-------------------------------!
  !---------------------------------------------------------------!


  !initialize the random number generator
  CALL INIT_RANDOM_SEED()

  IF (IARGC().LT.1) THEN
     file(1:5) = 'sdss6'
  ELSE
     CALL GETARG(1,file)
  ENDIF
  IF (IARGC().GT.1) THEN
     tag(1:1)='_'
     CALL GETARG(2,tag(2:))
  ENDIF

  IF (nmax.GT.nmcmc) THEN
     WRITE(*,*) 'ERROR: nmax > nmcmc'
     STOP
  ENDIF

  CALL date_and_time(TIME=time)
  WRITE(*,*) 
  WRITE(*,*) 'Start Time '//time(1:2)//':'//time(3:4)

  !read in the SSPs and bandpass filters
  CALL SFSETUP()
  lam = sspgrid%lam

  !read in the data to be fit
  CALL READ_DATA(file)
  !masked regions have wgt=0.0.  We'll use wgt as an error
  !array in contnormspec, so turn those into large numbers
  data%wgt = MIN(1/data%wgt,huge_number)
  !fold the masked regions into the errors
  data%err = data%err * data%wgt

  !determine step size based on S/N
  vv  = locate(data(1:datmax)%lam,5000.d0)
  !estimate S/N at 5000A based on a 50A window
  s2n = SUM(data(vv-20:vv+20)%flx/data(vv-20:vv+20)%err) / 51.d0
  !this is just an empirical relation between mean S/N and step size
  !that gives us the desired facc~25%
  mcstep = MIN(3E-5*(1.5E3/s2n)**(3./4),1E-3)

  !set initial params, step sizes, and prior ranges
  CALL SETUP_PARAMS(opos,dstep,prlo,prhi)

  !set wavelength limits
  CALL SET_WAVE_LIMITS(file)

  !convert the structures into their equivalent arrays
  CALL STR2ARR(1,dstep,dsteparr) !str->arr
  CALL STR2ARR(1,prlo,prloarr)   !str->arr
  CALL STR2ARR(1,prhi,prhiarr)   !str->arr
  CALL STR2ARR(1,opos,oposarr)   !str->arr

  !---------------------------------------------------------------!
  !---------------------Powell minimization-----------------------!
  !---------------------------------------------------------------!

  !make an initial estimate of the redshift
  !we have to do this otherwise Powell minimization sometimes
  !chooses incorrect solutions
  velz = getvelz()
  IF (file(1:5).EQ.'usher')   velz = 0.0
  IF (file(1:7).EQ.'mosfire') velz = 0.0
  
  IF (dopowell.EQ.1) THEN 

     DO j=1,2
        xi=0.0
        DO i=1,npar
           xi(i,i) = 1E-2
        ENDDO
        fret = huge_number
        CALL SETUP_PARAMS(opos,dstep,prlo,prhi,velz=velz)
        CALL STR2ARR(1,opos,oposarr) !str->arr
        CALL POWELL(oposarr,xi,ftol,iter,fret)
        CALL STR2ARR(2,opos,oposarr) !arr->str
        IF (fret.LT.bret) THEN
           bposarr = oposarr
           bpos    = opos
           bret    = fret
        ENDIF
     ENDDO

     !use the best-fit Powell position for the first MCMC position
     CALL STR2ARR(2,opos,bposarr) !arr->str
     !Powell doesn't seem to do a good first-guess for these
     !params, so set them to defaults
     opos%imf1 = 1.3
     opos%imf2 = 2.3
     opos%logage = LOG10(12.0)
     opos%cofe = 0.0
     opos%crfe = 0.0
     opos%mnfe = 0.0
     opos%nife = 0.0
     CALL STR2ARR(1,opos,oposarr) !str->arr

  ENDIF
 
  IF (maskem.EQ.1) THEN
     !now that we have a good guess of the redshift and velocity dispersion, 
     !mask out regions where emission line contamination may be a problem
     CALL MASKEMLINES(opos%velz,opos%sigma)
  ENDIF

  !---------------------------------------------------------------!
  !-----------------------------MCMC------------------------------!
  !---------------------------------------------------------------!

  !open output file
  OPEN(12,FILE=TRIM(SPECFIT_HOME)//TRIM(OUTDIR)//&
       TRIM(file)//TRIM(tag)//'.mcmc',STATUS='REPLACE')
   
  !run the chain
  DO j=1,nmcmc+nburn

     !take a step
     CALL GASDEV(granarr)
     nposarr = oposarr + granarr * dsteparr
     CALL STR2ARR(2,npos,nposarr) !arr->str

     IF (mwimf.EQ.1) THEN
        npos%imf1 = 1.3
        npos%imf2 = 2.3
        CALL STR2ARR(1,npos,nposarr) !str->arr
     ENDIF

     !get a new model and compute chi^2
     npos%chi2 = func(nposarr,spec=mspec)

     !keep the model with the lowest chi2
     IF (npos%chi2.LT.opos%chi2) bpos = npos
     deltachi2 = EXP(-(npos%chi2-opos%chi2)/2.)
     !accept the step?
     IF (myran().LT.deltachi2) THEN
        opos    = npos
        oposarr = nposarr
        IF (j.GT.nburn) totacc  = totacc+1
     ENDIF

     !write chain element to file, subsampling the main chain
     IF (MOD(j,sample).EQ.0.AND.j.GT.nburn) THEN
        !compute the main sequence turn-off mass
        msto=10**( msto_fit0 + msto_fit1*opos%logage )
        CALL GETM2L(msto,lam,mspec,opos,m2l) ! compute M/L
        CALL GETMODEL(opos,mspecmw,mw=1)     ! get spectra for MW IMF
        CALL GETM2L(msto,lam,mspecmw,opos,m2lmw,mw=1) !compute M/L_MW
        WRITE(12,'(ES11.5,99(F9.4,1x))') opos%chi2,oposarr,m2l,m2lmw
        CALL FLUSH(12)
    ENDIF

    !here we want to use all of the chain elements after Nburn
    !for a more reliable estimate of the variance, etc.
    IF (j.GT.nburn) THEN
        runtot(1,:)      = runtot(1,:)+1.
        runtot(2,1:npar) = runtot(2,1:npar) + oposarr
        runtot(3,1:npar) = runtot(3,1:npar) + oposarr**2
        runtot(2,npar+1:npar+nfil) = runtot(2,npar+1:npar+nfil)+m2l
        runtot(3,npar+1:npar+nfil) = runtot(3,npar+1:npar+nfil)+m2l**2
        runtot(2,npar+nfil+1:npar+2*nfil) = &
             runtot(2,npar+nfil+1:npar+2*nfil)+m2lmw
        runtot(3,npar+nfil+1:npar+2*nfil) = &
             runtot(3,npar+nfil+1:npar+2*nfil)+m2lmw**2
     ENDIF

  ENDDO

  CLOSE(12)

  CALL date_and_time(TIME=time)
  WRITE(*,*) '  End Time '//time(1:2)//':'//time(3:4)
  WRITE(*,*) 
  WRITE(*,'("  facc: ",F4.2)') REAL(totacc)/REAL(nmcmc)
  WRITE(*,*) 

  !---------------------------------------------------------------!
  !--------------------Write results to file----------------------!
  !---------------------------------------------------------------!

  !NB: the model written to file has the lowest chi^2
  CALL GETMODEL(bpos,mspec)  !get the model
  tlam      = data%lam / (1+bpos%velz/clight*1E5) !de-redshift the data
  idata%flx = linterp(tlam,data%flx,lam)
  idata%err = linterp(tlam,data%err,lam)
  idata%wgt = linterp(tlam,data%wgt,lam)

  !write best-fit model and data (cont-norm) to file
  OPEN(13,FILE=TRIM(SPECFIT_HOME)//TRIM(OUTDIR)//&
       TRIM(file)//TRIM(tag)//'.bestspec',STATUS='REPLACE')
  DO i=1,nlint
     IF (MAXVAL(tlam(1:datmax)).LT.l2(i)) CYCLE
     IF (MINVAL(tlam(1:datmax)).GT.l1(i)) CYCLE
     CALL CONTNORMSPEC(lam,idata%flx,idata%err,l1(i),l2(i),dflx)
     CALL CONTNORMSPEC(lam,mspec,idata%wgt,l1(i),l2(i),mflx)
     i1 = MIN(MAX(locate(lam,l1(i)),1),nl-1)
     i2 = MIN(MAX(locate(lam,l2(i)),2),nl)
     WRITE(*,'("  rms:",F5.2,"%")') &
          SQRT( SUM( (dflx(i1:i2)/mflx(i1:i2)-1)**2 ) / (i2-i1+1) )*100
     DO j=i1,i2
        WRITE(13,'(F9.2,3ES12.4)') lam(j),mflx(j),dflx(j),&
             idata(j)%flx/idata(j)%err
     ENDDO
  ENDDO
  CLOSE(13)
  
  !write best-fit parameters
  !here, "best-fit" means the mean of the posterior distributions
  OPEN(14,FILE=TRIM(SPECFIT_HOME)//TRIM(OUTDIR)//&
       TRIM(file)//TRIM(tag)//'.bestp',STATUS='REPLACE')
  WRITE(14,'(ES11.5,99(F9.3,1x))') bpos%chi2, runtot(2,:)/runtot(1,:)
  CLOSE(14)

  !write 1 sigma errors on parameters
  OPEN(14,FILE=TRIM(SPECFIT_HOME)//TRIM(OUTDIR)//&
       TRIM(file)//TRIM(tag)//'.errp',STATUS='REPLACE')
  WRITE(14,'(ES11.5,99(F9.3,1x))') 0.0, &
       SQRT( runtot(3,:)/runtot(1,:) - runtot(2,:)**2/runtot(1,:)**2 )
  CLOSE(14)

  WRITE(*,*)

END PROGRAM SPECFIT
