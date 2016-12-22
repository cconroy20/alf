MODULE ALF_VARS

  ! module to set up most arrays and variables 

  IMPLICIT NONE
  SAVE

!use the VCJ SSPs, else use CvD (only for imf1 option)
#define VCJ 1

  !directory for results (mcmc, bestspec, etc.)
  CHARACTER(100) :: OUTDIR='results/'

  INTEGER, PARAMETER :: DP = KIND(1.0D0)
  INTEGER, PARAMETER :: SP = KIND(1.0)

  !-------------------set various parameters---------------------!

  !flags for the user to choose.  These can also be set at the 
  !very beginning of alf.f90

  !0: fit the full model (IMF, all abundances, nuisance params, etc)
  !1: only fit velz, sigma, SSP age, Z, Fe,C,N,O,Mg,Si,Ca,Ti,Na
  !2: only fit velz, sigma, SSP age, Z
  INTEGER :: fit_type=0


  !turn on the use of age-dependent response functions
  INTEGER :: use_age_dep_resp_fcns=1
  !if above is set to 0, fix the response functions to this age (Gyr)
  REAL(DP) :: fix_age_dep_resp_fcns=10.0

  !turn on the use of Z-dependent response functions
  INTEGER :: use_z_dep_resp_fcns=1
  !if above is set to 0, fix the response functions to this [Z/H]
  REAL(DP) :: fix_z_dep_resp_fcns=0.0

  !flag to include transmission spectrum in fitting
  !even if flag is set, only included in full model
  INTEGER :: fit_trans=1

  !force the IMF to be a MW IMF if set
  !this is automatically assumed if fit_type=1,2
  INTEGER :: mwimf=0

  !flag to fit either a double-power law IMF or power-law + cutoff
  !0 = single power-law
  !1 = double power-law
  !2 = power-law + cutoff
  !3 = double power-law + cutoff
  INTEGER :: imf_type=1

  !are the data in the original observed frame?
  INTEGER :: observed_frame=1

  !IMF used to compute the element response functions
  CHARACTER(4), PARAMETER :: atlas_imf='krpa'  !'salp'

  !extra smoothing (km/s) of the transmission spectrum
  !if the input spectrum has been smoothed by an amount more than
  !the instrumental resolution, set the parameter below to that value
  REAL(DP) :: smooth_trans=0.0

  !flag used to tell the code if we are fitting in powell mode or not
  !this is set internally in the code 
  INTEGER :: powell_fitting = 0

  !IMF power-law slope within each bin for non-paramtric IMF
  REAL(DP) :: nonpimf_alpha = 2.3

  !fit two-component SFH. Also requires fit_type=0
  INTEGER :: fit_two_ages=1

  !regularize the non-parametric IMF
  INTEGER :: nonpimf_regularize=1

  !fit indices (only when using if.exe)
  INTEGER :: fit_indices=1

  !--------------------------------------------------------------!
  !  the options below have not been tested/used in a long time  !
  !  and so are effectively deprecated                           !
  !--------------------------------------------------------------!

  !fit a polynomial to the ratio of model and data
  !if zero, then both data and model are continuum divided
  INTEGER :: fit_poly=1
  !mask emission lines? (if 0, then the em lines are incl in the fit)
  INTEGER :: maskem=0
  !apply template error function? (only works for SDSS stacks)
  INTEGER :: apply_temperrfcn=0
  !flag to implement fake element response functions
  INTEGER :: fake_response=0
  !Turn off the IMF sensitivity at <7000A if this parameter is =1
  INTEGER :: blueimf_off=0
  !if set, compute velocity broadening via a simple method
  !rather than the proper convolution in log_lambda space
  !don't turn this on - the "correct" version is just as fast
  INTEGER :: velbroad_simple=0

  !--------------------------------------------------------------!
  !    the parameters below should not be modified unless you    !
  !    really know what you are doing!                           !
  !--------------------------------------------------------------!
  
#if (VCJ)
  !VCJ models
  CHARACTER(3), PARAMETER :: ssp_type='vcj'
  INTEGER, PARAMETER :: nzmet3 = 3
#else
  !CvD models
  CHARACTER(3), PARAMETER :: ssp_type='cvd'
  INTEGER, PARAMETER :: nzmet3 = 1
#endif

  !nstart and nend allow us to use only a subset of 
  !the full wavelength array
  INTEGER, PARAMETER :: nstart = 100 !100   ! 0.36 um
  INTEGER, PARAMETER :: nend   = 5830 !10566 ! 5830  ! 1.10 um
  !number of spectral elements in SSPs
  INTEGER, PARAMETER :: nl = nend-nstart+1
  !number actually used over the range to be fit
  INTEGER :: nl_fit=nl
  !(max) number of wavelength intervals
  INTEGER, PARAMETER :: nlint_max = 10
  !actual number of wavelength intervals, determined at run time
  INTEGER :: nlint = 0
  !total number of emission lines
  INTEGER, PARAMETER :: neml = 11
  !number of parameters
  INTEGER, PARAMETER :: npar = 43
  !number of ages in the empirical SSP grid
  INTEGER, PARAMETER :: nage = 7
  !number of metallicities in the empirical SSP grid
  INTEGER, PARAMETER :: nzmet = 5
  !number of parameters used when fitting in Powell model
  !or in the super-simple mode (fit_type=2)
  INTEGER, PARAMETER :: npowell = 4
  !number of ages in the response functions
  INTEGER, PARAMETER :: nage_rfcn = 5
  !number of IMF values in the SSP grid
  INTEGER, PARAMETER :: nimf_full=16, nmcut=8, nimfoff=2, nimfnp=9
  INTEGER, PARAMETER :: nimf=nimf_full-nimfoff
  !max degree of polynomial used for continuum fitting
  INTEGER, PARAMETER :: npolymax = 20
  !wavelength interval used to determine polynomial degree
  REAL, PARAMETER :: poly_dlam = 100.0
  !max number of data wavelength points
  INTEGER, PARAMETER :: ndat = 30000
  !total number of parameters in the simple model
  INTEGER, PARAMETER :: nparsimp = 14
  !number of indices defined in allindices.dat
  INTEGER, PARAMETER :: nindx=20
  !number of filters
  INTEGER, PARAMETER :: nfil=3
  !number of hot stars
  INTEGER, PARAMETER :: nhot=4
  !mag of sun in r,I,K filters (AB mag)
  REAL(DP), PARAMETER, DIMENSION(3) :: magsun = (/4.64,4.52,5.14/)
  !mag of sun in r,I,J filters (AB mag)
  !REAL(DP), PARAMETER, DIMENSION(3) :: magsun = (/4.64,4.52,4.56/)
  !lower and upper limits for the IMF,
  !except when the IMF cutoff parameterization is used
  REAL(DP), PARAMETER :: imflo=0.08,imfhi=100.0
  !power-law slopes for a Kroupa IMF
  REAL(DP), PARAMETER :: krpa_imf1=1.3,krpa_imf2=2.3,krpa_imf3=2.3
  !log(imf5) for non-parametric IMF
  REAL(DP), PARAMETER :: imf5 = 0.0
  !linear fit to log(age) vs. log(MS TO mass)
  REAL(DP), PARAMETER :: msto_t0=0.33250847,msto_t1=-0.29560944
  REAL(DP), PARAMETER :: msto_z0=0.95402521,msto_z1=0.21944863,&
       msto_z2=0.070565820
  INTEGER, PARAMETER :: nimfnp5=5
  !mass boundaries for non-para IMF (starting at imflo, and ending at imfhi)
  !REAL(DP), DIMENSION(nimfnp5) :: mbin_nimf = (/0.2,0.4,0.6,0.8,1.0/)
  REAL(DP), DIMENSION(nimfnp+1) :: &
       mbin_nimf9 = (/0.08,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0/)
  REAL(DP), DIMENSION(nimfnp) :: &
       corr_salp_flat = (/2.12,1.44,1.08,0.82,0.61,0.44,0.30,0.17,0.054/)

  !----------Setup a common block of arrays and vars-------------!

  !length of input data
  INTEGER :: datmax=0
  !index in lam array at 7000A
  INTEGER :: lam7=1
  !indices for the fiducial IMF (set in setup.f90)
  INTEGER :: imfr1=1,imfr2=1,imfr3=1

  !common array for filters
  REAL(DP), DIMENSION(nl,nfil)   :: filters=0.0
  !common array for wavelength intervals
  REAL(DP), DIMENSION(nlint_max) :: l1=0.,l2=0.
  !arrays containing the upper and lower prior limits
  REAL(DP), DIMENSION(npar) :: prloarr=0.,prhiarr=0.

  !array for the template error function
  REAL(DP), DIMENSION(nl) :: temperrfcn=1.0

  !array for wavelengths of emission lines
  REAL(DP), DIMENSION(neml) :: emlines=0.

  !variables used in velbroad.f90 routine
  REAL(DP) :: dlstep=0.
  REAL(DP), DIMENSION(nl) :: lnlam=0.

  !variable pointing to the specfit home dir
  !set this environment variable in your .cshrc file
  CHARACTER(250) :: ALF_HOME=''

  !arrays holding the sky emission lines
  INTEGER, PARAMETER :: nskylines = 39324
  REAL(DP), DIMENSION(nskylines) :: lsky,fsky

  !array of index definitions
  REAL(DP), DIMENSION(7,nindx) :: indxdef=0.
  INTEGER, DIMENSION(nindx) :: indx2fit=0

  !---------------------Physical Constants-----------------------!
  !---------------in cgs units where applicable------------------!

  !pi
  REAL(DP), PARAMETER :: mypi   = 3.14159265
  !speed of light (cm/s)
  REAL(DP), PARAMETER :: clight = 2.9979E10
  !Solar mass in grams
  REAL(DP), PARAMETER :: msun   = 1.989E33
  !Solar luminosity in erg/s
  REAL(DP), PARAMETER :: lsun   = 3.839E33
  !cm in a pc
  REAL(DP), PARAMETER :: pc2cm  = 3.08568E18

  !define small and large numbers
  REAL(DP), PARAMETER :: huge_number = 1E33
  REAL(DP), PARAMETER :: tiny_number = 1E-33
  
  !-------------------Define TYPE structures---------------------!
  
  !structure for the set of parameters necessary to generate a model
  TYPE PARAMS
     REAL(DP) :: velz=0.0,sigma=0.0,logage=1.0,zh=0.0,feh=0.0,ah=0.0,&
          nhe=0.0,ch=0.0,nh=0.0,nah=0.0,mgh=0.0,sih=0.0,kh=0.0,&
          cah=0.0,tih=0.0,vh=0.0,crh=0.0,mnh=0.0,coh=0.0,nih=0.0,&
          cuh=0.0,srh=0.0,bah=0.0,euh=0.0,teff=0.0,imf1=1.3,imf2=2.3,&
          logfy=-4.0,sigma2=0.0,velz2=0.0,logm7g=-4.0,hotteff=20.0,&
          loghot=-4.0,fy_logage=0.3,logtrans=-4.0,logemline_h=-4.0,&
          logemline_oiii=-4.0,logemline_sii=-4.0,logemline_ni=-4.0,&
          logemline_nii=-4.0,jitter=1.0,imf3=2.0,logsky=-4.0,imf4=0.0
     REAL(DP) :: chi2=huge_number
  END TYPE PARAMS
  
  !structure for the models
  TYPE SSP
     REAL(DP), DIMENSION(nl) :: lam,m7g
     REAL(DP), DIMENSION(nl,nage_rfcn,nzmet) :: solar,nap,nam,cap,cam,&
          fep,fem,cp,cm,ap,np,nm,tip,tim,mgp,mgm,sip,sim,crp,mnp,bap,bam,&
          nip,cup,cop,eup,srp,kp,vp,teffp,teffm,nap6,nap9
     REAL(DP), DIMENSION(nage_rfcn)     :: logagegrid_rfcn
     REAL(DP), DIMENSION(nage)          :: logagegrid
     REAL(DP), DIMENSION(nzmet)         :: logzgrid
     REAL(DP), DIMENSION(nzmet3)        :: logzgrid2     
     REAL(DP), DIMENSION(nl,nimf,nimf,nage,nzmet)        :: logssp
     REAL(DP), DIMENSION(nl,nimf,nimf,nage,nmcut,nzmet3) :: logsspm
     REAL(DP), DIMENSION(nl,nimfnp,nage,nzmet)         :: sspnp
     REAL(DP), DIMENSION(nimf)          :: imfx1,imfx2
     REAL(DP), DIMENSION(nmcut)         :: imfx3
     REAL(DP), DIMENSION(nl,nhot)       :: hotspec
     REAL(DP), DIMENSION(nl)            :: atm_trans_h2o,atm_trans_o2
     REAL(DP), DIMENSION(nhot)          :: teffarrhot
  END TYPE SSP

  !structure for the data
  TYPE TDATA
     REAL(DP) :: lam=1E6,flx=0.0,err=0.0,wgt=0.0,ires=0.0,lam0=1E6,sky=0.0
  END TYPE TDATA

  !structure for the indices measured from the data
  TYPE IDATA
     REAL(DP) :: indx=0.0, err=99.0
  END TYPE IDATA

  !define the actual variable holding the index data
  TYPE(IDATA), DIMENSION(nindx) :: data_indx
  !define the actual SSP grid to be shared between the routines
  TYPE(SSP) :: sspgrid
  !define the object for the raw data array
  TYPE(TDATA), DIMENSION(ndat) :: data
  !define the object for the data interpolated to the model arr
  TYPE(TDATA), DIMENSION(nl)  :: idata


END MODULE ALF_VARS
