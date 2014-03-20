MODULE SFVARS

  ! module to set up most arrays and variables 

  IMPLICIT NONE
  SAVE

  !directory for results (mcmc, bestspec, etc.)
  CHARACTER(100) :: OUTDIR='results/'

  INTEGER, PARAMETER :: DP = KIND(1.0D0)
  INTEGER, PARAMETER :: SP = KIND(1.0)

  !-----------------Define various parameters--------------------!

  !flags for the user to choose:

  !fit only a subset of the full model parameters 
  !e.g., no IMF, no nuisance parameters, no "exotic" elements
  INTEGER, PARAMETER :: fitsimple=0
  !force [Na/H]=[Mg/H]
  INTEGER, PARAMETER :: force_nah=0
  !force the IMF to be a MW IMF if set
  !this is automatically assumed if fitsimple=1
  INTEGER, PARAMETER :: mwimf=0

  !fit a polynomial to the ratio of model and data
  !if zero, then both data and model are continuum divided
  INTEGER, PARAMETER :: fitpoly=1
  !mask emission lines? (if 0, then the em lines are incl in the fit)
  INTEGER, PARAMETER :: maskem=0
  !apply template error function? (only works for SDSS stacks)
  INTEGER, PARAMETER :: apply_temperrfcn=0
  !if set, compute velocity broadening via a simple method
  !rather than the proper convolution in log_lambda space
  !don't turn this on - the "correct" version is just as fast
  INTEGER, PARAMETER :: velbroad_simple=0
  !turn on the use of age-dependent response functions
  INTEGER, PARAMETER :: use_age_dep_resp_fcns=1

  !--------------------------------------------------------!
  ! the parameters below should not be modified unless you !
  ! know what you are doing!                               !
  !--------------------------------------------------------!

  !nstart and nend allow us to use only a subset of 
  !the full wavelength array
  INTEGER, PARAMETER :: nstart = 2100  ! 0.39 um
  INTEGER, PARAMETER :: nend   = 14125 ! 2.4 um
  !number of spectral elements in SSPs
  INTEGER, PARAMETER :: nl = nend-nstart+1
  !number actually used over the range to be fit
  INTEGER :: nl_fit=nl
  !(max) number of wavelength intervals
  INTEGER, PARAMETER :: nlint_max = 10
  !actual number of wavelength intervals
  INTEGER :: nlint = 0  
  !number of emission lines to fit
  INTEGER, PARAMETER :: neml = 13
  !number of coefficients for the polynomial fitting
  INTEGER, PARAMETER :: ncoeff = 30
  !number of parameters (minus em lines and polynomial)
  INTEGER, PARAMETER :: npar1 = 36
  !number of ages in the empirical SSP grid
  INTEGER, PARAMETER :: nage = 7
  !number of ages in the response functions
  INTEGER, PARAMETER :: nage_rfcn = 5
  !number of IMF values in the SSP grid
  INTEGER, PARAMETER :: nimf = 35
  !max number of data wavelength points
  INTEGER, PARAMETER :: ndat = 1E5
  !total number of model parameters
  INTEGER, PARAMETER :: npar = npar1 + neml !+ ncoeff
  !number of filters
  INTEGER, PARAMETER :: nfil=3
  !mag of sun in r,I,K filters (AB mag)
  REAL(DP), PARAMETER, DIMENSION(3) :: magsun = (/4.64,4.52,5.14/)
  !lower and upper limits for the IMF
  REAL(DP), PARAMETER :: imflo=0.08,imfhi=100.0
  !power-law slopes for a Kroupa IMF
  REAL(DP), PARAMETER :: krpa_imf1=1.3,krpa_imf2=2.3,krpa_imf3=2.3
  !linear fit to log(age) vs. log(MS TO mass)
  REAL(DP), PARAMETER :: msto_fit0=0.290835,msto_fit1=-0.301566
  !length of input data
  INTEGER :: datmax=0

  !----------Setup a common block of arrays and vars-------------!

  !common array for filters
  REAL(DP), DIMENSION(nfil,nl) :: fil=0.0
  !common array for wavelength intervals
  REAL(DP), DIMENSION(nlint_max)   :: l1,l2
  !arrays containing the upper and lower prior limits
  REAL(DP), DIMENSION(npar) :: prloarr=0.,prhiarr=0.

  !array for the template error function
  REAL(DP), DIMENSION(nl) :: temperrfcn=1.0

  !variables used in velbroad.f90 routine
  REAL(DP) :: dlstep
  REAL(DP), DIMENSION(nl) :: lnlam

  !Kroupa-like IMF slopes
  REAL(DP), DIMENSION(3) :: imf_alpha=(/1.3,2.3,2.3/)
  REAL(DP) :: imf_flag=0

  !array of central wavelengths for emission lines
  REAL(DP), DIMENSION(neml) :: emlines=0.0

  !variable pointing to the specfit home dir
  !set this environment variable in your .cshrc file
  CHARACTER(250) :: SPECFIT_HOME=''

  !indices where x=1.3,x=2.3 in the IMF array
  INTEGER :: i13,i23

  !-------------Physical Constants---------------!
  !-------in cgs units where applicable----------!

  !pi
  REAL(DP), PARAMETER :: mypi    = 3.14159265
  !speed of light (cm/s)
  REAL(DP), PARAMETER :: clight  = 2.9979E10
  !Solar mass in grams
  REAL(DP), PARAMETER :: msun    = 1.989E33
  !Solar luminosity in erg/s
  REAL(DP), PARAMETER :: lsun    = 3.839E33
  !cm in a pc
  REAL(DP), PARAMETER :: pc2cm   = 3.08568E18
  
  !define small and large numbers
  REAL(DP), PARAMETER :: huge_number = 1E33
  REAL(DP), PARAMETER :: tiny_number = 1E-33
  
  !------------Define TYPE structures-------------!
  
  !structure for the set of parameters necessary to generate a model
  TYPE PARAMS
     REAL(DP) :: velz=0.0,sigma=0.0,logage=1.0,feh=0.0,ah=0.0,&
          nhe=0.0,ch=0.0,nh=0.0,nah=0.0,mgh=0.0,sih=0.0,kh=0.0,&
          cah=0.0,tih=0.0,vh=0.0,crh=0.0,mnh=0.0,coh=0.0,nih=0.0,&
          cuh=0.0,rbh=0.0,srh=0.0,yh=0.0,zrh=0.0,bah=0.0,euh=0.0,&
          teff=0.0,imf1=1.3,imf2=2.3,logfy=-5.0,sigma2=0.0,velz2=0.0,&
          logm7g=-5.0,hotteff=20.0,loghot=-5.0,fy_logage=0.0
     REAL(DP), DIMENSION(neml)   :: logemnorm=-5.0
     REAL(DP), DIMENSION(ncoeff) :: logcoeff=-10.0
     REAL(DP) :: chi2=huge_number
  END TYPE PARAMS
  
  !structure for the model SSPs
  TYPE SSP
     REAL(DP), DIMENSION(nl) :: lam,m7g
     REAL(DP), DIMENSION(nage_rfcn,nl) :: solar,hep,hem,nap,nam,cap,cam,&
          fep,fem,cp,cm,ap,np,nm,tip,tim,mgp,mgm,sip,sim,crp,mnp,bap,bam,&
          nip,cup,cop,eup,srp,kp,vp,yp,zp,zm,zrp,rbp,teffp,teffm,nap6,nap9
     REAL(DP), DIMENSION(nage_rfcn)   :: logagegrid_rfcn
     REAL(DP), DIMENSION(nage)    :: logagegrid
     REAL(DP), DIMENSION(nage,nl) :: logfkrpa
     REAL(DP), DIMENSION(nimf,nimf,nl) :: imf
     REAL(DP), DIMENSION(nimf) :: imfx
     REAL(DP), DIMENSION(4,nl) :: hotspec
     REAL(DP), DIMENSION(4)    :: teffarrhot
  END TYPE SSP

  !structure for the data
  TYPE TDATA
     REAL(DP) :: lam=1E6,flx=0.,err=0.,wgt=0.0
  END TYPE TDATA

  !define the actual SSP grid to be shared between the routines
  TYPE(SSP) :: sspgrid
  !define the actual variable for the raw data array
  TYPE(TDATA), DIMENSION(ndat) :: data

END MODULE SFVARS
