MODULE ALF_UTILS

  INTERFACE
     SUBROUTINE ADD_RESPONSE(spec,pos,range,dr,vr,solar,plus,minus)
       USE alf_vars
       INTEGER, INTENT(in) :: vr
       REAL(SP), INTENT(in) :: range
       REAL(DP), INTENT(in) :: pos,dr
       REAL(DP), DIMENSION(nl), INTENT(inout) :: spec
       REAL(DP), DIMENSION(nl,nage_rfcn), INTENT(in) :: plus,solar
       REAL(DP), DIMENSION(nl,nage_rfcn), INTENT(in), OPTIONAL :: minus
     END SUBROUTINE ADD_RESPONSE
  END INTERFACE

  INTERFACE
     SUBROUTINE CONTNORMSPEC(lam,flx,err,il1,il2,flxout,coeff)
       USE alf_vars
       REAL(DP), DIMENSION(:), INTENT(in) :: lam,flx,err
       REAL(DP), INTENT(in) :: il1,il2
       REAL(DP), DIMENSION(:), INTENT(inout) :: flxout
       REAL(DP), DIMENSION(ncoeff), OPTIONAL :: coeff
     END SUBROUTINE CONTNORMSPEC
  END INTERFACE

  INTERFACE
     SUBROUTINE EMCEE_ADVANCE(ndim,nwalkers,a,pin,lpin,&
          pout,lpout,accept)
       USE alf_vars
       INTEGER, INTENT(in) :: ndim, nwalkers
       REAL(DP), INTENT(in) :: a
       REAL(DP), INTENT(in), DIMENSION(ndim,nwalkers) :: pin
       REAL(DP), INTENT(in), DIMENSION(nwalkers) :: lpin  
       REAL(DP), INTENT(out), DIMENSION(ndim,nwalkers) :: pout
       REAL(DP), INTENT(out), DIMENSION(nwalkers) :: lpout
       INTEGER, INTENT(out), DIMENSION(nwalkers) :: accept
     END SUBROUTINE EMCEE_ADVANCE
  END INTERFACE

  INTERFACE
     SUBROUTINE EMCEE_ADVANCE_MPI(ndim,nwalkers,a,pin,lpin,&
          pout,lpout,accept,nworkers)
       USE alf_vars
       INTEGER, INTENT(in) :: ndim, nwalkers,nworkers
       REAL(DP), INTENT(in) :: a
       REAL(DP), INTENT(in), DIMENSION(ndim,nwalkers) :: pin
       REAL(DP), INTENT(in), DIMENSION(nwalkers) :: lpin  
       REAL(DP), INTENT(out), DIMENSION(ndim,nwalkers) :: pout
       REAL(DP), INTENT(out), DIMENSION(nwalkers) :: lpout
       INTEGER, INTENT(out), DIMENSION(nwalkers) :: accept
     END SUBROUTINE EMCEE_ADVANCE_MPI
  END INTERFACE

  INTERFACE
     SUBROUTINE FREE_WORKERS(nworkers)
       INTEGER, INTENT(in) :: nworkers
     END SUBROUTINE FREE_WORKERS
  END INTERFACE

  INTERFACE
     FUNCTION FUNC(nposarr,spec,funit)
       USE alf_vars
       REAL(DP), DIMENSION(:), INTENT(in) :: nposarr
       REAL(DP), DIMENSION(nl), OPTIONAL :: spec
       INTEGER, INTENT(in), OPTIONAL :: funit
       REAL(DP) :: func
     END FUNCTION FUNC
  END INTERFACE

  INTERFACE
     SUBROUTINE FUNCTION_PARALLEL_MAP(ndim, nk, nworkers, pos, lnpout)
       USE alf_vars
       INTEGER, INTENT(in) :: ndim, nk, nworkers
       REAL(DP), INTENT(in), DIMENSION(ndim,nk) :: pos
       REAL(DP), INTENT(out), DIMENSION(nk) :: lnpout
    END SUBROUTINE FUNCTION_PARALLEL_MAP
  END INTERFACE

  INTERFACE
     SUBROUTINE GETM2L(msto,lam,spec,pos,m2l,mw)
       USE alf_vars
       REAL(DP), DIMENSION(nl), INTENT(in) :: lam,spec
       REAL(DP), INTENT(in) :: msto
       TYPE(PARAMS), INTENT(in)   :: pos
       REAL(DP), DIMENSION(nfil), INTENT(out) :: m2l
       INTEGER, OPTIONAL :: mw
     END SUBROUTINE GETM2L
  END INTERFACE

  INTERFACE
     FUNCTION GETMASS(mto,imf1,imf2,imf3)
       USE alf_vars
       REAL(DP), INTENT(in) :: mto,imf1,imf2,imf3
       REAL(DP) :: getmass
     END FUNCTION GETMASS
  END INTERFACE

  INTERFACE
     SUBROUTINE GETMODEL(pos,spec,mw)
       USE alf_vars
       TYPE(PARAMS), INTENT(in) :: pos
       REAL(DP), DIMENSION(nl), INTENT(out) :: spec
       INTEGER, OPTIONAL :: mw
     END SUBROUTINE GETMODEL
  END INTERFACE

  INTERFACE
     FUNCTION GETVELZ()
       USE alf_vars
       REAL(DP) :: getvelz
     END FUNCTION GETVELZ
  END INTERFACE

  INTERFACE
     FUNCTION IMF(mass)
       USE alf_vars
       REAL(DP), DIMENSION(:), INTENT(in) :: mass
       REAL(DP), DIMENSION(size(mass)) :: imf
     END FUNCTION IMF
  END INTERFACE 

  INTERFACE
     FUNCTION LINTERP(xin,yin,xout)
       USE alf_vars
       REAL(DP), DIMENSION(:), INTENT(in) :: xin,yin
       REAL(DP), INTENT(in), DIMENSION(:)  :: xout
       REAL(DP), DIMENSION(SIZE(xout)) :: linterp
     END FUNCTION LINTERP
  END INTERFACE

  INTERFACE
     SUBROUTINE LINTERP3(xin,yin1,yin2,yin3,xout,yout1,yout2,yout3)
       USE alf_vars
       REAL(DP), DIMENSION(:), INTENT(in) :: xin,yin1,yin2,yin3,xout
       REAL(DP), DIMENSION(:), INTENT(inout) :: yout1,yout2,yout3
     END SUBROUTINE LINTERP3
  END INTERFACE

  INTERFACE
     SUBROUTINE MASKEMLINES(zred,sigma)
       USE alf_vars
       REAL(DP), INTENT(in) :: zred,sigma
     END SUBROUTINE MASKEMLINES
  END INTERFACE

  INTERFACE
     FUNCTION MYRAN()
       USE alf_vars
       USE nr, ONLY : ran1
       REAL(DP) :: myran
     END FUNCTION MYRAN
  END INTERFACE

  INTERFACE
     SUBROUTINE READ_DATA(file)
       USE alf_vars
       CHARACTER(50), INTENT(in)  :: file
     END SUBROUTINE READ_DATA
  END INTERFACE

  INTERFACE
     SUBROUTINE SET_PINIT_PRIORS(pos,prlo,prhi,velz)
       USE alf_vars
       TYPE(PARAMS), INTENT(inout) :: pos,prlo,prhi
       REAL(DP), OPTIONAL :: velz
     END SUBROUTINE SET_PINIT_PRIORS
  END INTERFACE

  INTERFACE
     SUBROUTINE SETUP()
       USE alf_vars
     END SUBROUTINE SETUP
  END INTERFACE

  INTERFACE
     SUBROUTINE STR2ARR(switch,str,arr)
       USE alf_vars
       TYPE(PARAMS), INTENT(inout) :: str
       REAL(DP), DIMENSION(npar), INTENT(inout) :: arr
       INTEGER, INTENT(in) :: switch
     END SUBROUTINE STR2ARR
  END INTERFACE

  INTERFACE
     FUNCTION TSUM(xin,yin)
       USE alf_vars       
       REAL(DP), DIMENSION(:), INTENT(in) :: xin,yin
       REAL(DP) :: tsum
     END FUNCTION TSUM
  END INTERFACE

  INTERFACE
     SUBROUTINE UPDATE_RUNTOT(runtot,inarr,m2l,m2lmw)
       USE alf_vars
       REAL(DP), INTENT(inout), DIMENSION(3,npar+2*nfil) :: runtot
       REAL(DP), INTENT(in), DIMENSION(nfil) :: m2l,m2lmw
       REAL(DP), INTENT(in), DIMENSION(npar) :: inarr
     END SUBROUTINE UPDATE_RUNTOT
  END INTERFACE

  INTERFACE
     SUBROUTINE VELBROAD(lambda,spec,sigma,minl,maxl,ires)
       USE alf_vars
       REAL(DP), INTENT(in), DIMENSION(:) :: lambda
       REAL(DP), INTENT(inout), DIMENSION(:) :: spec
       REAL(DP), INTENT(in) :: sigma,minl,maxl
       REAL(DP), INTENT(in), DIMENSION(:), OPTIONAL :: ires
     END SUBROUTINE VELBROAD
  END INTERFACE


END MODULE ALF_UTILS
