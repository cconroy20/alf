MODULE SFUTILS

  INTERFACE
     SUBROUTINE CONTNORMSPEC(lam,flx,err,il1,il2,flxout)
       USE sfvars
       REAL(DP), DIMENSION(nl), INTENT(in) :: lam,flx,err
       REAL(DP), INTENT(in) :: il1,il2
       REAL(DP), DIMENSION(nl), INTENT(inout) :: flxout
     END SUBROUTINE CONTNORMSPEC
  END INTERFACE

  INTERFACE
     FUNCTION FUNC(nposarr,spec)
       USE sfvars
       REAL(DP), DIMENSION(:), INTENT(in) :: nposarr
       REAL(DP), DIMENSION(:), OPTIONAL :: spec
       REAL(DP) :: func
     END FUNCTION FUNC
  END INTERFACE

  INTERFACE
     SUBROUTINE GETM2L(lam,spec,pos,m2l,mw)
       USE sfvars
       REAL(DP), DIMENSION(nl), INTENT(in) :: lam,spec
       TYPE(PARAMS), INTENT(in)   :: pos
       REAL(DP), DIMENSION(nfil), INTENT(out) :: m2l
       INTEGER, OPTIONAL :: mw
     END SUBROUTINE GETM2L
  END INTERFACE

  INTERFACE
     FUNCTION GETMASS(mto,imf1,imf2,imf3)
       USE sfvars
       REAL(DP), INTENT(in) :: mto
       REAL(DP), INTENT(in) :: imf1,imf2,imf3
       REAL(DP) :: getmass
     END FUNCTION GETMASS
  END INTERFACE

  INTERFACE
     SUBROUTINE GETMODEL(pos,spec,mw)
       USE sfvars
       TYPE(PARAMS), INTENT(in) :: pos
       REAL(DP), DIMENSION(nl), INTENT(out) :: spec
       INTEGER, OPTIONAL :: mw
     END SUBROUTINE GETMODEL
  END INTERFACE

  INTERFACE
     FUNCTION GETVELZ()
       USE sfvars
       REAL(DP) :: getvelz
     END FUNCTION GETVELZ
  END INTERFACE

  INTERFACE
     FUNCTION IMF(mass)
       USE sfvars
       REAL(DP), DIMENSION(:), INTENT(in) :: mass
       REAL(DP), DIMENSION(size(mass)) :: imf
     END FUNCTION IMF
  END INTERFACE 

  INTERFACE
     FUNCTION LINTERP(xin,yin,xout)
       USE sfvars
       REAL(DP), DIMENSION(:), INTENT(in) :: xin,yin
       REAL(DP), INTENT(in), DIMENSION(:)  :: xout
       REAL(DP), DIMENSION(SIZE(xout)) :: linterp
     END FUNCTION LINTERP
  END INTERFACE

  INTERFACE
     SUBROUTINE MASKEMLINES(zred,sigma)
       USE sfvars
       REAL(DP), INTENT(in) :: zred,sigma
     END SUBROUTINE MASKEMLINES
  END INTERFACE

  INTERFACE
     FUNCTION MYRAN()
       USE sfvars
       USE nr, ONLY : ran1
       REAL(DP) :: myran
     END FUNCTION MYRAN
  END INTERFACE

  INTERFACE
     SUBROUTINE READ_DATA(file)
       USE sfvars
       CHARACTER(50), INTENT(in)  :: file
     END SUBROUTINE READ_DATA
  END INTERFACE

  INTERFACE
     SUBROUTINE SET_WAVE_LIMITS(file)
       USE sfvars
       CHARACTER(50), INTENT(in) :: file
     END SUBROUTINE SET_WAVE_LIMITS
  END INTERFACE

  INTERFACE
     SUBROUTINE SETUP_PARAMS(pos,dstep,prlo,prhi,velz)
       USE sfvars
       TYPE(PARAMS), INTENT(inout) :: pos,prlo,prhi,dstep
       REAL(DP), OPTIONAL :: velz
     END SUBROUTINE SETUP_PARAMS
  END INTERFACE

  INTERFACE
     SUBROUTINE SFSETUP()
       USE sfvars
     END SUBROUTINE SFSETUP
  END INTERFACE

  INTERFACE
     SUBROUTINE STR2ARR(switch,str,arr)
       USE sfvars
       TYPE(PARAMS), INTENT(inout) :: str
       REAL(DP), DIMENSION(npar), INTENT(inout) :: arr
       INTEGER, INTENT(in) :: switch
     END SUBROUTINE STR2ARR
  END INTERFACE

  INTERFACE
     FUNCTION TSUM(xin,yin)
       USE sfvars       
       REAL(DP), DIMENSION(:), INTENT(in) :: xin,yin
       REAL(DP) :: tsum
     END FUNCTION TSUM
  END INTERFACE

  INTERFACE
     SUBROUTINE VELBROAD(lambda,spec,sigma,minl,maxl)
       USE sfvars
       REAL(DP), INTENT(in), DIMENSION(:) :: lambda
       REAL(DP), INTENT(inout), DIMENSION(:) :: spec
       REAL(DP), INTENT(in) :: sigma,minl,maxl
     END SUBROUTINE VELBROAD
  END INTERFACE


END MODULE SFUTILS
