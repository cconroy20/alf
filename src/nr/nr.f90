MODULE nr

  INTERFACE
     FUNCTION brent(ax,bx,cx,func,tol,xmin)
       USE nrtype
       REAL(SP), INTENT(IN) :: ax,bx,cx,tol
       REAL(SP), INTENT(OUT) :: xmin
       REAL(SP) :: brent
       INTERFACE
          FUNCTION func(x)
            USE nrtype
            REAL(SP), INTENT(IN) :: x
            REAL(SP) :: func
          END FUNCTION func
       END INTERFACE
     END FUNCTION brent
  END INTERFACE
  
  INTERFACE
     SUBROUTINE covsrt(covar,maska)
       USE nrtype
       REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: covar
       LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: maska
     END SUBROUTINE covsrt
  END INTERFACE
  
  INTERFACE gasdev
     SUBROUTINE gasdev_s(harvest)
       USE nrtype
       REAL(SP), INTENT(OUT) :: harvest
     END SUBROUTINE gasdev_s
     !BL
     SUBROUTINE gasdev_v(harvest)
       USE nrtype
       REAL(SP), DIMENSION(:), INTENT(OUT) :: harvest
     END SUBROUTINE gasdev_v
  END INTERFACE gasdev
 
  INTERFACE
     SUBROUTINE gaussj(a,b)
       USE nrtype
       REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a,b
     END SUBROUTINE gaussj
  END INTERFACE

  INTERFACE
     SUBROUTINE lfit(x,y,sig,a,maska,covar,chisq,funcs)
       USE nrtype
       IMPLICIT NONE
       REAL(SP), DIMENSION(:), INTENT(IN) :: x,y,sig
       REAL(SP), DIMENSION(:), INTENT(INOUT) :: a
       LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: maska
       REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: covar
       REAL(SP), INTENT(OUT) :: chisq
       INTERFACE
          SUBROUTINE funcs(x,arr)
            USE nrtype
            IMPLICIT NONE
            REAL(SP),INTENT(IN) :: x
            REAL(SP), DIMENSION(:), INTENT(OUT) :: arr
          END SUBROUTINE funcs
       END INTERFACE
     END SUBROUTINE lfit
  END INTERFACE

  INTERFACE
     SUBROUTINE linmin(p,xi,fret)
       USE nrtype
       REAL(SP), INTENT(OUT) :: fret
       REAL(SP), DIMENSION(:), TARGET, INTENT(INOUT) :: p,xi
     END SUBROUTINE linmin
  END INTERFACE

  INTERFACE
     FUNCTION locate(xx,x)
       USE nrtype
       REAL(SP), DIMENSION(:), INTENT(IN) :: xx
       REAL(SP), INTENT(IN) :: x
       INTEGER(I4B) :: locate
     END FUNCTION locate
  END INTERFACE

  INTERFACE
     SUBROUTINE mnbrak(ax,bx,cx,fa,fb,fc,func)
       USE nrtype
       REAL(SP), INTENT(INOUT) :: ax,bx
       REAL(SP), INTENT(OUT) :: cx,fa,fb,fc
       INTERFACE
          FUNCTION func(x)
            USE nrtype
            REAL(SP), INTENT(IN) :: x
            REAL(SP) :: func
          END FUNCTION func
       END INTERFACE
     END SUBROUTINE mnbrak
  END INTERFACE

  INTERFACE
     SUBROUTINE powell(p,xi,ftol,iter,fret)
       USE nrtype
       REAL(SP), DIMENSION(:), INTENT(INOUT) :: p
       REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: xi
       INTEGER(I4B), INTENT(OUT) :: iter
       REAL(SP), INTENT(IN) :: ftol
       REAL(SP), INTENT(OUT) :: fret
     END SUBROUTINE powell
  END INTERFACE

  INTERFACE ran1
     SUBROUTINE ran1_s(harvest)
       USE nrtype
       REAL(SP), INTENT(OUT) :: harvest
     END SUBROUTINE ran1_s
     !BL
     SUBROUTINE ran1_v(harvest)
       USE nrtype
       REAL(SP), DIMENSION(:), INTENT(OUT) :: harvest
     END SUBROUTINE ran1_v
  END INTERFACE ran1

  INTERFACE
     SUBROUTINE sort(arr)
       USE nrtype
       REAL(SP), DIMENSION(:), INTENT(INOUT) :: arr
     END SUBROUTINE sort
  END INTERFACE
  
END MODULE nr
