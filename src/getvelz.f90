FUNCTION GETVELZ()

  !function to estimate the recession velocity
  !this routine is used to get a first-guess at the velocity
  !so that the subsequent Powell minimization (in alf)
  !coverges faster.  Uses  the first two wavelength segments

  USE alf_vars; USE nr, ONLY : locate
  USE alf_utils, ONLY : linterp,contnormspec,velbroad
  IMPLICIT NONE

  REAL(DP) :: getvelz, chi2,lo,hi
  !delta_chi2 tolerance to accept a non-zero redshift
  REAL(DP), PARAMETER :: max_dvelz   = 5E3
  REAL(DP), PARAMETER :: delchi2_tol = 0.5   !0.5, 0.2
  REAL(DP), PARAMETER :: max_zred    = 0.18  !0.18, 0.03
  INTEGER, PARAMETER  :: nv=5000
  INTEGER :: i,i1,i2,j,ni,ng,k
  REAL(DP), DIMENSION(nl) :: mflx,dflx
  REAL(DP), DIMENSION(nv) :: tvz,tchi2,tvza
  TYPE(TDATA), DIMENSION(nl) :: iidata

  !------------------------------------------------------!
 
  chi2    = huge_number
  getvelz = 0.0
  tchi2   = 0.0

  !use ni=2 wavelength segments unless only 1 segment exists
  IF (nlint.GE.2) THEN
     ni = 2
  ELSE
     ni = 1
  ENDIF

  DO i=1,nv

     tvz(i) = REAL(i)/REAL(nv)*(max_zred*3E5+1E3)-1E3

     !de-redshift the data and interpolate to model wave array
     !NB: this is the old way of doing things, compare with func.f90
     data%lam0 = data%lam / (1+tvz(i)/clight*1E5)
     iidata(1:nl_fit)%flx = linterp(data(1:datmax)%lam0,&
          data(1:datmax)%flx,sspgrid%lam(1:nl_fit))
     iidata(1:nl_fit)%err = linterp(data(1:datmax)%lam0,&
          data(1:datmax)%err,sspgrid%lam(1:nl_fit))
     
     !only use the first ni wavelength segments
     DO j=1,ni

        lo = MAX(l1(j),data(1)%lam0)+50
        hi = MIN(l2(j),data(datmax)%lam0)-50
        !dont use the near-IR in redshift fitting
        IF (lo.GT.9500.) CYCLE
        IF (lo.GE.hi) THEN
           WRITE(*,*) 'GETVELZ ERROR:, lo>hi'
           STOP
        ENDIF

        !NB: this is the old way of doing things, compare with func.f90
        CALL CONTNORMSPEC(sspgrid%lam,iidata%flx,iidata%err,lo,hi,dflx)
        !use a 5 Gyr Zsol SSP
        CALL CONTNORMSPEC(sspgrid%lam,10**sspgrid%logssp(:,imfr1,imfr2,3,nzmet-1),&
             SQRT(10**sspgrid%logssp(:,imfr1,imfr2,3,nzmet-1)),lo,hi,mflx)

      !  CALL VELBROAD(sspgrid%lam,mflx,300.d0,lo,hi)

        i1 = MIN(MAX(locate(sspgrid%lam,lo),1),nl_fit-1)
        i2 = MIN(MAX(locate(sspgrid%lam,hi),2),nl_fit)

        !only count pixels with non-zero weights
        ng = 0
        DO k=i1,i2
           IF (iidata(k)%err.LT.(huge_number/2.)) ng=ng+1
        ENDDO

        tchi2(i) = SUM(iidata(i1:i2)%flx**2/iidata(i1:i2)%err**2*&
             (dflx(i1:i2)-mflx(i1:i2))**2) / ng
    !    tchi2(i) = SUM(iidata(i1:i2)%flx**2/iidata(i1:i2)%err**2*&
    !         (dflx(i1:i2)-mflx(i1:i2))**2) / (i2-i1)

        IF (tchi2(i).LT.chi2) THEN
           chi2    = tchi2(i)
           getvelz = tvz(i)
        ENDIF

     ENDDO

  ENDDO

  !test to see if the solution is good
  !we take all points with delta(chi2/dof)<delchi2_tol and
  !ask how large is the range in velocities.
  !If the range in velz is >max_dvelz then we've failed
  tchi2 = tchi2 - MINVAL(tchi2)
  tvza  = getvelz
  DO i=1,nv
     IF (tchi2(i).LT.delchi2_tol) tvza(i)=tvz(i)
  ENDDO
  IF ((MAXVAL(tvza)-MINVAL(tvza)).GT.max_dvelz) THEN
     WRITE(*,'("   Failed to find a redshift solution, setting velz=0.0")')
     WRITE(*,'("    delta(velz|chi2<0.5) = ",F8.2)') MAXVAL(tvza)-MINVAL(tvza)

     DO i=1,nv
        write(99,*) tvz(i),tchi2(i)
     ENDDO

     getvelz = 0.0
  ENDIF



END FUNCTION GETVELZ
