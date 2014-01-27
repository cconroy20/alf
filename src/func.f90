FUNCTION FUNC(nposarr,spec,funit)

  !routine to get a new model and compute chi^2.  Optionally,
  !the model spectrum is returned (spec)

  !TBD:
  !solve for the polynomial that matches the data to the model
  !i.e., solve for *one* polynomial, rather than two, and 
  !put the coefficients into the MCMC just like every other parameter
  !its a more robust way to do continuum matching; see the PPXF paper
  !or one of Dan Kelson's ancient papers

  USE sfvars; USE nr, ONLY : locate
  USE sfutils, ONLY : linterp,contnormspec,getmass,str2arr,getmodel
  IMPLICIT NONE

  REAL(DP), DIMENSION(:), INTENT(inout) :: nposarr
  REAL(DP), DIMENSION(:), OPTIONAL :: spec
  INTEGER, INTENT(in), OPTIONAL :: funit
  REAL(DP) :: func,pr,tchi2,ml
  REAL(DP), DIMENSION(nl)   :: mspec,mflx,dflx,poly
  REAL(DP), DIMENSION(ndat) :: tlam
  REAL(DP), DIMENSION(ncoeff) :: tcoeff
  INTEGER  :: i,i1,i2,j,npow,tpow
  TYPE(PARAMS)   :: npos
  TYPE(TDATA), DIMENSION(nl) :: idata

  !------------------------------------------------------!

  func = 0.0
  tpow = 0

  CALL STR2ARR(2,npos,nposarr) !arr->str

  !get a new model spectrum
  CALL GETMODEL(npos,mspec)

  IF (PRESENT(spec)) THEN
     spec = mspec
  ENDIF

  !de-redshift the data and interpolate to model wave array
  tlam      = data%lam / (1+npos%velz/clight*1E5)
  idata%flx = linterp(tlam,data%flx,sspgrid%lam)
  idata%err = linterp(tlam,data%err,sspgrid%lam)
  idata%wgt = linterp(tlam,data%wgt,sspgrid%lam)

  !compute chi2, looping over wavelength intervals
  DO i=1,nlint

     !if wavelength interval exceeds data range, then skip
     IF ((tlam(datmax).LT.l2(i)).OR.(tlam(1).GT.l1(i))) CYCLE
     i1 = MIN(MAX(locate(sspgrid%lam,l1(i)),1),nl-1)
     i2 = MIN(MAX(locate(sspgrid%lam,l2(i)),2),nl)
     ml = (l1(i)+l2(i))/2.0

     IF (fitpoly.EQ.1) THEN
        CALL CONTNORMSPEC(sspgrid%lam,idata%flx/mspec,idata%err,&
             l1(i),l2(i),mflx,coeff=tcoeff)
        poly = 0.0
        npow = MIN(NINT((l2(i)-l1(i))/100.0),14)
        DO j=1,npow+1 
           poly = poly + tcoeff(j)*(sspgrid%lam-ml)**(j-1)
        ENDDO
        mflx = mspec * poly
        dflx = idata%flx
        tchi2 = SUM((dflx(i1:i2)-mflx(i1:i2))**2/idata(i1:i2)%err**2)
     ELSE
        CALL CONTNORMSPEC(sspgrid%lam,idata%flx,idata%err,l1(i),l2(i),dflx)
        CALL CONTNORMSPEC(sspgrid%lam,mspec,idata%wgt*SQRT(mspec),&
             l1(i),l2(i),mflx)
        tchi2 = SUM(idata(i1:i2)%flx**2/idata(i1:i2)%err**2*&
             (dflx(i1:i2)-mflx(i1:i2))**2)
    ENDIF

     func  = func + tchi2

     IF (PRESENT(funit)) THEN
        !write final results to screen and file
        WRITE(*,'("  rms:",F5.2,"%")') &
             SQRT( SUM( (dflx(i1:i2)/mflx(i1:i2)-1)**2 )/(i2-i1+1) )*100
        DO j=i1,i2
           WRITE(funit,'(F9.2,3ES12.4)') sspgrid%lam(j),mflx(j),&
                dflx(j),idata(j)%flx/idata(j)%err
        ENDDO
     ENDIF

  ENDDO

  
  !include priors
  pr = 1.0
  DO i=1,npar
     IF (nposarr(i).GT.prhiarr(i)) &
          pr = pr*EXP(-(nposarr(i)-prhiarr(i))**2/2/0.01)
     IF (nposarr(i).LT.prloarr(i)) &
          pr = pr*EXP(-(nposarr(i)-prloarr(i))**2/2/0.01)
  ENDDO
  IF (pr.LT.tiny_number) THEN 
     func = huge_number 
  ELSE 
     func = func -2*LOG(pr)
  ENDIF
   
  !for testing purposes
  IF (1.EQ.0) THEN
     WRITE(*,'(2ES10.3,2F12.2,99F7.3)') &
          func,pr,npos%velz,npos%sigma,npos%logage,npos%feh,npos%ah,&
          npos%nhe,npos%ch,npos%nh,npos%nah,npos%mgh,npos%sih,npos%kh,&
          npos%cah,npos%tih,npos%vh,npos%crh,npos%mnh,npos%coh,npos%nih,&
          npos%rbh,npos%srh,npos%yh,npos%zrh,npos%bah,npos%euh,npos%teff,&
          npos%imf1,npos%imf2,npos%logfy,npos%velz,npos%logm7g,npos%hotteff,&
          npos%loghot
  ENDIF

END FUNCTION FUNC
