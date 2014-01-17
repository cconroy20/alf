FUNCTION FUNC(nposarr,spec)

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
  REAL(DP) :: func,pr,tchi2
  REAL(DP), DIMENSION(nl)   :: mspec,mflx,dflx
  REAL(DP), DIMENSION(ndat) :: tlam
  INTEGER  :: i,i1,i2
  TYPE(PARAMS)   :: npos
  TYPE(TDATA), DIMENSION(nl) :: idata

  !------------------------------------------------------!

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

  !compute chi2
  func = 0.0
  DO i=1,nlint
     !if wavelength interval exceeds data range, then skip
     IF (MAXVAL(tlam(1:datmax)).LT.l2(i)) CYCLE
     IF (MINVAL(tlam(1:datmax)).GT.l1(i)) CYCLE
     CALL CONTNORMSPEC(sspgrid%lam,idata%flx,idata%err,l1(i),l2(i),dflx)
     CALL CONTNORMSPEC(sspgrid%lam,mspec,idata%wgt,l1(i),l2(i),mflx)
     i1 = MIN(MAX(locate(sspgrid%lam,l1(i)),1),nl-1)
     i2 = MIN(MAX(locate(sspgrid%lam,l2(i)),2),nl)
     tchi2 = SUM(idata(i1:i2)%flx**2/idata(i1:i2)%err**2*&
          (dflx(i1:i2)-mflx(i1:i2))**2)
     func  = func + tchi2
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
     WRITE(*,'(2ES10.3,F12.2,99F7.3)') &
          func,pr,npos%sigma,npos%velz,npos%logage,npos%feh,npos%afe,&
          npos%nhe,npos%cfe,npos%nfe,npos%nafe,npos%mgfe,npos%sife,npos%kfe,&
          npos%cafe,npos%tife,npos%vfe,npos%crfe,npos%mnfe,npos%cofe,npos%nife,&
          npos%rbfe,npos%srfe,npos%yfe,npos%zrfe,npos%bafe,npos%eufe,npos%teff,&
          npos%imf1,npos%imf2,npos%logfy,npos%velz,npos%logm7g,npos%hotteff,&
          npos%loghot
  ENDIF

END FUNCTION FUNC
