FUNCTION FUNC(nposarr,spec,funit)

  !routine to get a new model and compute chi^2.  Optionally,
  !the model spectrum is returned (spec).  The model priors
  !are computed in this routine.

  USE alf_vars; USE nr, ONLY : locate
  USE alf_utils, ONLY : linterp3,contnormspec,getmass,&
       str2arr,getmodel
  IMPLICIT NONE

  REAL(DP), DIMENSION(:), INTENT(inout) :: nposarr
  REAL(DP), DIMENSION(nl), OPTIONAL :: spec
  INTEGER, INTENT(in), OPTIONAL :: funit
  REAL(DP) :: func,pr,tchi2,ml,tl1,tl2
  REAL(DP), DIMENSION(nl)   :: mspec,mflx,dflx,poly
  REAL(DP), DIMENSION(npar) :: tposarr=0.0
  REAL(DP), DIMENSION(ncoeff) :: tcoeff
  INTEGER  :: i,i1,i2,j,npow,tpow
  TYPE(PARAMS)   :: npos

  !------------------------------------------------------!

  func = 0.0
  tpow = 0

  IF (SIZE(nposarr).LT.npar) THEN
     tposarr(1:5) = nposarr(1:5)
  ELSE
     tposarr = nposarr
  ENDIF

  CALL STR2ARR(2,npos,tposarr) !arr->str

  !compute priors (don't count all the priors if fitting
  !in (super) simple mode or in powell fitting mode)
  pr = 1.0
  DO i=1,npar
     IF (i.GT.npowell.AND.(powell_fitting.EQ.1.OR.fit_type.EQ.2)) CYCLE
     IF (fit_type.EQ.1.AND.i.GT.nparsimp) CYCLE
     IF (nposarr(i).GT.prhiarr(i)) &
          pr = pr*EXP(-(nposarr(i)-prhiarr(i))**2/2/0.001)
     IF (nposarr(i).LT.prloarr(i)) &
          pr = pr*EXP(-(nposarr(i)-prloarr(i))**2/2/0.001)
  ENDDO

  !only compute the model and chi2 if the priors are >0.0
  IF (pr.GT.tiny_number) THEN

     !get a new model spectrum
     CALL GETMODEL(npos,mspec)
     
     IF (PRESENT(spec)) THEN
        spec = mspec
     ENDIF

     !de-redshift the data and interpolate to model wavelength array
     data%lam0 = data%lam / (1+npos%velz/clight*1E5)
     CALL LINTERP3(data(1:datmax)%lam0,data(1:datmax)%flx,&
          data(1:datmax)%err,data(1:datmax)%wgt,&
          sspgrid%lam(1:nl_fit),idata(1:nl_fit)%flx,&
          idata(1:nl_fit)%err,idata(1:nl_fit)%wgt)

     !compute chi2, looping over wavelength intervals
     DO i=1,nlint
        
        tl1 = MAX(l1(i),data(1)%lam0)
        tl2 = MIN(l2(i),data(datmax)%lam0)
        !if wavelength interval falls completely outside 
        !of the range of the data, then skip
        IF (tl1.GE.tl2) CYCLE
        
        i1 = MIN(MAX(locate(sspgrid%lam,tl1),1),nl-1)
        i2 = MIN(MAX(locate(sspgrid%lam,tl2),2),nl)
        ml = (tl1+tl2)/2.0
        
        !fit a polynomial to the ratio of model and data
        IF (fitpoly.EQ.1) THEN
           CALL CONTNORMSPEC(sspgrid%lam,idata%flx/mspec,idata%err,&
                tl1,tl2,mflx,coeff=tcoeff)
           poly = 0.0
           npow = MIN(NINT((tl2-tl1)/100.0),npolymax)
           DO j=1,npow+1 
              poly = poly + tcoeff(j)*(sspgrid%lam-ml)**(j-1)
           ENDDO
           mflx  = mspec * poly
           dflx  = idata%flx
           tchi2 = SUM((dflx(i1:i2)-mflx(i1:i2))**2/idata(i1:i2)%err**2)
        ELSE
           CALL CONTNORMSPEC(sspgrid%lam,idata%flx,idata%err,tl1,tl2,dflx)
           CALL CONTNORMSPEC(sspgrid%lam,mspec,idata%wgt*SQRT(mspec),&
                tl1,tl2,mflx)
           tchi2 = SUM(idata(i1:i2)%flx**2/idata(i1:i2)%err**2*&
                (dflx(i1:i2)-mflx(i1:i2))**2)
        ENDIF
        
        !error checking
        IF (isnan(tchi2)) THEN
           WRITE(*,'(" FUNC ERROR: chi2 returned a NaN")') 
           WRITE(*,'(" error occured at wavelength interval: ",I1)') i
           WRITE(*,'("data:",2000F14.2)') dflx(i1:i2)
           WRITE(*,*)
           WRITE(*,'("model:",2000F14.2)') mspec(i1:i2)*1E3
           WRITE(*,*)
           WRITE(*,'("errors:",2000F14.2)') idata(i1:i2)%err
           STOP
        ENDIF
        
        func  = func + tchi2
        
        IF (PRESENT(funit)) THEN
           !write final results to screen and file
           WRITE(*,'(2x,F5.2,"um-",F5.2,"um:"," rms:",F5.1,"%, ","Chi2/dof:",F5.1)') &
                tl1/1E4,tl2/1E4,SQRT( SUM( (dflx(i1:i2)/mflx(i1:i2)-1)**2 )/&
                (i2-i1+1) )*100, tchi2/(i2-i1)
           DO j=i1,i2
              WRITE(funit,'(F9.2,4ES12.4)') sspgrid%lam(j),mflx(j),&
                   dflx(j),idata(j)%flx/idata(j)%err,poly(j)
           ENDDO
        ENDIF

     ENDDO
  
  ENDIF

  !include priors
  IF (pr.LE.tiny_number) THEN 
     func = huge_number 
  ELSE 
     func = func -2*LOG(pr)
  ENDIF

  !for testing purposes
  IF (1.EQ.0) THEN
     IF (powell_fitting.EQ.1) THEN
        WRITE(*,'(2ES10.3,2F12.2,99F7.3)') &
             func,pr,npos%velz,npos%sigma,10**npos%logage,npos%feh
     ELSE
        WRITE(*,'(2ES10.3,2F12.2,99F7.3)') &
             func,pr,npos%velz,npos%sigma,10**npos%logage,npos%feh,npos%ah,&
             npos%nhe,npos%ch,npos%nh,npos%nah,npos%mgh,npos%sih,npos%kh,&
             npos%cah,npos%tih!,npos%vh,npos%crh,npos%mnh,npos%coh,npos%nih !,&
        !     npos%rbh,npos%srh,npos%yh,npos%zrh,npos%bah,npos%euh,npos%teff,&
        !     npos%imf1,npos%imf2,npos%logfy,npos%velz,npos%logm7g,npos%hotteff,&
        !     npos%loghot
     ENDIF
  ENDIF

END FUNCTION FUNC
