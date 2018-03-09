FUNCTION FUNC(nposarr,spec,funit)

  !routine to get a new model and compute chi^2.  Optionally,
  !the model spectrum is returned (spec).  The model priors
  !are computed in this routine.

  USE alf_vars; USE nr, ONLY : locate
  USE alf_utils, ONLY : linterp3,contnormspec,getmass,&
       str2arr,getmodel,linterp,getindx,getm2l
  IMPLICIT NONE

  REAL(DP), DIMENSION(:), INTENT(inout) :: nposarr
  REAL(DP), DIMENSION(nl), OPTIONAL :: spec
  INTEGER, INTENT(in), OPTIONAL :: funit
  REAL(DP) :: func,pr,tchi2,ml,tl1,tl2,oneplusz,tmps
  REAL(DP), DIMENSION(nfil) :: mlpr,mlalf
  REAL(DP), DIMENSION(nl)   :: mspec
  REAL(DP), DIMENSION(ndat) :: mflx,poly,zmspec,terr
  REAL(DP), DIMENSION(npar) :: tposarr=0.0
  REAL(DP), DIMENSION(npolymax+1) :: tcoeff
  REAL(DP), DIMENSION(nindx) :: mindx
  INTEGER  :: i,i1,i2,j,npow,tpow
  TYPE(PARAMS)   :: npos

  !------------------------------------------------------!

  func = 0.0
  tpow = 0

  !this is for Powell minimization
  IF (SIZE(nposarr).LT.npar) THEN
     !copy over the default parameters first
     CALL STR2ARR(1,npos,tposarr) !str->arr
     !only copy the first four params (velz,sigma,age,[Z/H])
     tposarr(1:npowell) = nposarr(1:npowell)
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
     IF (fit_indices.EQ.1.AND.i.LE.2) CYCLE
     IF ( (nposarr(i).GT.prhiarr(i)).OR.(nposarr(i).LT.prloarr(i)) ) pr=0.0
  ENDDO

  !regularize the non-parametric IMF
  !the IMF cannot be convex (U shaped)
  IF (imf_type.EQ.4.AND.nonpimf_regularize.EQ.1) THEN
     IF ( (npos%imf2-npos%imf1+corr_bin_weight(3)-corr_bin_weight(1)).LT.0.0.AND.&
          (npos%imf3-npos%imf2+corr_bin_weight(5)-corr_bin_weight(3)).GT.0.0 ) pr=0.0
     IF ( (npos%imf3-npos%imf2+corr_bin_weight(5)-corr_bin_weight(3)).LT.0.0.AND.&
          (npos%imf4-npos%imf3+corr_bin_weight(7)-corr_bin_weight(5)).GT.0.0 ) pr=0.0
     IF ( (npos%imf4-npos%imf3+corr_bin_weight(7)-corr_bin_weight(5)).LT.0.0.AND.&
          (0.0-npos%imf4+corr_bin_weight(9)-corr_bin_weight(7)).GT.0.0 ) pr=0.0
  ENDIF

  !only compute the model and chi2 if the priors are >0.0
  IF (pr.GT.tiny_number) THEN

     !get a new model spectrum
     CALL GETMODEL(npos,mspec)
     
     IF (PRESENT(spec)) THEN
        spec = mspec
     ENDIF
     
     !include external M/L prior 
     IF (extmlpr.EQ.1) THEN
        CALL GETM2L(sspgrid%lam,mspec,npos,mlalf)
        mlpr = linterp(mlprtab(1:nmlprtabmax,1),mlprtab(1:nmlprtabmax,2),mlalf)
     ENDIF

     
     IF (fit_indices.EQ.0) THEN

        !redshift the model and interpolate to data wavelength array
        oneplusz = (1+npos%velz/clight*1E5)
        zmspec   = 0.0
        zmspec(1:datmax) = LINTERP(sspgrid%lam(1:nl_fit)*oneplusz,&
             mspec(1:nl_fit),data(1:datmax)%lam)

        !compute chi2, looping over wavelength intervals
        DO i=1,nlint
        
           tl1 = MAX(l1(i)*oneplusz,data(1)%lam)
           tl2 = MIN(l2(i)*oneplusz,data(datmax)%lam)
           ml  = (tl1+tl2)/2.0
           !if wavelength interval falls completely outside 
           !of the range of the data, then skip
           IF (tl1.GE.tl2) CYCLE

           i1 = MIN(MAX(locate(data(1:datmax)%lam,tl1),1),datmax-1)
           i2 = MIN(MAX(locate(data(1:datmax)%lam,tl2),2),datmax)

           !fit a polynomial to the ratio of model and data
           CALL CONTNORMSPEC(data(1:datmax)%lam,&
                data(1:datmax)%flx/zmspec(1:datmax),&
                data(1:datmax)%err/zmspec(1:datmax),tl1,tl2,&
                mflx(1:datmax),coeff=tcoeff)
           poly = 0.0
           npow = MIN(NINT((tl2-tl1)/poly_dlam),npolymax)
           DO j=1,npow+1 
              poly = poly + tcoeff(j)*(data%lam-ml)**(j-1)
           ENDDO
           mflx  = zmspec * poly
           !compute chi^2
           IF (fit_type.EQ.0) THEN
              !include jitter term
              tchi2 = SUM( (data(i1:i2)%flx-mflx(i1:i2))**2 / &
                   (data(i1:i2)%err**2*npos%jitter**2+&
                   (10**npos%logsky*data(i1:i2)%sky)**2) + &
                   LOG(2*mypi*(data(i1:i2)%err**2*npos%jitter**2+&
                   (10**npos%logsky*data(i1:i2)%sky)**2)) )
           ELSE
              !no jitter in simple mode
              tchi2 = SUM( (data(i1:i2)%flx-mflx(i1:i2))**2/data(i1:i2)%err**2 )
           ENDIF
           
           !error checking
           IF (isnan(tchi2)) THEN
              WRITE(*,'(" FUNC ERROR: chi2 returned a NaN")') 
              WRITE(*,'(" error occured at wavelength interval: ",I1)') i
              WRITE(*,*) 'lam  data   err   model   poly'
              DO j=i1,i2
                 WRITE(*,'(F7.1,10ES12.3)') data(j)%lam,data(j)%flx,&
                      SQRT(data(j)%err**2*npos%jitter**2+&
                      (10**npos%logsky*data(j)%sky)**2),mflx(j),poly(j)
              ENDDO
              WRITE(*,*)
              WRITE(*,'("params:",100F14.2)') tposarr
              STOP
           ENDIF
        
           func  = func + tchi2
        
           IF (PRESENT(funit)) THEN
              !write final results to screen and file
              WRITE(*,'(2x,F5.2,"um-",F5.2,"um:"," rms:",F5.1,"%,'//&
                   ' ","Chi2/dof:",F6.1)') &
                   tl1/1E4/oneplusz,tl2/1E4/oneplusz,&
                   SQRT( SUM( (data(i1:i2)%flx/mflx(i1:i2)-1)**2 )/&
                   (i2-i1+1) )*100,tchi2/(i2-i1)
              IF (fit_type.EQ.0) THEN
                 terr = SQRT(data%err**2*npos%jitter**2+&
                      (10**npos%logsky*data%sky)**2)
              ELSE
                 terr = data%err
              ENDIF
              DO j=i1,i2
                 WRITE(funit,'(F9.2,8ES12.4)') data(j)%lam,mflx(j),&
                      data(j)%flx,data(j)%flx/terr(j),poly(j),data(j)%err
              ENDDO
           ENDIF

        ENDDO

     ELSE

        !compute indices
        CALL GETINDX(sspgrid%lam,mspec,mindx)
        !compute chi^2 for the indices
        func = SUM( (data_indx%indx-mindx)**2/data_indx%err**2 )

        IF (PRESENT(funit)) THEN
           DO j=1,nindx
              WRITE(funit,'(F8.2,3F9.4)') (indxdef(1,j)+indxdef(2,j))/2.,&
                   mindx(j),data_indx(j)%indx,data_indx(j)%err
           ENDDO
        ENDIF

     ENDIF
  
  ENDIF

  !include priors
  IF (pr.LE.tiny_number) THEN 
     func = huge_number 
  ELSE 
     func = func - 2*LOG(pr)
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
