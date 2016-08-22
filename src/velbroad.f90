SUBROUTINE VELBROAD(lambda,spec,sigma,minl,maxl,ires)

  !routine to compute velocity broadening of an input spectrum
  !the PSF kernel has a width of m*sigma, where m=4
  !If optional input ires is present, then the spectrum will be
  !smoothed by a wavelength dependent velocity dispersion in the
  !'vebroad_simple=1' mode

  USE alf_vars; USE nr, ONLY : locate
  USE alf_utils, ONLY : linterp,tsum
  IMPLICIT NONE
  
  REAL(DP), INTENT(in), DIMENSION(:) :: lambda
  REAL(DP), INTENT(inout), DIMENSION(:) :: spec
  REAL(DP), INTENT(in) :: sigma,minl,maxl
  REAL(DP), INTENT(in), DIMENSION(:), OPTIONAL :: ires
  REAL(DP), DIMENSION(40000) :: tspec,nspec,vel,func,psf
  REAL(DP) :: xmax,fwhm,psig,sigmal
  INTEGER :: i,il,ih,m=4,grange,nn,n2

  !---------------------------------------------------------------!
  !---------------------------------------------------------------!
   
  nn = SIZE(lambda)

  IF (sigma.LE.10.0) RETURN

  IF (sigma.GE.1E4) THEN
     WRITE(*,*) "VELBROAD ERROR: sigma>1E4 km/s - you've "//&
          "probably done something wrong..."
     STOP
  ENDIF

  !compute smoothing the slightly less accurate way
  !but the only way in the case of wave-dep smoothing
  IF (velbroad_simple.EQ.1.OR.PRESENT(ires)) THEN

     RETURN

     tspec(1:nn) = spec(1:nn)
     spec(1:nn)  = 0.0

     DO i=1,nn

        IF (lambda(i).LT.minl.OR.lambda(i).GT.maxl) THEN
           spec(i)=tspec(i)
           CYCLE
        ENDIF
        
        IF (PRESENT(ires)) THEN
           sigmal = ires(i)
           IF (sigmal.LE.tiny_number) CYCLE
        ELSE
           sigmal = sigma
        ENDIF

        xmax = lambda(i)*(m*sigmal/clight*1E5+1)
        ih   = MIN(locate(lambda(1:nn),xmax),nn)
        il   = MAX(2*i-ih,1)
                
        IF (il.EQ.ih) THEN
           spec(i) = tspec(i)
        ELSE
           vel(il:ih)  = (lambda(i)/lambda(il:ih)-1)*clight/1E5
           func(il:ih) =  1/SQRT(2*mypi)/sigmal * &
                EXP(-vel(il:ih)**2/2./sigmal**2)
           !normalize the weights to integrate to unity
           func(il:ih) = func(il:ih) / TSUM(vel(il:ih),func(il:ih))
           spec(i) = TSUM(vel(il:ih),func(il:ih)*tspec(il:ih))
        ENDIF
         
     ENDDO

  !compute smoothing the correct way
  ELSE

     !fancy footwork to allow for input spectra of either length
     IF (nn.EQ.nl) THEN 
        n2 = nl_fit
     ELSE
        n2 = nl
     ENDIF

     fwhm   = sigma*2.35482/clight*1E5/dlstep
     psig   = fwhm/2.d0/SQRT(-2.d0*LOG(0.5d0)) ! equivalent sigma for kernel
     grange = FLOOR(m*psig) ! range for kernel (-range:range)

     IF (grange.GT.1) THEN 

        tspec = 0.0
        nspec = 0.0
        tspec(1:n2) = linterp(LOG(lambda(1:n2)),spec(1:n2),lnlam(1:n2))
 
        DO i=1,2*grange+1
           psf(i) = 1.d0/SQRT(2.d0*mypi)/psig*EXP(-((i-grange-1)/psig)**2/2.d0)
        ENDDO
        psf(1:2*grange+1) = psf(1:2*grange+1) / SUM(psf(1:2*grange+1))
      
        DO i=grange+1,n2-grange
           nspec(i) = SUM( psf(1:2*grange+1)*tspec(i-grange:i+grange) )
        ENDDO
        nspec(n2-grange+1:n2) = spec(n2-grange+1:n2)

        !interpolate back to the main array
        spec(1:n2) = linterp(EXP(lnlam(1:n2)),nspec(1:n2),lambda(1:n2))
 
     ENDIF

  ENDIF

END SUBROUTINE VELBROAD
