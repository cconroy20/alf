SUBROUTINE VELBROAD(lambda,spec,sigma,minl,maxl)

  !routine to compute velocity broadening of an input spectrum
  !the PSF kernel has a width of m*sigma, where m=4

  USE sfvars; USE nr, ONLY : locate; USE sfutils, ONLY : linterp,tsum
  IMPLICIT NONE
  
  REAL(DP), INTENT(in), DIMENSION(nl) :: lambda
  REAL(DP), INTENT(inout), DIMENSION(nl) :: spec
  REAL(DP), INTENT(in) :: sigma,minl,maxl
  REAL(DP), DIMENSION(nl) :: tspec,nspec,vel,func,gauss,psf
  REAL(DP) :: xmax,xmin,fwhm,psig
  INTEGER :: i,j,il,ih,m=4,grange

  !---------------------------------------------------------------!
  !---------------------------------------------------------------!
  
  IF (sigma.EQ.0) RETURN

  !compute smoothing the fast (and slightly less accurate) way
  IF (velbroad_simple.EQ.1) THEN

     tspec = spec
     spec  = 0.0

     DO i=1,nl
        
        IF (lambda(i).LT.minl.OR.lambda(i).GT.maxl) THEN
           spec(i)=tspec(i)
           CYCLE
        ENDIF
        
        xmax = lambda(i)*(m*sigma/clight*1E5+1)
        ih   = MIN(locate(lambda(1:nl),xmax),nl)
        il   = MAX(2*i-ih,1)
                
        IF (il.EQ.ih) THEN
           spec(i) = tspec(i)
        ELSE
           vel(il:ih)  = (lambda(i)/lambda(il:ih)-1)*clight/1E5
           func(il:ih) =  1/SQRT(2*mypi)/sigma * &
                EXP(-vel(il:ih)**2/2./sigma**2)
           !normalize the weights to integrate to unity
           func(il:ih) = func(il:ih) / TSUM(vel(il:ih),func(il:ih))
           spec(i) = TSUM(vel(il:ih),func(il:ih)*tspec(il:ih))
        ENDIF
         
     ENDDO

  !compute smoothing the correct way
  !actually, this way is faster than the one above!
  ELSE

     tspec = linterp(LOG(lambda(1:nl)),spec(1:nl),lnlam)
     
     fwhm   = sigma*2.35482/clight*1E5/dlstep
     psig   = fwhm/2.d0/SQRT(-2.d0*LOG(0.5d0)) ! equivalent sigma for kernel
     grange = FLOOR(m*psig)	               ! range for kernel (-range:range)
     
     DO i=1,2*grange+1
        psf(i) = 1.d0/SQRT(2.d0*mypi)/psig*EXP(-((i-grange-1)/psig)**2/2.d0)
     ENDDO
     psf(1:2*grange+1) = psf(1:2*grange+1) / SUM(psf(1:2*grange+1))
     
     DO i=grange+1,nl-grange
        nspec(i) = SUM( psf(1:2*grange+1)*tspec(i-grange:i+grange) )
     ENDDO
     
     !interpolate back to the main array
     spec = linterp(EXP(lnlam(1:nl)),nspec(1:nl),lambda)
  
  ENDIF


END SUBROUTINE VELBROAD
