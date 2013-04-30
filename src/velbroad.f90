SUBROUTINE VELBROAD(lambda,spec,sigma)

  !routine to compute velocity broadening of an input spectrum
  !the PSF kernel has a width of m*sigma, where m=4

  USE sfvars; USE nr, ONLY : locate; USE sfutils, ONLY : linterp
  IMPLICIT NONE
  
  REAL(DP), INTENT(in), DIMENSION(:) :: lambda
  REAL(DP), INTENT(inout), DIMENSION(:) :: spec
  REAL(DP), INTENT(in) :: sigma
  REAL(DP), DIMENSION(100000) :: tspec,nspec,vel,func,gauss,psf
  REAL(DP) :: cg,xmax,xmin,fwhm,psig
  INTEGER :: i,j,il,ih,m=4,grange
  !convolve at fixed sigma_velocity
  !if set to 0, then the convolution is at fixed sigma_wavelength
  INTEGER, PARAMETER :: velocity=1

  !---------------------------------------------------------------!
  !---------------------------------------------------------------!
 

  tspec(1:nl) = linterp(LOG(lambda(1:nl)),spec(1:nl),lnlam(1:nl))
     
  fwhm   = sigma*2.35482/clight*1E5/dlstep
  psig   = fwhm/2.d0/SQRT(-2.d0*LOG(0.5d0))  ! equivalent sigma for kernel
  grange = FLOOR(m*psig)	             ! range for kernel (-range:range)
  DO i=1,2*grange+1
     psf(i) = 1.d0/SQRT(2.d0*mypi)/psig*EXP(-((i-grange-1)/psig)**2/2.d0)
  ENDDO
  psf(1:2*grange+1) = psf(1:2*grange+1) / SUM(psf(1:2*grange+1))
  
  DO i=grange+1,nl-grange
     nspec(i) = SUM( psf(1:2*grange+1)*tspec(i-grange:i+grange) )
  ENDDO
  
  !interpolate back to the main array
  spec(1:nl) = linterp(EXP(lnlam(1:nl)),nspec(1:nl),lambda(1:nl))
  

END SUBROUTINE VELBROAD
