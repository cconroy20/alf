SUBROUTINE SETUP_PARAMS(pos,prlo,prhi,velz)

  !define the first position (pos), and the lower and upper bounds 
  !on the priors (prlo, prhi).  The priors are defined in such a way
  !that if the user defines a prior limit that is **different from
  !the default parameter set**, then that value overrides the defaults below

  USE sfvars; USE sfutils, ONLY : myran
  IMPLICIT NONE

  TYPE(PARAMS), INTENT(inout) :: pos,prlo,prhi
  TYPE(PARAMS) :: test
  REAL(DP), OPTIONAL :: velz
  INTEGER :: i
  
  !---------------------------------------------------------------!
  !---------------------------------------------------------------!

  !setup the first position
  pos%logage    = myran()*0.2+0.8
  pos%feh       = myran()*1.0-0.5
  pos%ah        = myran()*1.0-0.5
  pos%nhe       = myran()*1.0-0.5
  pos%ch        = myran()*1.0-0.5
  pos%nh        = myran()*1.0-0.5
  pos%nah       = myran()*1.0-0.5
  pos%mgh       = myran()*1.0-0.5
  pos%sih       = myran()*1.0-0.5
  pos%kh        = myran()*1.0-0.5
  pos%cah       = myran()*1.0-0.5
  pos%tih       = myran()*1.0-0.5
  pos%vh        = myran()*1.0-0.5
  pos%crh       = myran()*1.0-0.5
  pos%mnh       = myran()*1.0-0.5
  pos%coh       = myran()*1.0-0.5
  pos%nih       = myran()*1.0-0.5
  pos%cuh       = myran()*1.0-0.5
  pos%rbh       = myran()*1.0-0.5
  pos%srh       = myran()*1.0-0.5
  pos%yh        = myran()*1.0-0.5
  pos%zrh       = myran()*1.0-0.5
  pos%bah       = myran()*1.0-0.5
  pos%euh       = myran()*1.0-0.5
  pos%teff      = myran()*80-40
  pos%imf1      = myran()*0.6-0.3 + 1.3
  pos%imf2      = myran()*0.6-0.3 + 2.3
  pos%logfy     = myran()*2-3
  pos%fy_logage = myran()*0.6-0.2
  pos%logm7g    = myran()*2-4
  pos%hotteff   = myran()*5+15
  pos%loghot    = myran()*2-3
  pos%chi2      = huge_number
  pos%sigma     = myran()*100+100
  pos%sigma2    = myran()*100+100
  pos%velz2     = myran()*10-5
  DO i=1,neml
     pos%logemnorm(i) = myran()*2-4
  ENDDO
  IF (PRESENT(velz)) THEN
     pos%velz  = velz + (myran()*10-5)
  ELSE
     pos%velz  = myran()*10-5
  ENDIF


  !priors (low)
  IF (prlo%logage.EQ.test%logage) prlo%logage = LOG10(0.5)
  IF (prlo%feh.EQ.test%feh) prlo%feh          = -1.0
  IF (prlo%ah.EQ.test%ah) prlo%ah             = -1.0
  IF (prlo%nhe.EQ.test%nhe) prlo%nhe          = -1.0
  IF (prlo%ch.EQ.test%ch) prlo%ch             = -1.0
  IF (prlo%nh.EQ.test%nh) prlo%nh             = -1.0
  IF (prlo%nah.EQ.test%nah) prlo%nah          = -1.0
  IF (prlo%mgh.EQ.test%mgh) prlo%mgh          = -1.0
  IF (prlo%sih.EQ.test%sih) prlo%sih          = -1.0
  IF (prlo%kh.EQ.test%kh) prlo%kh             = -1.0
  IF (prlo%cah.EQ.test%cah) prlo%cah          = -1.0
  IF (prlo%tih.EQ.test%tih) prlo%tih          = -1.0
  IF (prlo%vh.EQ.test%vh) prlo%vh             = -1.0
  IF (prlo%crh.EQ.test%crh) prlo%crh          = -1.0
  IF (prlo%mnh.EQ.test%mnh) prlo%mnh          = -1.0
  IF (prlo%coh.EQ.test%coh) prlo%coh          = -1.0
  IF (prlo%nih.EQ.test%nih) prlo%nih          = -1.0
  IF (prlo%cuh.EQ.test%cuh) prlo%cuh          = -1.0
  IF (prlo%rbh.EQ.test%rbh) prlo%rbh          = -1.0
  IF (prlo%srh.EQ.test%srh) prlo%srh          = -1.0
  IF (prlo%yh.EQ.test%yh) prlo%yh             = -1.0
  IF (prlo%zrh.EQ.test%zrh) prlo%zrh          = -1.0
  IF (prlo%feh.EQ.test%feh) prlo%bah          = -1.0
  IF (prlo%euh.EQ.test%euh) prlo%euh          = -1.0
  IF (prlo%teff.EQ.test%teff) prlo%teff       = -100.0
  IF (prlo%imf1.EQ.test%imf1) prlo%imf1       = 0.5
  IF (prlo%imf2.EQ.test%imf2) prlo%imf2       = 0.5
  IF (prlo%logfy.EQ.test%logfy) prlo%logfy    = -5.0
  IF (prlo%fy_logage.EQ.test%fy_logage) prlo%fy_logage = LOG10(0.5)
  IF (prlo%logm7g.EQ.test%logm7g) prlo%logm7g   = -5.0
  IF (prlo%hotteff.EQ.test%hotteff) prlo%hotteff= 8.0
  IF (prlo%loghot.EQ.test%loghot) prlo%loghot   = -5.0
  IF (prlo%sigma.EQ.test%sigma) prlo%sigma      = 20.0
  IF (prlo%sigma2.EQ.test%sigma2) prlo%sigma2   = 20.0
  IF (prlo%velz.EQ.test%velz) prlo%velz         = -1E3
  IF (prlo%velz2.EQ.test%velz2) prlo%velz2      = -1E3
  DO i=1,neml 
     IF (prlo%logemnorm(i).EQ.test%logemnorm(i)) prlo%logemnorm(i) = -8.0
  ENDDO

  !priors (high)
  IF (prhi%logage.EQ.test%logage) prhi%logage = LOG10(13.7)
  IF (prhi%feh.EQ.test%feh) prhi%feh          = 1.0
  IF (prhi%ah.EQ.test%ah) prhi%ah             = 1.0
  IF (prhi%nhe.EQ.test%nhe) prhi%nhe          = 1.0
  IF (prhi%ch.EQ.test%ch) prhi%ch             = 1.0
  IF (prhi%nh.EQ.test%nh) prhi%nh             = 1.0
  IF (prhi%nah.EQ.test%nah) prhi%nah          = 1.0
  IF (prhi%mgh.EQ.test%mgh) prhi%mgh          = 1.0
  IF (prhi%sih.EQ.test%sih) prhi%sih          = 1.0
  IF (prhi%kh.EQ.test%kh) prhi%kh             = 1.0
  IF (prhi%cah.EQ.test%cah) prhi%cah          = 1.0
  IF (prhi%tih.EQ.test%tih) prhi%tih          = 1.0
  IF (prhi%vh.EQ.test%vh) prhi%vh             = 1.0
  IF (prhi%crh.EQ.test%crh) prhi%crh          = 1.0
  IF (prhi%mnh.EQ.test%mnh) prhi%mnh          = 1.0
  IF (prhi%coh.EQ.test%coh) prhi%coh          = 1.0
  IF (prhi%nih.EQ.test%nih) prhi%nih          = 1.0
  IF (prhi%cuh.EQ.test%cuh) prhi%cuh          = 1.0
  IF (prhi%rbh.EQ.test%rbh) prhi%rbh          = 1.0
  IF (prhi%srh.EQ.test%srh) prhi%srh          = 1.0
  IF (prhi%yh.EQ.test%yh) prhi%yh             = 1.0
  IF (prhi%zrh.EQ.test%zrh) prhi%zrh          = 1.0
  IF (prhi%feh.EQ.test%feh) prhi%bah          = 1.0
  IF (prhi%euh.EQ.test%euh) prhi%euh          = 1.0
  IF (prhi%teff.EQ.test%teff) prhi%teff       = 100.0
  IF (prhi%imf1.EQ.test%imf1) prhi%imf1       = 3.5
  IF (prhi%imf2.EQ.test%imf2) prhi%imf2       = 3.5
  IF (prhi%logfy.EQ.test%logfy) prhi%logfy    = -0.3
  IF (prhi%fy_logage.EQ.test%fy_logage) prhi%fy_logage = LOG10(5.0)
  IF (prhi%logm7g.EQ.test%logm7g) prhi%logm7g   = -0.5
  IF (prhi%hotteff.EQ.test%hotteff) prhi%hotteff= 30.0
  IF (prhi%loghot.EQ.test%loghot) prhi%loghot   = -0.1
  IF (prhi%sigma.EQ.test%sigma) prhi%sigma      = 1E3
  IF (prhi%sigma2.EQ.test%sigma2) prhi%sigma2   = 1E3
  IF (prhi%velz.EQ.test%velz) prhi%velz         = 1E4
  IF (prhi%velz2.EQ.test%velz2) prhi%velz2      = 1E3
  DO i=1,neml 
     IF (prhi%logemnorm(i).EQ.test%logemnorm(i)) prhi%logemnorm(i) = 2.0
  ENDDO

END SUBROUTINE SETUP_PARAMS
