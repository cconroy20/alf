SUBROUTINE STR2ARR(switch,str,arr)

  !routine to dump the information in the parameteter 
  !structure (str) into an array (arr), or visa versa, depending
  !on the value of the switch

  !this is actually the location where the parameters to be fit
  !in the full model are specified.  If a parameter is left out
  !of this list (e.g., [Rb/H]) then it does not affect Chi^2

  USE alf_vars
  IMPLICIT NONE

  TYPE(PARAMS), INTENT(inout) :: str
  REAL(DP), DIMENSION(npar), INTENT(inout) :: arr
  INTEGER, INTENT(in) :: switch

  !---------------------------------------------------------------!
  !---------------------------------------------------------------!

  IF (switch.EQ.1) THEN 

     !str->arr

     arr(1) = str%velz
     arr(2) = str%sigma
     arr(3) = str%logage
     arr(4) = str%zh
     !end of the super-simple and Powell-mode parameters

     arr(5) = str%feh
     arr(6) = str%ah
     arr(7) = str%ch
     arr(8) = str%nh
     arr(9) = str%nah
     arr(10) = str%mgh
     arr(11) = str%sih
     arr(12) = str%kh
     arr(13) = str%cah
     arr(14) = str%tih
     !end of the simple model parameters

     arr(15) = str%vh
     arr(16) = str%crh
     arr(17) = str%mnh
     arr(18) = str%coh
     arr(19) = str%nih
     arr(20) = str%cuh
     arr(21) = str%srh
     arr(22) = str%bah
     arr(23) = str%euh

     arr(24) = str%teff
     arr(25) = str%imf1
     arr(26) = str%imf2
     arr(27) = str%logfy
     arr(28) = str%sigma2
     arr(29) = str%velz2
     arr(30) = str%logm7g
     arr(31) = str%hotteff
     arr(32) = str%loghot
     arr(33) = str%fy_logage
     arr(34) = str%logtrans
     arr(35) = str%logemline_h
     arr(36) = str%logemline_oiii
     arr(37) = str%logemline_sii
     arr(38) = str%logemline_ni
     arr(39) = str%logemline_nii
     arr(40) = str%jitter
     arr(41) = str%imf3
     arr(42) = str%logsky
     arr(43) = str%imf4
     arr(44) = str%h3
     arr(45) = str%h4
     arr(46) = str%yzh

  ELSE

     !arr->str

     str%velz   = arr(1)
     str%sigma  = arr(2)
     str%logage = arr(3)
     str%zh     = arr(4)
     !end of the super-simple and Powell-mode parameters

     str%feh    = arr(5)
     str%ah     = arr(6)
     str%ch     = arr(7)
     str%nh     = arr(8)
     str%nah    = arr(9)
     str%mgh    = arr(10)
     str%sih    = arr(11)
     str%kh     = arr(12)
     str%cah    = arr(13)
     str%tih    = arr(14)
     !end of the simple model parameters

     str%vh     = arr(15)
     str%crh    = arr(16)
     str%mnh    = arr(17)
     str%coh    = arr(18)
     str%nih    = arr(19)
     str%cuh    = arr(20)
     str%srh    = arr(21)
     str%bah    = arr(22)
     str%euh    = arr(23)

     str%teff      = arr(24)
     str%imf1      = arr(25)
     str%imf2      = arr(26)
     str%logfy     = arr(27)
     str%sigma2    = arr(28)
     str%velz2     = arr(29)
     str%logm7g    = arr(30)
     str%hotteff   = arr(31)
     str%loghot    = arr(32)
     str%fy_logage = arr(33)
     str%logtrans  = arr(34)
     str%logemline_h    = arr(35)
     str%logemline_oiii = arr(36)
     str%logemline_sii  = arr(37)
     str%logemline_ni   = arr(38)
     str%logemline_nii  = arr(39)
     str%jitter = arr(40)
     str%imf3   = arr(41)
     str%logsky = arr(42)
     str%imf4   = arr(43)
     str%h3     = arr(44)
     str%h4     = arr(45)
     str%yzh    = arr(46)
     
  ENDIF


END SUBROUTINE STR2ARR
