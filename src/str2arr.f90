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
     arr(4) = str%feh
     arr(5) = str%ah
     arr(6) = str%ch
     arr(7) = str%nh
     arr(8) = str%nah
     arr(9) = str%mgh
     arr(10) = str%sih
     arr(11) = str%kh
     arr(12) = str%cah
     arr(13) = str%tih
     !end of the simple model parameters

     arr(14) = str%vh
     arr(15) = str%crh
     arr(16) = str%mnh
     arr(17) = str%coh
     arr(18) = str%nih
     arr(19) = str%cuh
     arr(20) = str%srh
     arr(21) = str%bah
     arr(22) = str%euh

     arr(23) = str%teff
     arr(24) = str%imf1
     arr(25) = str%imf2
     arr(26) = str%logfy
     arr(27) = str%sigma2
     arr(28) = str%velz2
     arr(29) = str%logm7g
     arr(30) = str%hotteff
     arr(31) = str%loghot
     arr(32) = str%fy_logage
     arr(33) = str%logtrans
     arr(34) = str%logemline_h
     arr(35) = str%logemline_oiii
     arr(36) = str%logemline_sii
     arr(37) = str%logemline_ni
     arr(38) = str%logemline_nii
 
  ELSE

     !arr->str

     str%velz   = arr(1)
     str%sigma  = arr(2)
     str%logage = arr(3)
     str%feh    = arr(4)
     str%ah     = arr(5)
     str%ch     = arr(6)
     str%nh     = arr(7)
     str%nah    = arr(8)
     str%mgh    = arr(9)
     str%sih    = arr(10)
     str%kh     = arr(11)
     str%cah    = arr(12)
     str%tih    = arr(13)
     !end of the simple model parameters

     str%vh     = arr(14)
     str%crh    = arr(15)
     str%mnh    = arr(16)
     str%coh    = arr(17)
     str%nih    = arr(18)
     str%cuh    = arr(19)
     str%srh    = arr(20)
     str%bah    = arr(21)
     str%euh    = arr(22)

     str%teff      = arr(23)
     str%imf1      = arr(24)
     str%imf2      = arr(25)
     str%logfy     = arr(26)
     str%sigma2    = arr(27)
     str%velz2     = arr(28)
     str%logm7g    = arr(29)
     str%hotteff   = arr(30)
     str%loghot    = arr(31)
     str%fy_logage = arr(32)
     str%logtrans  = arr(33)
     str%logemline_h    = arr(34)
     str%logemline_oiii = arr(35)
     str%logemline_sii  = arr(36)
     str%logemline_ni   = arr(37)
     str%logemline_nii  = arr(38)

  ENDIF


END SUBROUTINE STR2ARR
