SUBROUTINE STR2ARR(switch,str,arr)

  !simple routine to dump the information in the parameteter 
  !structure (str) into an array (arr), or visa versa, depending
  !on the value of the switch

  USE sfvars
  IMPLICIT NONE

  TYPE(PARAMS), INTENT(inout) :: str
  REAL(DP), DIMENSION(npar), INTENT(inout) :: arr
  INTEGER, INTENT(in) :: switch
  INTEGER :: i

  !---------------------------------------------------------------!
  !---------------------------------------------------------------!

  IF (switch.EQ.1) THEN 

     !str->arr

     arr(1) = str%velz
     arr(2) = str%sigma
     arr(3) = str%logage
     arr(4) = str%feh
     arr(5) = str%ah
     arr(6) = str%nhe
     arr(7) = str%ch
     arr(8) = str%nh
     IF (force_nah.EQ.1) THEN
        arr(9) = str%mgh
     ELSE
        arr(9) = str%nah
     ENDIF
     arr(10) = str%mgh
     arr(11) = str%sih
     arr(12) = str%kh
     arr(13) = str%cah
     arr(14) = str%tih
     arr(15) = str%vh
     arr(16) = str%crh
     arr(17) = str%mnh
     arr(18) = str%coh
     arr(19) = str%nih
     arr(20) = str%cuh
     arr(21) = str%rbh
     arr(22) = str%srh
     arr(23) = str%yh
     arr(24) = str%zrh
     arr(25) = str%bah
     arr(26) = str%euh
     arr(27) = str%teff
     arr(28) = str%imf1
     arr(29) = str%imf2
     arr(30) = str%logfy
     arr(31) = str%sigma2
     arr(32) = str%velz2
     arr(33) = str%logm7g
     arr(34) = str%hotteff
     arr(35) = str%loghot
     arr(36) = str%fy_logage
     DO i=1,neml
        arr(npar1+i) = str%logemnorm(i)
     ENDDO
     !DO i=1,ncoeff
     !   arr(npar1+neml+i) = str%logcoeff(i)
     !ENDDO

  ELSE

     !arr->str

     str%velz   = arr(1)
     str%sigma  = arr(2)
     str%logage = arr(3)
     str%feh    = arr(4)
     str%ah     = arr(5)
     str%nhe    = arr(6)
     str%ch     = arr(7)
     str%nh     = arr(8)
     IF (force_nah.EQ.1) THEN
        str%nah = arr(10)
     ELSE
        str%nah = arr(9)
     ENDIF
     str%mgh    = arr(10)
     str%sih    = arr(11)
     str%kh     = arr(12)
     str%cah    = arr(13)
     str%tih    = arr(14)
     str%vh     = arr(15)
     str%crh    = arr(16)
     str%mnh    = arr(17)
     str%coh    = arr(18)
     str%nih    = arr(19)
     str%cuh    = arr(20)
     str%rbh    = arr(21)
     str%srh    = arr(22)
     str%yh     = arr(23)
     str%zrh    = arr(24)
     str%bah    = arr(25)
     str%euh    = arr(26)
     str%teff   = arr(27)
     str%imf1   = arr(28)
     str%imf2   = arr(29)
     str%logfy  = arr(30)
     str%sigma2 = arr(31)
     str%velz2  = arr(32)
     str%logm7g = arr(33)
     str%hotteff = arr(34)
     str%loghot = arr(35)
     str%fy_logage = arr(36)
     DO i=1,neml
        str%logemnorm(i) = arr(npar1+i)
     ENDDO
     !DO i=1,ncoeff
     !    str%logcoeff(i) = arr(npar1+neml+i)
     !ENDDO

  ENDIF


END SUBROUTINE STR2ARR
