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

     arr(1)=str%logage
     arr(2)=str%feh
     arr(3)=str%ah
     arr(4)=str%nhe
     arr(5)=str%ch
     arr(6)=str%nh
     IF (force_nah.EQ.1) THEN
        arr(7)=str%mgh
     ELSE
        arr(7)=str%nah
     ENDIF
     arr(8)=str%mgh
     arr(9)=str%sih
     arr(10)=str%kh
     arr(11)=str%cah
     arr(12)=str%tih
     arr(13)=str%vh
     arr(14)=str%crh
     arr(15)=str%mnh
     arr(16)=str%coh
     arr(17)=str%nih
     arr(18)=str%rbh
     arr(19)=str%srh
     arr(20)=str%yh
     arr(21)=str%zrh
     arr(22)=str%bah
     arr(23)=str%euh
     arr(24)=str%teff
     arr(25)=str%imf1
     arr(26)=str%imf2
     arr(27)=str%logfy
     arr(28)=str%sigma
     arr(29)=str%sigma2
     arr(30)=str%velz
     arr(31)=str%velz2
     arr(32)=str%logm7g
     arr(33)=str%hotteff
     arr(34)=str%loghot
  !   arr(35)=str%cuh
     DO i=1,neml
        arr(34+i) = str%logemnorm(i)
     ENDDO

  ELSE

     !arr->str

     str%logage=arr(1)
     str%feh=arr(2)
     str%ah=arr(3)
     str%nhe=arr(4)
     str%ch=arr(5)
     str%nh=arr(6)
     IF (force_nah.EQ.1) THEN
        str%nah=arr(8)
     ELSE
        str%nah=arr(7)
     ENDIF
     str%mgh=arr(8)
     str%sih=arr(9)
     str%kh=arr(10)
     str%cah=arr(11)
     str%tih=arr(12)
     str%vh=arr(13)
     str%crh=arr(14)
     str%mnh=arr(15)
     str%coh=arr(16)
     str%nih=arr(17)
     str%rbh=arr(18)
     str%srh=arr(19)
     str%yh=arr(20)
     str%zrh=arr(21)
     str%bah=arr(22)
     str%euh=arr(23)
     str%teff=arr(24)
     str%imf1=arr(25)
     str%imf2=arr(26)
     str%logfy=arr(27)
     str%sigma=arr(28)
     str%sigma2=arr(29)
     str%velz=arr(30)
     str%velz2=arr(31)
     str%logm7g=arr(32)
     str%hotteff=arr(33)
     str%loghot=arr(34)
 !    str%cuh=arr(35)
     DO i=1,neml
        str%logemnorm(i) = arr(34+i)
     ENDDO

  ENDIF


END SUBROUTINE STR2ARR
