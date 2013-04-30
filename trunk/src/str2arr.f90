SUBROUTINE STR2ARR(switch,str,arr)

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

     arr(1)=str%age
     arr(2)=str%feh
     arr(3)=str%afe
     arr(4)=str%nhe
     arr(5)=str%cfe
     arr(6)=str%nfe
     IF (force_nafe.EQ.1) THEN
        arr(7)=str%mgfe
     ELSE
        arr(7)=str%nafe
     ENDIF
     arr(8)=str%mgfe
     arr(9)=str%sife
     arr(10)=str%kfe
     arr(11)=str%cafe
     arr(12)=str%tife
     arr(13)=str%vfe
     arr(14)=str%crfe
     arr(15)=str%mnfe
     arr(16)=str%cofe
     arr(17)=str%nife
     arr(18)=str%rbfe
     arr(19)=str%srfe
     arr(20)=str%yfe
     arr(21)=str%zrfe
     arr(22)=str%bafe
     arr(23)=str%eufe
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
  !   arr(35)=str%cufe
     DO i=1,neml
        arr(34+i) = str%logemnorm(i)
     ENDDO

  ELSE

     !arr->str

     str%age=arr(1)
     str%feh=arr(2)
     str%afe=arr(3)
     str%nhe=arr(4)
     str%cfe=arr(5)
     str%nfe=arr(6)
     IF (force_nafe.EQ.1) THEN
        str%nafe=arr(8)
     ELSE
        str%nafe=arr(7)
     ENDIF
     str%mgfe=arr(8)
     str%sife=arr(9)
     str%kfe=arr(10)
     str%cafe=arr(11)
     str%tife=arr(12)
     str%vfe=arr(13)
     str%crfe=arr(14)
     str%mnfe=arr(15)
     str%cofe=arr(16)
     str%nife=arr(17)
     str%rbfe=arr(18)
     str%srfe=arr(19)
     str%yfe=arr(20)
     str%zrfe=arr(21)
     str%bafe=arr(22)
     str%eufe=arr(23)
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
 !    str%cufe=arr(35)
     DO i=1,neml
        str%logemnorm(i) = arr(34+i)
     ENDDO

  ENDIF


END SUBROUTINE STR2ARR
