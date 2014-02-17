SUBROUTINE READ_DATA(file)

  !routine to read in the data that will be used in the fit
  !returns a structure for the data and an integer specifying
  !the length of the data array

  USE sfvars
  IMPLICIT NONE

  CHARACTER(50), INTENT(in)  :: file
  INTEGER :: stat,i
  
  !---------------------------------------------------------------!
  !---------------------------------------------------------------!

  OPEN(10,FILE=TRIM(SPECFIT_HOME)//'/indata/'//TRIM(file)//'.dat',&
       STATUS='OLD',iostat=stat,ACTION='READ')
  IF (stat.NE.0) THEN
     WRITE(*,*) 'READ_DATA ERROR: file not found'
     STOP
  ENDIF

  DO i=1,ndat
     READ(10,*,IOSTAT=stat) data(i)%lam,data(i)%flx,data(i)%err,data(i)%wgt
     IF (data(i)%err.LE.tiny_number) data(i)%err = huge_number
     IF (stat.NE.0) GOTO 20
  ENDDO
20 CONTINUE
  CLOSE(10)

  IF (i.GT.ndat) THEN
     WRITE(*,*) 'READ_DATA ERROR: data file is too long, returning'
     STOP
  ENDIF

  datmax = i-1

END SUBROUTINE READ_DATA
