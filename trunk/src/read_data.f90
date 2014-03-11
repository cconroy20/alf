SUBROUTINE READ_DATA(file)

  !routine to read in the data that will be used in the fit
  !returns a structure for the data and an integer specifying
  !the length of the data array

  USE sfvars
  IMPLICIT NONE

  CHARACTER(50), INTENT(in)  :: file
  INTEGER :: stat,i
  CHARACTER(1) :: char
  REAL(DP) :: ll1,ll2
  
  !---------------------------------------------------------------!
  !---------------------------------------------------------------!

  OPEN(10,FILE=TRIM(SPECFIT_HOME)//'/indata/'//TRIM(file)//'.dat',&
       STATUS='OLD',iostat=stat,ACTION='READ')
  IF (stat.NE.0) THEN
     WRITE(*,*) 'READ_DATA ERROR: file not found'
     STOP
  ENDIF

  !--------Read in the wavelength boundaries, which are in the header------!

  char='#'
  nlint = 0
  DO WHILE (char.EQ.'#') 
     READ(10,*) char,ll1,ll2
     IF (char.EQ.'#') THEN 
        nlint = nlint+1
        IF (nlint.GT.nlint_max) THEN 
           WRITE(*,*) 'READ_DATA ERROR: number of wavelength '//&
                'intervals exceeds nlint_max'
           STOP
        ENDIF
        IF (ll1.GE.ll2) THEN 
           WRITE(*,*) 'READ_DATA ERROR: l1>l2!  returning...'
           STOP
        ENDIF
        l1(nlint) = ll1
        l2(nlint) = ll2
     ENDIF
  ENDDO
  BACKSPACE(10)

  IF (nlint.EQ.0) THEN
     WRITE(*,*) 'no wavelength boundaries specified, using default' 
     nlint = 2
     l1(1) = 0.40
     l1(2) = 0.47
     l2(1) = 0.47
     l2(2) = 0.55
  ENDIF

  !convert from um to A.
  l1 = l1*1E4
  l2 = l2*1E4

  !--------now read in the input spectrum, errors, and weights----------!

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
