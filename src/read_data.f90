SUBROUTINE READ_DATA(file,sigma,velz)

  !routine to read in the data that will be used in the fit
  !returns a structure for the data and an integer specifying
  !the length of the data array

  USE alf_vars
  IMPLICIT NONE

  CHARACTER(50), INTENT(in)  :: file
  INTEGER :: stat,i
  CHARACTER(1) :: char
  REAL(DP) :: ll1,ll2
  REAL(DP), OPTIONAL :: sigma, velz
  
  !---------------------------------------------------------------!
  !---------------------------------------------------------------!

  CALL GETENV('ALF_HOME',ALF_HOME)
  IF (TRIM(ALF_HOME).EQ.'') THEN
     WRITE(*,*) 'ALF ERROR: ALF_HOME environment variable not set!'
     STOP
  ENDIF

  IF (fit_indices.EQ.1) THEN

     OPEN(10,FILE=TRIM(ALF_HOME)//'/indata/'//TRIM(file)//'.indx',&
          STATUS='OLD',iostat=stat,ACTION='READ')
     IF (stat.NE.0) THEN
        WRITE(*,*) 'READ_DATA ERROR: file not found'
        STOP
     ENDIF
     
     DO i=1,nindx
        READ(10,'(I1)') indx2fit(i)
     ENDDO
     READ(10,*) velz,sigma

  ELSE

     OPEN(10,FILE=TRIM(ALF_HOME)//'/indata/'//TRIM(file)//'.dat',&
          STATUS='OLD',iostat=stat,ACTION='READ')
     IF (stat.NE.0) THEN
        WRITE(*,*) 'READ_DATA ERROR: file not found'
        STOP
     ENDIF

     !-----Read in the wavelength boundaries, which are in the header-----!

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

  ENDIF


  !--------now read in the input spectrum, errors, and weights----------!

  DO i=1,ndat

     READ(10,*,IOSTAT=stat) data(i)%lam,data(i)%flx,data(i)%err,&
          data(i)%wgt,data(i)%ires
     IF (stat.NE.0) GOTO 20
     
     IF (data(i)%lam.LT.1E3.OR.data(i)%lam.GT.5E4) THEN
        WRITE(*,*) 'READ_DATA ERROR: Input lambda  <1E3 or >5E4: ',&
             data(i)%lam,data(i)%wgt
        STOP
     ENDIF
     IF (data(i)%wgt.LT.0.0.OR.data(i)%wgt.GT.1.0) THEN
        WRITE(*,*) 'READ_DATA ERROR: Input weight <0 or >1: ',&
             data(i)%lam,data(i)%wgt
        STOP
     ENDIF
     IF (data(i)%ires.LT.0.0.OR.data(i)%ires.GT.1E4) THEN
        WRITE(*,*) 'READ_DATA ERROR: Input instr. res <0 or >1E4: ',&
             data(i)%lam,data(i)%ires
        STOP
     ENDIF

     !NB: err<0.0 set to large number
     IF (data(i)%err.LE.tiny_number) data(i)%err = huge_number
 
  ENDDO

20 CONTINUE
  CLOSE(10)

  IF (i.GT.ndat) THEN
     WRITE(*,*) 'READ_DATA ERROR: data file length exceeds ndat, returning'
     STOP
  ENDIF

  datmax = i-1

END SUBROUTINE READ_DATA
