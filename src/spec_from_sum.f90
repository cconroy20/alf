PROGRAM SPEC_FROM_SUM

  !takes a *sum file as input and returns the corresponding 
  !model spectrum associated with min(chi^2)

  USE alf_vars; USE alf_utils
  USE nr, ONLY : gasdev,locate,powell,ran1
  USE ran_state, ONLY : ran_seed,ran_init

  IMPLICIT NONE

  INTEGER  :: i,stat
  REAL(DP), DIMENSION(nl) :: mspec,lam,zmspec
  REAL(DP) :: d1,oneplusz
  CHARACTER(50)  :: infile=''
  TYPE(PARAMS) :: pos
  CHARACTER(1) :: char
  REAL(DP), DIMENSION(npar) :: posarr=0.0
  REAL(DP), DIMENSION(6)    :: mlx2=0.0

  !-----------------------------------------------------------!
  !-----------------------------------------------------------!

  imf_type = 1
  l1(1) = 0.0
  nlint = 1
  l2(nlint) = 1E5

  IF (IARGC().LT.1) THEN
     WRITE(*,*) 'ERROR: no input file!'
     STOP
  ELSE
     CALL GETARG(1,infile)
  ENDIF

  CALL GETENV('ALF_HOME',ALF_HOME)

  !initialize the random number generator
  CALL INIT_RANDOM_SEED()

  !read in parameters from the *sum file
  OPEN(11,FILE=TRIM(ALF_HOME)//'results/'//TRIM(infile)//'.sum',&
       STATUS='OLD',iostat=stat,ACTION='READ')
  IF (stat.NE.0) THEN
     WRITE(*,*) 'ERROR, file not found: ', infile
     STOP
  ENDIF
  
  char = '#'
  DO WHILE (char.EQ.'#')
     READ(11,*,IOSTAT=stat) char
  ENDDO
  BACKSPACE(11)
  READ(11,*) !burn the row containing the mean parameters
  READ(11,*) d1,posarr,mlx2
  CLOSE(11)

  !copy the input parameter array into the structure
  CALL STR2ARR(2,pos,posarr) !arr->str

  !setup the models
  CALL SETUP()
  lam = sspgrid%lam

  !we turn off the emission lines, since they are often highly
  !unconstrained if they are not included in the wavelength range
  pos%logemline_h    = -8.0
  pos%logemline_oii  = -8.0
  pos%logemline_oiii = -8.0
  pos%logemline_nii  = -8.0
  pos%logemline_sii  = -8.0
  pos%logemline_ni   = -8.0


  !------------------------------------------------------------!
  !here is the place to make changes to the best-fit model,
  !if so desired

  pos%loghot = -8.0

  !------------------------------------------------------------!


  !get the model spectrum
  CALL GETMODEL(pos,mspec)
     
  !redshift the spectrum
  oneplusz = (1+pos%velz/clight*1E5)
  zmspec   = 0.0
  zmspec   = LINTERP(lam*oneplusz,mspec,lam)

  !print to file
  OPEN(12,FILE=TRIM(ALF_HOME)//'models/'//TRIM(infile)//'.bestspec2',&
       STATUS='REPLACE')
  DO i=1,nl
     WRITE(12,'(F10.3,2ES12.4)') lam(i),zmspec(i)
  ENDDO
  CLOSE(12)


END PROGRAM SPEC_FROM_SUM
