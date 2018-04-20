PROGRAM ML_FROM_MCMC

  ! takes an *mcmc file as input and returns M/L in many filters
  ! note that M/L is returned in the *observed* frame
  
  USE alf_vars; USE alf_utils
  USE nr, ONLY : locate

  IMPLICIT NONE

  INTEGER, PARAMETER :: nfil2=27
  INTEGER  :: i,j,k,stat,nwalkers,nchain,nsample
  REAL(DP), DIMENSION(nl) :: mspec=0.,lam=0.,zmspec=0.
  REAL(DP) :: d1=0.,oneplusz=0.,mass=0.,msto=0.
  CHARACTER(50)  :: infile='',line=''
  TYPE(PARAMS) :: pos
  CHARACTER(1) :: char=''
  REAL(DP), DIMENSION(npar)  :: posarr=0.0
  REAL(DP), DIMENSION(nl,nfil2) :: filters2=0.
  REAL(DP), DIMENSION(nfil2) :: m2l=0.0,mag=0.
  REAL(DP), DIMENSION(6)     :: mlx2=0.0

  !-----------------------------------------------------------!
  !-----------------------------------------------------------!

  nlint = 1
  l1(1) = 0.0
  l2(nlint) = 2.5E4

  IF (IARGC().LT.1) THEN
     WRITE(*,*) 'ERROR: no input file!'
     STOP
  ELSE
     CALL GETARG(1,infile)
  ENDIF

  CALL GETENV('ALF_HOME',ALF_HOME)

  !read in input from the *sum file
  OPEN(11,FILE=TRIM(ALF_HOME)//'results/'//TRIM(infile)//'.sum',&
       STATUS='OLD',iostat=stat,ACTION='READ')
  IF (stat.NE.0) THEN
     WRITE(*,*) 'ERROR, file not found: ', infile
     STOP
  ENDIF

  !open file for output
  OPEN(12,FILE=TRIM(ALF_HOME)//'models/'//TRIM(infile)//'.ml_all',&
       STATUS='REPLACE')

  !read in the header and set the relevant parameters
  char = '#'
  DO WHILE (char.EQ.'#')
     READ(11,'(A50)',IOSTAT=stat) line
     char = line(1:1)
     IF (index(line,'mwimf').GT.0) THEN
        read(line(index(line,'=')+2:index(line,'=')+3),*) mwimf
     ENDIF
     IF (index(line,'imf_type').GT.0) THEN
        read(line(index(line,'=')+2:index(line,'=')+3),*) imf_type
     ENDIF
     IF (index(line,'fit_type').GT.0) THEN
        read(line(index(line,'=')+2:index(line,'=')+3),*) fit_type
     ENDIF
     IF (index(line,'fit_two_ages').GT.0) THEN
        read(line(index(line,'=')+2:index(line,'=')+3),*) fit_two_ages
     ENDIF
     IF (index(line,'fit_hermite').GT.0) THEN
        read(line(index(line,'=')+2:index(line,'=')+3),*) fit_hermite
     ENDIF
     IF (index(line,'nonpimf').GT.0) THEN
        read(line(index(line,'=')+2:index(line,'=')+3),*) nonpimf_alpha
     ENDIF
     IF (index(line,'Nwalkers').GT.0) THEN
        read(line(index(line,'=')+2:index(line,'=')+8),*) nwalkers
     ENDIF
     IF (index(line,'Nchain').GT.0) THEN
        read(line(index(line,'=')+2:index(line,'=')+8),*) nchain
     ENDIF
  ENDDO
  CLOSE(11)

  !read in parameters from the *mcmc file
  OPEN(11,FILE=TRIM(ALF_HOME)//'results/'//TRIM(infile)//'.mcmc',&
       STATUS='OLD',iostat=stat,ACTION='READ')
  IF (stat.NE.0) THEN
     WRITE(*,*) 'ERROR, file not found: ', infile
     STOP
  ENDIF

  !setup the models
  CALL SETUP()
  lam = sspgrid%lam

  !read in filter transmission curves
  OPEN(15,FILE=TRIM(ALF_HOME)//'/infiles/filters2.dat',&
       STATUS='OLD',iostat=stat,ACTION='READ')
  IF (stat.NE.0) THEN
     WRITE(*,*) 'SETUP ERROR: filter curves not found'
     STOP
  ENDIF
  READ(15,*) magsun
  DO i=1,nstart-1
     READ(15,*) 
  ENDDO
  DO i=1,nl
     READ(15,*) d1,filters2(i,:)
  ENDDO
  CLOSE(15)
  
  !----------------------------------------------------------------------!
  !----------------------------------------------------------------------!

  
  DO j=1,nchain
     DO k=1,nwalkers

        READ(11,*) d1,posarr,mlx2

        !copy the input parameter array into the structure
        CALL STR2ARR(2,pos,posarr) !arr->str

        !turn off the emission lines
        pos%logemline_h    = -8.0
        pos%logemline_oii  = -8.0
        pos%logemline_oiii = -8.0
        pos%logemline_nii  = -8.0
        pos%logemline_sii  = -8.0
        pos%logemline_ni   = -8.0

        !get the model spectrum
        CALL GETMODEL(pos,mspec)
        
        !redshift the spectrum
        oneplusz = (1+pos%velz/clight*1E5)
        zmspec   = MAX(LINTERP(lam*oneplusz,mspec,lam),0.0)
        !convert to proper units
        zmspec  = zmspec*lsun/1E6*lam**2/clight/1E8/4/mypi/pc2cm**2

        !MSTO
        msto = 10**(msto_t0+msto_t1*pos%logage) * &
             ( msto_z0 + msto_z1*pos%zh + msto_z2*pos%zh**2 )

        !stellar mass
        IF (imf_type.EQ.0) THEN
           mass = getmass(imflo,msto,pos%imf1,pos%imf1,krpa_imf3)
        ELSE IF (imf_type.EQ.1) THEN
           mass = getmass(imflo,msto,pos%imf1,pos%imf2,krpa_imf3)
        ELSE IF (imf_type.EQ.2) THEN
           mass = getmass(pos%imf3,msto,pos%imf1,pos%imf1,krpa_imf3)
        ELSE IF (imf_type.EQ.3) THEN
           mass = getmass(pos%imf3,msto,pos%imf1,pos%imf2,krpa_imf3)
        ELSE IF (imf_type.EQ.4) THEN
           mass = getmass(imflo,msto,pos%imf1,pos%imf2,krpa_imf3,&
                pos%imf3,pos%imf4)
        ENDIF

        !loop over the filters
        DO i=1,nfil2
           mag(i) = tsum(lam,zmspec*filters2(:,i)/lam)
           IF (mag(i).LE.0.0) THEN
              m2l(i) = 0.0
           ELSE 
              mag(i) = -2.5*LOG10(mag(i))-48.60
              m2l(i) = mass/10**(2./5*(magsun2(i)-mag(i)))
              IF (m2l(i).GT.100.) m2l(i)=0.0
           ENDIF
        ENDDO


        WRITE(12,'(50F7.2)') m2l

     ENDDO
  ENDDO
        
  CLOSE(11)
  CLOSE(12)


END PROGRAM ML_FROM_MCMC
