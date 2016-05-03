SUBROUTINE SETUP()

  !read in and set up all the arrays

  USE alf_vars; USE nr, ONLY : locate; USE alf_utils, ONLY : linterp,velbroad
  IMPLICIT NONE
  
  REAL(DP) :: d1,d2,d3,l1um=1E4,t13=1.3,t23=2.3,m07=0.07,sig0=99.,lamlo,lamhi
  REAL(DP), DIMENSION(nimf*nimf) :: tmp
  REAL(DP), DIMENSION(nl) :: dumi,smooth=0.0,lam
  INTEGER :: stat,i,vv,j,k,t,z,ii,shift=100
  INTEGER, PARAMETER :: ntrans=22800
  REAL(DP), DIMENSION(ntrans) :: ltrans,ftrans,strans
  CHARACTER(4), DIMENSION(nzmet) :: charz
  CHARACTER(5), DIMENSION(nage) :: chart
  CHARACTER(20) :: imfstr

  !---------------------------------------------------------------!
  !---------------------------------------------------------------!

  CALL GETENV('SPECFIT_HOME',SPECFIT_HOME)
  IF (TRIM(SPECFIT_HOME).EQ.'') THEN
     WRITE(*,*) 'ALF ERROR: SPECFIT_HOME environment variable not set!'
     STOP
  ENDIF

  charz = (/'m1.5','m1.0','m0.5','p0.0','p0.3'/)
  chart = (/'t01.0','t03.0','t05.0','t07.0','t09.0','t11.0','t13.5'/)

  IF (fit_2ximf.EQ.1) THEN
     imfstr = 'varydoublex'
  ELSE
     imfstr = 'varymcut_varyx'
  ENDIF

  !if the data has not been read in, then we need to manually
  !define the lower and upper limits for the instrumental resolution
  !broadening.  Currently this is only triggered if write_a_model is 
  !being called (or if an error is made in READ_DATA).
  IF (nlint.EQ.0) THEN
     lamlo = 3.8E3
     lamhi = 2.4E4
  ELSE
     lamlo = l1(1)-500
     lamhi = l2(nlint)+500
  ENDIF

  !read in filter transmission curves
  OPEN(15,FILE=TRIM(SPECFIT_HOME)//'/infiles/filters.dat',&
       STATUS='OLD',iostat=stat,ACTION='READ')
  DO i=1,nstart-1
     READ(15,*) 
  ENDDO
  DO i=1,nl
     READ(15,*) d1,filters(i,1),filters(i,2),filters(i,3) !r,i,K filters
  ENDDO
  CLOSE(15)

  !read in the ATLAS SSPs
  DO k=1,nzmet

     DO j=1,nage_rfcn

        IF (j.EQ.1) THEN
           OPEN(20,FILE=TRIM(SPECFIT_HOME)//'/infiles/atlas_ssp_t01_Z'//&
                charz(k)//'.abund.krpa.s100',STATUS='OLD',iostat=stat,ACTION='READ')
        ELSE IF (j.EQ.2) THEN
           OPEN(20,FILE=TRIM(SPECFIT_HOME)//'/infiles/atlas_ssp_t03_Z'//&
                charz(k)//'.abund.krpa.s100',STATUS='OLD',iostat=stat,ACTION='READ')
        ELSE IF (j.EQ.3) THEN
           OPEN(20,FILE=TRIM(SPECFIT_HOME)//'/infiles/atlas_ssp_t05_Z'//&
                charz(k)//'.abund.krpa.s100',STATUS='OLD',iostat=stat,ACTION='READ')
        ELSE IF (j.EQ.4) THEN
           OPEN(20,FILE=TRIM(SPECFIT_HOME)//'/infiles/atlas_ssp_t09_Z'//&
                charz(k)//'.abund.krpa.s100',STATUS='OLD',iostat=stat,ACTION='READ')
        ELSE IF (j.EQ.5) THEN
           OPEN(20,FILE=TRIM(SPECFIT_HOME)//'/infiles/atlas_ssp_t13_Z'//&
                charz(k)//'.abund.krpa.s100',STATUS='OLD',iostat=stat,ACTION='READ')
        ENDIF

        IF (stat.NE.0) THEN
           WRITE(*,*) 'SETUP ERROR: ATLAS models not found'
           STOP
        ENDIF

        READ(20,*) !burn the header
        READ(20,*)
        DO i=1,nstart-1
           READ(20,*) 
        ENDDO
        DO i=1,nl
           READ(20,*) sspgrid%lam(i),sspgrid%solar(i,j,k),sspgrid%nap(i,j,k),&
                sspgrid%nam(i,j,k),sspgrid%cap(i,j,k),sspgrid%cam(i,j,k),sspgrid%fep(i,j,k),&
                sspgrid%fem(i,j,k),sspgrid%cp(i,j,k),sspgrid%cm(i,j,k),d1,&
                sspgrid%np(i,j,k),sspgrid%nm(i,j,k),sspgrid%ap(i,j,k),sspgrid%tip(i,j,k),&
                sspgrid%tim(i,j,k),sspgrid%mgp(i,j,k),sspgrid%mgm(i,j,k),sspgrid%sip(i,j,k),&
                sspgrid%sim(i,j,k),sspgrid%teffp(i,j,k),sspgrid%teffm(i,j,k),sspgrid%crp(i,j,k),&
                sspgrid%mnp(i,j,k),sspgrid%bap(i,j,k),sspgrid%bam(i,j,k),sspgrid%nip(i,j,k),&
                sspgrid%cop(i,j,k),sspgrid%eup(i,j,k),sspgrid%srp(i,j,k),sspgrid%kp(i,j,k),&
                sspgrid%vp(i,j,k),sspgrid%cup(i,j,k),sspgrid%nap6(i,j,k),sspgrid%nap9(i,j,k)
        ENDDO
        CLOSE(20)

     ENDDO
  ENDDO

  lam = sspgrid%lam
  sspgrid%logagegrid_rfcn = LOG10((/1.0,3.0,5.0,9.0,13.0/))


  !create fake response functions by shifting
  !the wavelengths by n pixels
  IF (fake_response.EQ.1) THEN

     WRITE(*,*) 'ERROR: this option is not currently supported'
     STOP

     !DO i=1,nage_rfcn

        !dumi = sspgrid%crp(:,i) / sspgrid%solar(:,i)
        !sspgrid%crp(:,i) = 1.0
        !sspgrid%crp(shift:nl-1,i) = dumi(1:nl-shift)
        !sspgrid%crp(:,i) = sspgrid%crp(:,i) * sspgrid%solar(:,i)
     
        !dumi = sspgrid%mnp(:,i) / sspgrid%solar(:,i)
        !sspgrid%mnp = 1.0
        !sspgrid%mnp(shift:nl-1,i) = dumi(1:nl-shift)
        !sspgrid%mnp(:,i) = sspgrid%mnp(:,i) * sspgrid%solar(:,i)
  
        !dumi = sspgrid%cop(:,i) / sspgrid%solar(:,i)
        !sspgrid%cop = 1.0
        !sspgrid%cop(shift:nl-1,i) = dumi(1:nl-shift)
        !sspgrid%cop(:,i) = sspgrid%cop(:,i) * sspgrid%solar(:,i)
 
     !ENDDO

  ENDIF


  !read in empirical spectra as a function of age and metallicity
  DO j=1,nzmet
     OPEN(21,FILE=TRIM(SPECFIT_HOME)//'/infiles/CvD_krpaIMF_Z'//charz(j)//&
          '.ssp.s100',STATUS='OLD',iostat=stat,ACTION='READ')
     IF (stat.NE.0) THEN
        WRITE(*,*) 'SETUP ERROR: empirical models not found'
        STOP
     ENDIF
     DO i=1,nstart-1
        READ(21,*) 
     ENDDO
     DO i=1,nl
        !order of the logfkrpa array is lam,age,zmet
        READ(21,*) d1,sspgrid%logfkrpa(i,1,j),sspgrid%logfkrpa(i,2,j),&
             sspgrid%logfkrpa(i,3,j),sspgrid%logfkrpa(i,4,j),sspgrid%logfkrpa(i,5,j),&
             sspgrid%logfkrpa(i,6,j),sspgrid%logfkrpa(i,7,j)
     ENDDO
     CLOSE(21)
  ENDDO

  sspgrid%logfkrpa   = LOG10(sspgrid%logfkrpa+tiny_number)
  sspgrid%logagegrid = LOG10((/1.0,3.0,5.0,7.0,9.0,11.0,13.5/))
  sspgrid%logzgrid   = (/-1.5,-1.0,-0.5,0.0,0.3/)

  DO z=1,nzmet
     DO t=1,nage
        !read in two parameter IMF models
        OPEN(22,FILE=TRIM(SPECFIT_HOME)//'/infiles/CvD.v3_'//chart(t)//'_Z'//&
             charz(z)//'.ssp.'//'imf_'//TRIM(imfstr)//'.s100',STATUS='OLD',&
             iostat=stat,ACTION='READ')
        IF (stat.NE.0) THEN
           WRITE(*,*) 'SETUP ERROR: IMF models not found'
           STOP
        ENDIF
        DO i=1,nstart-1
           READ(22,*) 
        ENDDO
        DO i=1,nl
           READ(22,*) d1,tmp
           ii=1
           DO j=1,nimf
              DO k=1,nimf
                 sspgrid%imf(i,j,k,t,z) = tmp(ii)
                 ii=ii+1
              ENDDO
           ENDDO
        ENDDO
        CLOSE(22)
     ENDDO
  ENDDO

  !values of IMF parameters at the 16 grid points
  DO i=1,nimf
     sspgrid%imfx1(i) = 0.d5+REAL(i-1)/5.d0 
  ENDDO
  IF (fit_2ximf.EQ.1) THEN
     DO i=1,nimf
        sspgrid%imfx2(i) = 0.d5+REAL(i-1)/5.d0 
     ENDDO
     imfr1 = locate(sspgrid%imfx1,t13+1E-3)
     imfr2 = locate(sspgrid%imfx2,t23+1E-3)
  ELSE
     sspgrid%imfx2 = (/0.07,0.10,0.15,0.2,0.25,0.3,0.35,0.4,&
          0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8/)
     imfr1 = locate(sspgrid%imfx1,t23+1E-3)
     imfr2 = locate(sspgrid%imfx2,m07+1E-3)
  ENDIF

  !read in M7III star, normalized to a 13 Gyr SSP at 1um
  OPEN(23,FILE=TRIM(SPECFIT_HOME)//'/infiles/M7III.spec.s100',&
       STATUS='OLD',iostat=stat,ACTION='READ')
  DO i=1,nstart-1
     READ(23,*) 
  ENDDO
  DO i=1,nl
     READ(23,*) d1,sspgrid%m7g(i)
  ENDDO
  CLOSE(23)

  !read in hot stars 
  OPEN(24,FILE=TRIM(SPECFIT_HOME)//'/infiles/ap00t08000g4.00at12.spec.s100',&
       STATUS='OLD',iostat=stat,ACTION='READ')
  DO i=1,nstart-1
     READ(24,*) 
  ENDDO
  DO i=1,nl
     READ(24,*) d1,sspgrid%hotspec(i,1)
  ENDDO
  CLOSE(24)
  OPEN(25,FILE=TRIM(SPECFIT_HOME)//'/infiles/ap00t10000g4.00at12.spec.s100',&
       STATUS='OLD',iostat=stat,ACTION='READ')
  DO i=1,nstart-1
     READ(25,*) 
  ENDDO
  DO i=1,nl
     READ(25,*) d1,sspgrid%hotspec(i,2)
  ENDDO
  CLOSE(25)
  OPEN(26,FILE=TRIM(SPECFIT_HOME)//'/infiles/ap00t20000g4.00at12.spec.s100',&
       STATUS='OLD',iostat=stat,ACTION='READ')
  DO i=1,nstart-1
     READ(26,*) 
  ENDDO
  DO i=1,nl
     READ(26,*) d1,sspgrid%hotspec(i,3)
  ENDDO
  CLOSE(26)
  OPEN(27,FILE=TRIM(SPECFIT_HOME)//'/infiles/ap00t30000g4.00at12.spec.s100',&
       STATUS='OLD',iostat=stat,ACTION='READ')
  DO i=1,nstart-1
     READ(27,*) 
  ENDDO
  DO i=1,nl
     READ(27,*) d1,sspgrid%hotspec(i,4)
  ENDDO
  CLOSE(27)
 
  !normalize to a 13 Gyr SSP at 1um (same norm for the M7III param)
  !NB: this normalization was changed on 7/20/15.  Also, a major
  !bug was found in which the indices of the array were reversed.
  vv = locate(sspgrid%lam(1:nl),l1um)
  DO i=1,nhot
     sspgrid%hotspec(:,i) = sspgrid%hotspec(:,i)/sspgrid%hotspec(vv,i)*&
          10**sspgrid%logfkrpa(vv,nage,nzmet-1)
  ENDDO
  !hot star Teff in kK
  sspgrid%teffarrhot = (/8.0,10.,20.,30./)

  !define central wavelengths of emission lines (in vacuum)
  !these wavelengths come from NIST
  emlines(1)  = 4102.89  ! Hd
  emlines(2)  = 4341.69  ! Hy
  emlines(3)  = 4862.71  ! Hb
  emlines(4)  = 4960.30  ! [OIII]
  emlines(5)  = 5008.24  ! [OIII]
  emlines(6)  = 5203.05  ! [NI]
  emlines(7)  = 6549.86  ! [NII]
  emlines(8)  = 6564.61  ! Ha
  emlines(9)  = 6585.27  ! [NII]
  emlines(10) = 6718.29  ! [SII]
  emlines(11) = 6732.67  ! [SII]

  !read in template error function (computed from SDSS stacks)
  !NB: this hasn't been used in years!
  OPEN(28,FILE=TRIM(SPECFIT_HOME)//'/infiles/temperrfcn.s350',&
       STATUS='OLD',iostat=stat,ACTION='READ')
  DO i=1,nstart-1
     READ(28,*) 
  ENDDO
  DO i=1,nl
     READ(28,*) d1,temperrfcn(i)
  ENDDO
  CLOSE(28)

  !read in the atm transmission function
  OPEN(29,FILE=TRIM(SPECFIT_HOME)//'/infiles/atm_trans.dat',&
       STATUS='OLD',iostat=stat,ACTION='READ')
  DO i=1,ntrans
     READ(29,*) ltrans(i),ftrans(i)
  ENDDO
  CLOSE(29)

  !smooth the trans curve here before interpolation to the main grid
  IF (MAXVAL(data(1:datmax)%ires).GT.0.0) THEN
     strans = linterp(data(1:datmax)%lam,data(1:datmax)%ires,ltrans)
     strans = MIN(MAX(strans,0.0),MAXVAL(data(1:datmax)%ires))
  ELSE
     !force the instrumental resolution to 100 km/s if not explicitly set
     !only done here b/c the transmission function is tabulated at high res
     strans = 100.0 
  ENDIF
  CALL VELBROAD(ltrans,ftrans,sig0,lamlo,lamhi,strans)
  
  !interpolate onto the main wavelength grid.  Force transmission
  !to be 1.0 outside of the bounds of the tabulated function
  sspgrid%atm_trans=1.0
  sspgrid%atm_trans = linterp(ltrans,ftrans,lam)
  DO i=1,nl
     IF (lam(i).LT.ltrans(1).OR.lam(i).GT.ltrans(ntrans)) &
        sspgrid%atm_trans(i) = 1.0
  ENDDO

 
  !smooth the models to the input instrumental resolution
  IF (MAXVAL(data(1:datmax)%ires).GT.10.0) THEN

     !the interpolation here is a massive extrapolation beyond the range
     !of the data.  This *should't* matter since we dont use the model
     !beyond the range of the data, but I should double check this at some point
     smooth = linterp(data(1:datmax)%lam,data(1:datmax)%ires,sspgrid%lam)
     smooth = MIN(MAX(smooth,0.0),MAXVAL(data(1:datmax)%ires))

     DO k=1,nzmet
        DO j=1,nage_rfcn
           CALL VELBROAD(lam,sspgrid%solar(:,j,k),sig0,lamlo,lamhi,smooth)
           CALL VELBROAD(lam,sspgrid%nap(:,j,k),sig0,lamlo,lamhi,smooth)
           CALL VELBROAD(lam,sspgrid%nam(:,j,k),sig0,lamlo,lamhi,smooth)
           CALL VELBROAD(lam,sspgrid%cap(:,j,k),sig0,lamlo,lamhi,smooth)
           CALL VELBROAD(lam,sspgrid%cam(:,j,k),sig0,lamlo,lamhi,smooth)
           CALL VELBROAD(lam,sspgrid%fep(:,j,k),sig0,lamlo,lamhi,smooth)
           CALL VELBROAD(lam,sspgrid%fem(:,j,k),sig0,lamlo,lamhi,smooth)
           CALL VELBROAD(lam,sspgrid%cp(:,j,k),sig0,lamlo,lamhi,smooth)
           CALL VELBROAD(lam,sspgrid%cm(:,j,k),sig0,lamlo,lamhi,smooth)
           CALL VELBROAD(lam,sspgrid%np(:,j,k),sig0,lamlo,lamhi,smooth)
           CALL VELBROAD(lam,sspgrid%nm(:,j,k),sig0,lamlo,lamhi,smooth)
           CALL VELBROAD(lam,sspgrid%ap(:,j,k),sig0,lamlo,lamhi,smooth)
           CALL VELBROAD(lam,sspgrid%tip(:,j,k),sig0,lamlo,lamhi,smooth)
           CALL VELBROAD(lam,sspgrid%tim(:,j,k),sig0,lamlo,lamhi,smooth)
           CALL VELBROAD(lam,sspgrid%mgp(:,j,k),sig0,lamlo,lamhi,smooth)
           CALL VELBROAD(lam,sspgrid%mgm(:,j,k),sig0,lamlo,lamhi,smooth)
           CALL VELBROAD(lam,sspgrid%sip(:,j,k),sig0,lamlo,lamhi,smooth)
           CALL VELBROAD(lam,sspgrid%sim(:,j,k),sig0,lamlo,lamhi,smooth)
           CALL VELBROAD(lam,sspgrid%teffp(:,j,k),sig0,lamlo,lamhi,smooth)
           CALL VELBROAD(lam,sspgrid%teffm(:,j,k),sig0,lamlo,lamhi,smooth)
           CALL VELBROAD(lam,sspgrid%crp(:,j,k),sig0,lamlo,lamhi,smooth)
           CALL VELBROAD(lam,sspgrid%mnp(:,j,k),sig0,lamlo,lamhi,smooth)
           CALL VELBROAD(lam,sspgrid%bap(:,j,k),sig0,lamlo,lamhi,smooth)
           CALL VELBROAD(lam,sspgrid%bam(:,j,k),sig0,lamlo,lamhi,smooth)
           CALL VELBROAD(lam,sspgrid%nip(:,j,k),sig0,lamlo,lamhi,smooth)
           CALL VELBROAD(lam,sspgrid%cop(:,j,k),sig0,lamlo,lamhi,smooth)
           CALL VELBROAD(lam,sspgrid%eup(:,j,k),sig0,lamlo,lamhi,smooth)
           CALL VELBROAD(lam,sspgrid%srp(:,j,k),sig0,lamlo,lamhi,smooth)
           CALL VELBROAD(lam,sspgrid%kp(:,j,k),sig0,lamlo,lamhi,smooth)
           CALL VELBROAD(lam,sspgrid%vp(:,j,k),sig0,lamlo,lamhi,smooth)
           CALL VELBROAD(lam,sspgrid%cup(:,j,k),sig0,lamlo,lamhi,smooth)
           CALL VELBROAD(lam,sspgrid%nap6(:,j,k),sig0,lamlo,lamhi,smooth)
           CALL VELBROAD(lam,sspgrid%nap9(:,j,k),sig0,lamlo,lamhi,smooth)
        ENDDO
     ENDDO

     DO i=1,nhot
        CALL VELBROAD(lam,sspgrid%hotspec(:,i),sig0,lamlo,lamhi,smooth)
     ENDDO

     CALL VELBROAD(lam,sspgrid%m7g,sig0,lamlo,lamhi,smooth)

     DO z=1,nzmet
        DO t=1,nage
           DO j=1,nimf
              DO k=1,nimf
                 CALL VELBROAD(lam,sspgrid%imf(:,j,k,t,z),sig0,lamlo,lamhi,smooth)
              ENDDO
           ENDDO
        ENDDO
     ENDDO

     sspgrid%logfkrpa = 10**sspgrid%logfkrpa
     DO j=1,nzmet
        DO i=1,nage
           CALL VELBROAD(lam,sspgrid%logfkrpa(:,i,j),sig0,lamlo,lamhi,smooth)
        ENDDO
     ENDDO
     sspgrid%logfkrpa = LOG10(sspgrid%logfkrpa+tiny_number)

  ENDIF

  !locate where lam=7000A
  lam7 = locate(lam,7000.d0)


END SUBROUTINE SETUP
