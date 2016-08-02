SUBROUTINE SETUP()

  !read in and set up all the arrays

  USE alf_vars; USE nr, ONLY : locate; USE alf_utils, ONLY : linterp,velbroad
  IMPLICIT NONE
  
  REAL(DP) :: d1,l1um=1E4,t13=1.3,t23=2.3,sig0=99.,lamlo,lamhi
  REAL(DP), DIMENSION(nimf_full*nimf_full) :: tmp
  REAL(DP), DIMENSION(nl) :: dumi,smooth=0.0,lam
  INTEGER :: stat,i,vv,j,k,t,z,ii,shift=100,m
  INTEGER, PARAMETER :: ntrans=22800
  REAL(DP), DIMENSION(ntrans) :: ltrans,ftrans_h2o,ftrans_o2,strans
  CHARACTER(4), DIMENSION(nzmet) :: charz
  CHARACTER(4), DIMENSION(nmcut) :: charm
  CHARACTER(5), DIMENSION(nage)  :: chart
  CHARACTER(3), DIMENSION(nage_rfcn) :: chart2
  CHARACTER(20) :: imfstr

  !---------------------------------------------------------------!
  !---------------------------------------------------------------!

  CALL GETENV('ALF_HOME',ALF_HOME)
  IF (TRIM(ALF_HOME).EQ.'') THEN
     WRITE(*,*) 'ALF ERROR: ALF_HOME environment variable not set!'
     STOP
  ENDIF

  charz  = (/'m1.5','m1.0','m0.5','p0.0','p0.2'/)
  charm  = (/'0.08','0.10','0.15','0.20','0.25','0.30','0.35','0.40'/)
  chart  = (/'t01.0','t03.0','t05.0','t07.0','t09.0','t11.0','t13.5'/)
  chart2 = (/'t01','t03','t05','t09','t13'/)

  !this is the IMF flag for the "main" models
  IF (imf_type.EQ.2) THEN
     WRITE(*,*) 'ERROR: imf_type=2 is not currently supported'
     STOP
     imfstr = 'varymcut_varyx'     
  ELSE
     imfstr = 'varydoublex'
  ENDIF

  sspgrid%logssp  = tiny_number
  sspgrid%logsspm = tiny_number

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
  OPEN(15,FILE=TRIM(ALF_HOME)//'/infiles/filters.dat',&
       STATUS='OLD',iostat=stat,ACTION='READ')
  IF (stat.NE.0) THEN
     WRITE(*,*) 'SETUP ERROR: filter curves not found'
     STOP
  ENDIF
  DO i=1,nstart-1
     READ(15,*) 
  ENDDO
  DO i=1,nl
     READ(15,*) d1,filters(i,1),filters(i,2),filters(i,3) !r,i,K filters
  ENDDO
  CLOSE(15)


  !-------------------------------------------------------------------------!
  !-----------------read in the theoretical response functions--------------!
  !-------------------------------------------------------------------------!

  DO k=1,nzmet

     DO j=1,nage_rfcn

        OPEN(20,FILE=TRIM(ALF_HOME)//'/infiles/atlas_ssp_'//&
             chart2(j)//'_Z'//charz(k)//'.abund.'//atlas_imf//'.s100',&
             STATUS='OLD',iostat=stat,ACTION='READ')
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
           READ(20,*) sspgrid%lam(i),sspgrid%solar(i,j,k),&
                sspgrid%nap(i,j,k),sspgrid%nam(i,j,k),&
                sspgrid%cap(i,j,k),sspgrid%cam(i,j,k),sspgrid%fep(i,j,k),&
                sspgrid%fem(i,j,k),sspgrid%cp(i,j,k),sspgrid%cm(i,j,k),d1,&
                sspgrid%np(i,j,k),sspgrid%nm(i,j,k),sspgrid%ap(i,j,k),&
                sspgrid%tip(i,j,k),sspgrid%tim(i,j,k),sspgrid%mgp(i,j,k),&
                sspgrid%mgm(i,j,k),sspgrid%sip(i,j,k),sspgrid%sim(i,j,k),&
                sspgrid%teffp(i,j,k),sspgrid%teffm(i,j,k),sspgrid%crp(i,j,k),&
                sspgrid%mnp(i,j,k),sspgrid%bap(i,j,k),sspgrid%bam(i,j,k),&
                sspgrid%nip(i,j,k),sspgrid%cop(i,j,k),sspgrid%eup(i,j,k),&
                sspgrid%srp(i,j,k),sspgrid%kp(i,j,k),sspgrid%vp(i,j,k),&
                sspgrid%cup(i,j,k),sspgrid%nap6(i,j,k),sspgrid%nap9(i,j,k)
        ENDDO
        CLOSE(20)

     ENDDO
  ENDDO

  lam = sspgrid%lam
  sspgrid%logagegrid_rfcn = LOG10((/1.0,3.0,5.0,9.0,13.0/))


  !create fake response functions by shifting
  !the wavelengths by n pixels
  IF (fake_response.EQ.1) THEN

     WRITE(*,*) 'WARNING: this option has not been tested in a long time!!'

     DO i=1,nage_rfcn
        DO k=1,nzmet

           dumi = sspgrid%crp(:,i,k) / sspgrid%solar(:,i,k)
           sspgrid%crp(:,i,k) = 1.0
           sspgrid%crp(shift:nl-1,i,k) = dumi(1:nl-shift)
           sspgrid%crp(:,i,k) = sspgrid%crp(:,i,k) * sspgrid%solar(:,i,k)
     
           dumi = sspgrid%mnp(:,i,k) / sspgrid%solar(:,i,k)
           sspgrid%mnp = 1.0
           sspgrid%mnp(shift:nl-1,i,k) = dumi(1:nl-shift)
           sspgrid%mnp(:,i,k) = sspgrid%mnp(:,i,k) * sspgrid%solar(:,i,k)
           
           dumi = sspgrid%cop(:,i,k) / sspgrid%solar(:,i,k)
           sspgrid%cop = 1.0
           sspgrid%cop(shift:nl-1,i,k) = dumi(1:nl-shift)
           sspgrid%cop(:,i,k) = sspgrid%cop(:,i,k) * sspgrid%solar(:,i,k)
           
        ENDDO
     ENDDO

  ENDIF

  !-------------------------------------------------------------------------!
  !-----read in empirical spectra as a function of age and metallicity------!
  !-------------------------------------------------------------------------!

  sspgrid%logagegrid = LOG10((/1.0,3.0,5.0,7.0,9.0,11.0,13.5/))
  sspgrid%logzgrid   = (/-1.5,-1.0,-0.5,0.0,0.25/)
  sspgrid%logzgrid2  = (/-0.5,0.0,0.25/)

  !read in two parameter IMF models
  DO z=1,nzmet
     DO t=1,nage
        OPEN(22,FILE=TRIM(ALF_HOME)//'/infiles/VCJ_v1_mcut0.08_'//&
             chart(t)//'_Z'//charz(z)//'.ssp.imf_'//TRIM(imfstr)//'.s100',&
             STATUS='OLD',iostat=stat,ACTION='READ')
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
           DO j=1,nimf+nimfoff
              DO k=1,nimf+nimfoff
                 IF (k.GT.nimfoff.AND.j.GT.nimfoff) THEN
                    sspgrid%logssp(i,j-nimfoff,k-nimfoff,t,z) = tmp(ii)
                 ENDIF
                 ii=ii+1
              ENDDO
           ENDDO
        ENDDO
        CLOSE(22)
     ENDDO
  ENDDO

  !read in 3 parameter IMF models
  IF (imf_type.EQ.3) THEN 
     DO z=1,nzmet3
        DO m=1,nmcut
           DO t=1,nage
              OPEN(22,FILE=TRIM(ALF_HOME)//'/infiles/VCJ_v1_mcut'//&
                   charm(m)//'_'//chart(t)//'_Z'//charz(z+2)//&
                   '.ssp.imf_'//TRIM(imfstr)//'.s100',STATUS='OLD',&
                   iostat=stat,ACTION='READ')
              IF (stat.NE.0) THEN
                 WRITE(*,*) 'SETUP ERROR: IMF 3-part models not found'
                 STOP
              ENDIF
              DO i=1,nstart-1
                 READ(22,*) 
              ENDDO
              DO i=1,nl
                 READ(22,*) d1,tmp
                 ii=1
                 DO j=1,nimf+nimfoff
                    DO k=1,nimf+nimfoff
                       IF (k.GT.nimfoff.AND.j.GT.nimfoff) THEN
                          sspgrid%logsspm(i,j-nimfoff,k-nimfoff,t,m,z) = tmp(ii)
                       ENDIF
                       ii=ii+1
                    ENDDO
                 ENDDO
              ENDDO
              CLOSE(22)
           ENDDO
        ENDDO
     ENDDO
  ENDIF

  !read in non-parametric IMF models
  IF (imf_type.EQ.4) THEN 
     DO z=1,nzmet
        DO t=1,nage
           OPEN(22,FILE=TRIM(ALF_HOME)//'/infiles/VCJ_v1_'//&
                chart(t)//'_Z'//charz(z)//'.ssp.imf_nonpara_flat'//&
                '.s100',STATUS='OLD',iostat=stat,ACTION='READ')
           IF (stat.NE.0) THEN
              WRITE(*,*) 'SETUP ERROR: non-param IMF models not found'
              STOP
           ENDIF
           DO i=1,nstart-1
              READ(22,*) 
           ENDDO
           DO i=1,nl
              READ(22,*) d1,sspgrid%sspnp(i,1,t,z),sspgrid%sspnp(i,2,t,z),&
                   sspgrid%sspnp(i,3,t,z),sspgrid%sspnp(i,4,t,z),&
                   sspgrid%sspnp(i,5,t,z),sspgrid%sspnp(i,6,t,z),&
                   sspgrid%sspnp(i,7,t,z),sspgrid%sspnp(i,8,t,z),&
                   sspgrid%sspnp(i,9,t,z)
           ENDDO
        ENDDO
     ENDDO
  ENDIF


  !values of IMF parameters at the grid points
  DO i=1,nimf
     sspgrid%imfx1(i) = 0.5+REAL(i-1+nimfoff)/5.d0 
  ENDDO
  IF (imf_type.EQ.2) THEN
     !sspgrid%imfx2 = (/0.07,0.10,0.15,0.2,0.25,0.3,0.35,0.4,&
     !     0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8/)
     write(*,*) 'ERROR: imf_type=2 no longer supported'
     STOP
  ELSE
     sspgrid%imfx2 = sspgrid%imfx1
     imfr1 = locate(sspgrid%imfx1,t13+1E-3)
     imfr2 = locate(sspgrid%imfx2,t23+1E-3)
  ENDIF
  IF (imf_type.EQ.3) THEN
     sspgrid%imfx3 = (/0.08,0.10,0.15,0.2,0.25,0.3,0.35,0.4/)
     imfr3 = locate(sspgrid%imfx3,imflo+1E-3)
  ENDIF

  !-------------------------------------------------------------------------!
  !------------------------set up nuisance features-------------------------!
  !-------------------------------------------------------------------------!

  !read in hot stars 
  OPEN(24,FILE=TRIM(ALF_HOME)//'/infiles/ap00t08000g4.00at12.spec.s100',&
       STATUS='OLD',iostat=stat,ACTION='READ')
  DO i=1,nstart-1
     READ(24,*) 
  ENDDO
  DO i=1,nl
     READ(24,*) d1,sspgrid%hotspec(i,1)
  ENDDO
  CLOSE(24)
  OPEN(25,FILE=TRIM(ALF_HOME)//'/infiles/ap00t10000g4.00at12.spec.s100',&
       STATUS='OLD',iostat=stat,ACTION='READ')
  DO i=1,nstart-1
     READ(25,*) 
  ENDDO
  DO i=1,nl
     READ(25,*) d1,sspgrid%hotspec(i,2)
  ENDDO
  CLOSE(25)
  OPEN(26,FILE=TRIM(ALF_HOME)//'/infiles/ap00t20000g4.00at12.spec.s100',&
       STATUS='OLD',iostat=stat,ACTION='READ')
  DO i=1,nstart-1
     READ(26,*) 
  ENDDO
  DO i=1,nl
     READ(26,*) d1,sspgrid%hotspec(i,3)
  ENDDO
  CLOSE(26)
  OPEN(27,FILE=TRIM(ALF_HOME)//'/infiles/ap00t30000g4.00at12.spec.s100',&
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
          sspgrid%logssp(vv,imfr1,imfr2,nage,nzmet-1)
  ENDDO
  !hot star Teff in kK
  sspgrid%teffarrhot = (/8.0,10.,20.,30./)

  !read in M7III star, normalized to a 13 Gyr SSP at 1um
  OPEN(23,FILE=TRIM(ALF_HOME)//'/infiles/M7III.spec.s100',&
       STATUS='OLD',iostat=stat,ACTION='READ')
  DO i=1,nstart-1
     READ(23,*) 
  ENDDO
  DO i=1,nl
     READ(23,*) d1,sspgrid%m7g(i)
  ENDDO
  CLOSE(23)
  !normalization provided here as opposed to in external IDL routine on 5/13/16
  sspgrid%m7g = sspgrid%m7g/sspgrid%m7g(vv)*&
       sspgrid%logssp(vv,imfr1,imfr2,nage,nzmet-1)

  !-------------------------------------------------------------------------!
  !-------------------------------------------------------------------------!
  !-------------------------------------------------------------------------!

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

  IF (apply_temperrfcn.EQ.1) THEN

     !read in template error function (computed from SDSS stacks)
     !NB: this hasn't been used in years!
     WRITE(*,*) 'WARNING: this option has not been tested in a long time!!'

     OPEN(28,FILE=TRIM(ALF_HOME)//'/infiles/temperrfcn.s350',&
          STATUS='OLD',iostat=stat,ACTION='READ')
     IF (stat.NE.0) THEN
        WRITE(*,*) 'SETUP ERROR: template error function not found'
        STOP
     ENDIF
     DO i=1,nstart-1
        READ(28,*) 
     ENDDO
     DO i=1,nl
        READ(28,*) d1,temperrfcn(i)
     ENDDO
     CLOSE(28)

  ENDIF

  !-------------------------------------------------------------------------!
  !------------set up the atm transmission function & sky lines-------------!
  !-------------------------------------------------------------------------!

  OPEN(29,FILE=TRIM(ALF_HOME)//'/infiles/atm_trans.dat',&
       STATUS='OLD',iostat=stat,ACTION='READ')
  IF (stat.NE.0) THEN
     WRITE(*,*) 'SETUP ERROR: atm trans function not found'
     STOP
  ENDIF
  DO i=1,ntrans
     READ(29,*) ltrans(i),ftrans_h2o(i),ftrans_o2(i)
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
  !add all the terms in quad, including a floor of 10 km/s
  strans = SQRT(strans**2+smooth_trans**2+10.**2)
  !use the simple version which allows for arrays of arbitrary length
  d1 = velbroad_simple
  velbroad_simple = 1 
  CALL VELBROAD(ltrans,ftrans_h2o,sig0,lamlo,lamhi,strans)
  CALL VELBROAD(ltrans,ftrans_o2,sig0,lamlo,lamhi,strans)
  velbroad_simple = d1
  
  !interpolate onto the main wavelength grid.  Force transmission
  !to be 1.0 outside of the bounds of the tabulated function
  sspgrid%atm_trans_h2o = 1.0
  sspgrid%atm_trans_o2  = 1.0
  sspgrid%atm_trans_h2o = linterp(ltrans,ftrans_h2o,lam)
  sspgrid%atm_trans_o2  = linterp(ltrans,ftrans_o2,lam)
  DO i=1,nl
     IF (lam(i).LT.ltrans(1).OR.lam(i).GT.ltrans(ntrans)) THEN
        sspgrid%atm_trans_h2o(i) = 1.0
        sspgrid%atm_trans_o2(i)  = 1.0
     ENDIF
  ENDDO

  !sky lines
  OPEN(29,FILE=TRIM(ALF_HOME)//'/infiles/radiance_lines.dat',&
       STATUS='OLD',iostat=stat,ACTION='READ')
  IF (stat.NE.0) THEN
     WRITE(*,*) 'SETUP ERROR: sky lines file not found'
     STOP
  ENDIF
  DO i=1,nskylines
     READ(29,*) lsky(i),fsky(i)
  ENDDO
  CLOSE(29)

  !smooth by the instrumental resolution
  !use the simple version which allows for arrays of arbitrary length
  d1 = velbroad_simple
  velbroad_simple = 1 
  CALL VELBROAD(lsky,fsky,sig0,lamlo,lamhi,strans)
  velbroad_simple = d1
  fsky = fsky / MAXVAL(fsky)

  !-------------------------------------------------------------------------!
  !---------smooth the models to the input instrumental resolution----------!
  !-------------------------------------------------------------------------!

  IF (MAXVAL(data(1:datmax)%ires).GT.10.0) THEN

     !the interpolation here is a massive extrapolation beyond the range
     !of the data.  This *should't* matter since we dont use the model
     !beyond the range of the data, but I should double check this at some point
     smooth = linterp(data(1:datmax)%lam,data(1:datmax)%ires,sspgrid%lam)
     smooth = MIN(MAX(smooth,0.0),MAXVAL(data(1:datmax)%ires))

     !smooth the response functions
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

     !smooth the standard two-part power-law IMF models
     DO z=1,nzmet
        DO t=1,nage
           DO j=1,nimf
              DO k=1,nimf
                 CALL VELBROAD(lam,sspgrid%logssp(:,k,j,t,z),sig0,lamlo,lamhi,smooth)
              ENDDO
           ENDDO
        ENDDO
     ENDDO

     !smooth the 3-part IMF models
     IF (imf_type.EQ.3) THEN
        DO z=1,nzmet3
           DO m=1,nmcut
              DO t=1,nage
                 DO j=1,nimf
                    DO k=1,nimf
                       CALL VELBROAD(lam,sspgrid%logsspm(:,k,j,t,m,z),&
                            sig0,lamlo,lamhi,smooth)
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDIF

     !smooth the non-parametric IMF models
     IF (imf_type.EQ.4) THEN
        DO z=1,nzmet
           DO t=1,nage
              DO j=1,nimfnp
                 CALL VELBROAD(lam,sspgrid%sspnp(:,j,t,z),sig0,lamlo,lamhi,smooth)
              ENDDO
           ENDDO
        ENDDO
     ENDIF

  ENDIF

  sspgrid%logssp  = LOG10(sspgrid%logssp+tiny_number)
  sspgrid%logsspm = LOG10(sspgrid%logsspm+tiny_number)

  !locate where lam=7000A
  lam7 = locate(lam,7000.d0)


END SUBROUTINE SETUP
