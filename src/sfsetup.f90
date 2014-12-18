SUBROUTINE SFSETUP()

  !read in and set up all the arrays

  USE sfvars; USE nr, ONLY : locate
  IMPLICIT NONE
  
  REAL(DP) :: d1,l5000=5000.0, t13=1.3,t23=2.3
  REAL(DP), DIMENSION(nimf*nimf) :: tmp
  REAL(DP), DIMENSION(nl) :: test2
  INTEGER :: stat,i,vv,j,k,ii,shift

  !---------------------------------------------------------------!
  !---------------------------------------------------------------!

  CALL GETENV('SPECFIT_HOME',SPECFIT_HOME)
  IF (TRIM(SPECFIT_HOME).EQ.'') THEN
     WRITE(*,*) 'SFSETUP ERROR: SPECFIT_HOME environment variable not set!'
     STOP
  ENDIF

  !read in filter transmission curves
  OPEN(11,FILE=TRIM(SPECFIT_HOME)//'/infiles/filters.dat',&
       STATUS='OLD',iostat=stat,ACTION='READ')
  DO i=1,nstart-1
     READ(11,*) 
  ENDDO
  DO i=1,nl
     READ(11,*) d1,filters(i,1),filters(i,2),filters(i,3) !r,i,K filters
  ENDDO
  CLOSE(11)


  !read in the ATLAS SSPs
  DO j=1,nage_rfcn

     IF (j.EQ.1) THEN
        OPEN(20,FILE=TRIM(SPECFIT_HOME)//'/infiles/atlas_ssp_t01.abund.krpa.s100',&
             STATUS='OLD',iostat=stat,ACTION='READ')
     ELSE IF (j.EQ.2) THEN
        OPEN(20,FILE=TRIM(SPECFIT_HOME)//'/infiles/atlas_ssp_t03.abund.krpa.s100',&
             STATUS='OLD',iostat=stat,ACTION='READ')
      ELSE IF (j.EQ.3) THEN
        OPEN(20,FILE=TRIM(SPECFIT_HOME)//'/infiles/atlas_ssp_t05.abund.krpa.s100',&
             STATUS='OLD',iostat=stat,ACTION='READ')
     ELSE IF (j.EQ.4) THEN
        OPEN(20,FILE=TRIM(SPECFIT_HOME)//'/infiles/atlas_ssp_t09.abund.krpa.s100',&
             STATUS='OLD',iostat=stat,ACTION='READ')
     ELSE IF (j.EQ.5) THEN
        OPEN(20,FILE=TRIM(SPECFIT_HOME)//'/infiles/atlas_ssp_t13.abund.krpa.s100',&
             STATUS='OLD',iostat=stat,ACTION='READ')
     ENDIF

     READ(20,*) !burn the header
     READ(20,*)
     DO i=1,nstart-1
        READ(20,*) 
     ENDDO
     DO i=1,nl
        READ(20,*) sspgrid%lam(i),sspgrid%solar(i,j),sspgrid%nap(i,j),&
             sspgrid%nam(i,j),sspgrid%cap(i,j),sspgrid%cam(i,j),sspgrid%fep(i,j),&
             sspgrid%fem(i,j),sspgrid%cp(i,j),sspgrid%cm(i,j),d1,sspgrid%zp(i,j),&
             sspgrid%zm(i,j),sspgrid%np(i,j),sspgrid%nm(i,j),sspgrid%ap(i,j),&
             sspgrid%tip(i,j),sspgrid%tim(i,j),sspgrid%mgp(i,j),sspgrid%mgm(i,j),&
             sspgrid%sip(i,j),sspgrid%sim(i,j),sspgrid%hep(i,j),sspgrid%hem(i,j),&
             sspgrid%teffp(i,j),sspgrid%teffm(i,j),sspgrid%crp(i,j),sspgrid%mnp(i,j),&
             sspgrid%bap(i,j),sspgrid%bam(i,j),sspgrid%nip(i,j),sspgrid%cop(i,j),&
             sspgrid%eup(i,j),sspgrid%srp(i,j),sspgrid%kp(i,j),sspgrid%vp(i,j),&
             sspgrid%yp(i,j),sspgrid%zrp(i,j),sspgrid%rbp(i,j),&
             sspgrid%cup(i,j),sspgrid%nap6(i,j),sspgrid%nap9(i,j)
     ENDDO
     CLOSE(20)

  ENDDO

  sspgrid%logagegrid_rfcn = LOG10((/1.0,3.0,5.0,9.0,13.0/))

  IF (1.EQ.0) THEN

     !create fake response functions for Cr, Mn, and Co, by shifting
     !the wavelengths by n pixels
     !shift = 150
     !test2 = sspgrid%crp / sspgrid%solar
     !sspgrid%crp = 1.0
     !sspgrid%crp(shift:nl-1) = test2(1:nl-shift)
     !sspgrid%crp = sspgrid%crp * sspgrid%solar
     
     !test2 = sspgrid%mnp / sspgrid%solar
     !sspgrid%mnp = 1.0
     !sspgrid%mnp(shift:nl-1) = test2(1:nl-shift)
     !sspgrid%mnp = sspgrid%mnp * sspgrid%solar
  
     !test2 = sspgrid%cop / sspgrid%solar
     !sspgrid%cop = 1.0
     !sspgrid%cop(shift:nl-1) = test2(1:nl-shift)
     !sspgrid%cop = sspgrid%cop * sspgrid%solar
 
  ENDIF


  !read in empirical spectra as a function of age
  OPEN(21,FILE=TRIM(SPECFIT_HOME)//'/infiles/CvD_krpaIMF.ssp.s100',&
       STATUS='OLD',iostat=stat,ACTION='READ')
  DO i=1,nstart-1
     READ(21,*) 
  ENDDO
  DO i=1,nl
     READ(21,*) d1,sspgrid%logfkrpa(i,1),sspgrid%logfkrpa(i,2),&
          sspgrid%logfkrpa(i,3),sspgrid%logfkrpa(i,4),sspgrid%logfkrpa(i,5),&
          sspgrid%logfkrpa(i,6),sspgrid%logfkrpa(i,7)
  ENDDO
  CLOSE(21)
  sspgrid%logfkrpa   = LOG10(sspgrid%logfkrpa+tiny_number)
  sspgrid%logagegrid = LOG10((/1.0,3.0,5.0,7.0,9.0,11.0,13.5/))

  !vary two power-law slopes, from 0.1<M<0.5 and 0.5<M<1.0
  OPEN(26,FILE=TRIM(SPECFIT_HOME)//'/infiles/CvD_t13.5.ssp.'//&
       'imf_varydoublex.s100',STATUS='OLD',iostat=stat,ACTION='READ')
  DO i=1,nstart-1
     READ(26,*) 
  ENDDO
  DO i=1,nl
     READ(26,*) d1,tmp
     ii=1
     DO j=1,nimf
        DO k=1,nimf
           sspgrid%imf(i,j,k) = tmp(ii)
           ii=ii+1
        ENDDO
     ENDDO
  ENDDO
  CLOSE(26)
  !values of IMF slopes at the 35 grid points
  DO i=1,nimf
     sspgrid%imfx(i) = 0.d2+REAL(i-1)/10.d0 
  ENDDO
  i13 = locate(sspgrid%imfx,t13+1E-3)
  i23 = locate(sspgrid%imfx,t23+1E-3)

  !read in M7III star, normalized to a 13 Gyr SSP at 1um
  OPEN(22,FILE=TRIM(SPECFIT_HOME)//'/infiles/M7III.spec.s100',&
       STATUS='OLD',iostat=stat,ACTION='READ')
  DO i=1,nstart-1
     READ(22,*) 
  ENDDO
  DO i=1,nl
     READ(22,*) d1,sspgrid%m7g(i)
  ENDDO
  CLOSE(22)
 

  !read in hot stars 

  OPEN(23,FILE=TRIM(SPECFIT_HOME)//'/infiles/ap00t8000g4.00at12.spec.s100',&
       STATUS='OLD',iostat=stat,ACTION='READ')
  DO i=1,nstart-1
     READ(23,*) 
  ENDDO
  DO i=1,nl
     READ(23,*) d1,sspgrid%hotspec(i,1)
  ENDDO
  CLOSE(23)

  OPEN(23,FILE=TRIM(SPECFIT_HOME)//'/infiles/ap00t10000g4.00at12.spec.s100',&
       STATUS='OLD',iostat=stat,ACTION='READ')
  DO i=1,nstart-1
     READ(23,*) 
  ENDDO
  DO i=1,nl
     READ(23,*) d1,sspgrid%hotspec(i,2)
  ENDDO
  CLOSE(23)

  OPEN(24,FILE=TRIM(SPECFIT_HOME)//'/infiles/ap00t20000g4.00at12.spec.s100',&
       STATUS='OLD',iostat=stat,ACTION='READ')
  DO i=1,nstart-1
     READ(24,*) 
  ENDDO
  DO i=1,nl
     READ(24,*) d1,sspgrid%hotspec(i,3)
  ENDDO
  CLOSE(24)

  OPEN(25,FILE=TRIM(SPECFIT_HOME)//'/infiles/ap00t30000g4.00at12.spec.s100',&
       STATUS='OLD',iostat=stat,ACTION='READ')
  DO i=1,nstart-1
     READ(25,*) 
  ENDDO
  DO i=1,nl
     READ(25,*) d1,sspgrid%hotspec(i,4)
  ENDDO
  CLOSE(25)
 
  !normalize to 5000A
  vv = locate(sspgrid%lam(1:nl),l5000)
  DO i=1,nhot
     sspgrid%hotspec(i,:) = sspgrid%hotspec(i,:)/sspgrid%hotspec(vv,i)*&
          10**sspgrid%logfkrpa(vv,6)
  ENDDO
  !hot star Teff in kK
  sspgrid%teffarrhot = (/8.0,10.,20.,30./)

  !define central wavelengths of emission lines (in vacuum)
  emlines(1)  = 4102.92  ! [Hd]
  emlines(2)  = 4341.69  ! [Hg]
  emlines(3)  = 4862.69  ! [Hb]
  emlines(4)  = 4960.30  ! [OIII]
  emlines(5)  = 5008.24  ! [OIII]
  !emlines(6) = 5201.29  ! [NI] - this line was wrong!
  emlines(6)  = 5203.05  ! [NI]
  emlines(7)  = 6549.84  ! [NII]
  emlines(8)  = 6564.61  ! [Ha]
  emlines(9)  = 6585.23  ! [NII]
  emlines(10) = 6718.32  ! [SII]
  emlines(11) = 6732.71  ! [SII]

  
  !read in template error function (computed from SDSS stacks)
  OPEN(22,FILE=TRIM(SPECFIT_HOME)//'/infiles/temperrfcn.s350',&
       STATUS='OLD',iostat=stat,ACTION='READ')
  DO i=1,nstart-1
     READ(22,*) 
  ENDDO
  DO i=1,nl
     READ(22,*) d1,temperrfcn(i)
  ENDDO
  CLOSE(22)
 

END SUBROUTINE SFSETUP
