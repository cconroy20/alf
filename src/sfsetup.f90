SUBROUTINE SFSETUP()

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
     READ(11,*) d1,fil(1,i),fil(2,i),fil(3,i) !r,i,K filters
  ENDDO
  CLOSE(11)

  !read in the ATLAS SSPs
  OPEN(20,FILE=TRIM(SPECFIT_HOME)//'/infiles/atlas_ssp.abund.chab.s100',&
       STATUS='OLD',iostat=stat,ACTION='READ')
  READ(20,*) !burn the header
  READ(20,*)
  DO i=1,nstart-1
     READ(20,*) 
  ENDDO
  DO i=1,nl
     READ(20,*) sspgrid%lam(i),sspgrid%solar(i),sspgrid%nap(i),sspgrid%nam(i),&
          sspgrid%cap(i),sspgrid%cam(i),sspgrid%fep(i),sspgrid%fem(i),&
          sspgrid%cp(i),sspgrid%cm(i),d1,sspgrid%zp(i),sspgrid%zm(i),&
          sspgrid%np(i),sspgrid%nm(i),sspgrid%ap(i),sspgrid%tip(i),sspgrid%tim(i),&
          sspgrid%mgp(i),sspgrid%mgm(i),sspgrid%sip(i),sspgrid%sim(i),&
          sspgrid%hep(i),sspgrid%hem(i),sspgrid%teffp(i),sspgrid%teffm(i),&
          sspgrid%crp(i),sspgrid%mnp(i),sspgrid%bap(i),sspgrid%bam(i),&
          sspgrid%nip(i),sspgrid%cop(i),sspgrid%eup(i),sspgrid%srp(i),&
          sspgrid%kp(i),sspgrid%vp(i),sspgrid%yp(i),sspgrid%zrp(i),sspgrid%rbp(i),&
          sspgrid%cup(i)
  ENDDO
  CLOSE(20)

  IF (1.EQ.0) THEN

     !create fake response functions for Cr, Mn, Ni, and Co, by shifting
     !the wavelengths by n pixels
     shift = 150
     test2 = sspgrid%crp / sspgrid%solar
     sspgrid%crp = 1.0
     sspgrid%crp(shift:nl-1) = test2(1:nl-shift)
     sspgrid%crp = sspgrid%crp * sspgrid%solar
     
     test2 = sspgrid%mnp / sspgrid%solar
     sspgrid%mnp = 1.0
     sspgrid%mnp(shift:nl-1) = test2(1:nl-shift)
     sspgrid%mnp = sspgrid%mnp * sspgrid%solar
     
     test2 = sspgrid%nip / sspgrid%solar
     sspgrid%nip = 1.0
     sspgrid%nip(shift:nl-1) = test2(1:nl-shift)
     sspgrid%nip = sspgrid%nip * sspgrid%solar
     
     test2 = sspgrid%cop / sspgrid%solar
     sspgrid%cop = 1.0
     sspgrid%cop(shift:nl-1) = test2(1:nl-shift)
     sspgrid%cop = sspgrid%cop * sspgrid%solar
     
     test2 = sspgrid%vp / sspgrid%solar
     sspgrid%vp = 1.0
     sspgrid%vp(shift:nl-1) = test2(1:nl-shift)
     sspgrid%vp = sspgrid%vp * sspgrid%solar

  ENDIF

  !read in the extra ATLAS SSPs
  OPEN(20,FILE=TRIM(SPECFIT_HOME)//'/infiles/atlas_ssp.abundex.chab.s100',&
       STATUS='OLD',iostat=stat,ACTION='READ')
  READ(20,*) !burn the header
  READ(20,*)
  DO i=1,nstart-1
     READ(20,*) 
  ENDDO
  DO i=1,nl
     READ(20,*) d1,sspgrid%nap6(i),sspgrid%nap9(i)
  ENDDO
  CLOSE(20)


  !read in empirical spectra as a function of age
  OPEN(21,FILE=TRIM(SPECFIT_HOME)//'/infiles/CvD_chabIMF.ssp.s100',&
       STATUS='OLD',iostat=stat,ACTION='READ')
  DO i=1,nstart-1
     READ(21,*) 
  ENDDO
  DO i=1,nl
     READ(21,*) d1,sspgrid%logfchab(1,i),sspgrid%logfchab(2,i),&
          sspgrid%logfchab(3,i),sspgrid%logfchab(4,i),sspgrid%logfchab(5,i),&
          sspgrid%logfchab(6,i)
  ENDDO
  CLOSE(21)
  sspgrid%logfchab = LOG10(sspgrid%logfchab+tiny_number)
  sspgrid%agegrid  = (/3.0,5.0,7.0,9.0,11.0,13.5/)


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
           sspgrid%imf(j,k,i) = tmp(ii)
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


  !read in M7III star, normalized to an SSP at 1um
  OPEN(22,FILE=TRIM(SPECFIT_HOME)//'/infiles/M7III.spec.s100',&
       STATUS='OLD',iostat=stat,ACTION='READ')
  DO i=1,nstart-1
     READ(22,*) 
  ENDDO
  DO i=1,nl
     READ(22,*) d1,sspgrid%m7g(i)
  ENDDO
  CLOSE(22)
 

  !read in three hot stars 

  OPEN(23,FILE=TRIM(SPECFIT_HOME)//'/infiles/ap00t8000g4.00at12.spec.s100',&
       STATUS='OLD',iostat=stat,ACTION='READ')
  DO i=1,nstart-1
     READ(23,*) 
  ENDDO
  DO i=1,nl
     READ(23,*) d1,sspgrid%hotspec(1,i)
  ENDDO
  CLOSE(23)

  OPEN(23,FILE=TRIM(SPECFIT_HOME)//'/infiles/ap00t10000g4.00at12.spec.s100',&
       STATUS='OLD',iostat=stat,ACTION='READ')
  DO i=1,nstart-1
     READ(23,*) 
  ENDDO
  DO i=1,nl
     READ(23,*) d1,sspgrid%hotspec(2,i)
  ENDDO
  CLOSE(23)

  OPEN(24,FILE=TRIM(SPECFIT_HOME)//'/infiles/ap00t20000g4.00at12.spec.s100',&
       STATUS='OLD',iostat=stat,ACTION='READ')
  DO i=1,nstart-1
     READ(24,*) 
  ENDDO
  DO i=1,nl
     READ(24,*) d1,sspgrid%hotspec(3,i)
  ENDDO
  CLOSE(24)

  OPEN(25,FILE=TRIM(SPECFIT_HOME)//'/infiles/ap00t30000g4.00at12.spec.s100',&
       STATUS='OLD',iostat=stat,ACTION='READ')
  DO i=1,nstart-1
     READ(25,*) 
  ENDDO
  DO i=1,nl
     READ(25,*) d1,sspgrid%hotspec(4,i)
  ENDDO
  CLOSE(25)
 
  !normalize to 5000A
  vv = locate(sspgrid%lam(1:nl),l5000)
  DO i=1,3
     sspgrid%hotspec(i,:) = sspgrid%hotspec(i,:)/sspgrid%hotspec(i,vv)*&
          10**sspgrid%logfchab(6,vv)
  ENDDO
  !hot star Teff in kK
  sspgrid%teffarrhot = (/8.0,10.,20.,30./)

  !logarithmic wavelength grid used in velbroad.f90
  dlstep = (LOG(MAXVAL(sspgrid%lam(1:nl)))-LOG(MINVAL(sspgrid%lam(1:nl))))/nl
  DO i=1,nl
     lnlam(i) = i*dlstep+LOG(MINVAL(sspgrid%lam(1:nl)))
  ENDDO


  !define central wavelengths of emission lines (in vacuum)
  emlines(1)  = 3728.38  ! [OII]  the other line is at 3730.29
  emlines(2)  = 3868.00  ! [NeIII]
  emlines(3)  = 4102.89  ! [Hd]
  emlines(4)  = 4341.69  ! [Hg]
  emlines(5)  = 4862.71  ! [Hb]
  emlines(6)  = 4960.30  ! [OIII]
  emlines(7)  = 5008.24  ! [OIII]
  emlines(8)  = 5201.29  ! [NI]  another line is at 5203.05
  emlines(9)  = 6547.35  ! [NII]
  emlines(10) = 6564.60  ! [Ha]
  emlines(11) = 6584.42  ! [NII]
  emlines(12) = 6718.29  ! [SII]
  emlines(13) = 6734.67  ! [SII]

  
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
