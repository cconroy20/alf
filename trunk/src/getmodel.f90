SUBROUTINE GETMODEL(pos,spec,mw)

  USE sfvars; USE nr, ONLY : locate
  USE sfutils, ONLY : velbroad
  IMPLICIT NONE

  TYPE(PARAMS), INTENT(in) :: pos
  REAL(DP), DIMENSION(nl), INTENT(out) :: spec
  INTEGER, OPTIONAL :: mw
  REAL(DP), DIMENSION(nl) :: tmp
  INTEGER :: vt,vv1,vv2,i
  REAL(DP) :: dt,fy,dx1,dx2,lsig,vz

  !---------------------------------------------------------------!
  !---------------------------------------------------------------!

  !vary age
  vt   = MAX(MIN(locate(sspgrid%agegrid,pos%age),nage-1),1)
  dt   = (pos%age-sspgrid%agegrid(vt))/&
       (sspgrid%agegrid(vt+1)-sspgrid%agegrid(vt))
  spec = dt*sspgrid%logfkrpa(vt+1,:) + (1-dt)*sspgrid%logfkrpa(vt,:)
  spec = 10**spec


  !vary [Fe/H]
  IF (pos%feh.GT.0.0) THEN 
     tmp  = 1 + (sspgrid%fep/sspgrid%solar-1)*pos%feh/0.3
     spec = spec * tmp
  ELSE IF (pos%feh.LT.0.0) THEN     
     tmp  = 1 + (sspgrid%fem/sspgrid%solar-1)*ABS(pos%feh)/0.3
     spec = spec * tmp
  ENDIF

  !vary [O/Fe]
  tmp  = 1 + (sspgrid%ap/sspgrid%solar-1)*pos%afe/0.2
  spec = spec * tmp

  !vary [C/Fe]
  IF (pos%cfe.GT.0.0) THEN 
     tmp  = 1 + (sspgrid%cp/sspgrid%solar-1)*pos%cfe/0.15
     spec = spec * tmp
  ELSE IF (pos%cfe.LT.0.0) THEN     
     tmp  = 1 + (sspgrid%cm/sspgrid%solar-1)*ABS(pos%cfe)/0.15
     spec = spec * tmp
  ENDIF

  !vary [N/Fe]
  IF (pos%nfe.GT.0.0) THEN 
     tmp  = 1 + (sspgrid%np/sspgrid%solar-1)*pos%nfe/0.3
     spec = spec * tmp
  ELSE IF (pos%nfe.LT.0.0) THEN     
     tmp  = 1 + (sspgrid%nm/sspgrid%solar-1)*ABS(pos%nfe)/0.3
     spec = spec * tmp
  ENDIF

  !vary [Na/Fe]
  IF (pos%nafe.GT.0.0.AND.pos%nafe.LT.0.3) THEN 
     tmp  = 1 + (sspgrid%nap/sspgrid%solar-1)*pos%nafe/0.3
     spec = spec * tmp
  ELSE IF (pos%nafe.GE.0.3.AND.pos%nafe.LT.0.6) THEN 
     tmp  = sspgrid%nap + (sspgrid%nap6-sspgrid%nap)*(pos%nafe-0.3)/0.3
     spec = spec * tmp/sspgrid%solar
  ELSE IF (pos%nafe.GE.0.6) THEN 
     tmp  = sspgrid%nap6 + (sspgrid%nap9-sspgrid%nap6)*(pos%nafe-0.6)/0.6
     spec = spec * tmp/sspgrid%solar
  ELSE IF (pos%nafe.LT.0.0) THEN     
     tmp  = 1 + (sspgrid%nam/sspgrid%solar-1)*ABS(pos%nafe)/0.3
     spec = spec * tmp
  ENDIF

  !vary [Mg/Fe]
  IF (pos%mgfe.GT.0.0) THEN 
     tmp  = 1 + (sspgrid%mgp/sspgrid%solar-1)*pos%mgfe/0.3
     spec = spec * tmp
  ELSE IF (pos%mgfe.LT.0.0) THEN     
     tmp  = 1 + (sspgrid%mgm/sspgrid%solar-1)*ABS(pos%mgfe)/0.3
     spec = spec * tmp
  ENDIF

  !vary [Si/Fe]
  IF (pos%sife.GT.0.0) THEN 
     tmp  = 1 + (sspgrid%sip/sspgrid%solar-1)*pos%sife/0.3
     spec = spec * tmp
  ELSE IF (pos%sife.LT.0.0) THEN     
     tmp  = 1 + (sspgrid%sim/sspgrid%solar-1)*ABS(pos%sife)/0.3
     spec = spec * tmp
  ENDIF

  !vary [Ca/Fe]
  IF (pos%cafe.GT.0.0) THEN 
     tmp  = 1 + (sspgrid%cap/sspgrid%solar-1)*pos%cafe/0.15
     spec = spec * tmp
  ELSE IF (pos%cafe.LT.0.0) THEN     
     tmp  = 1 + (sspgrid%cam/sspgrid%solar-1)*ABS(pos%cafe)/0.15
     spec = spec * tmp
  ENDIF

  !vary [Ti/Fe]
  IF (pos%tife.GT.0.0) THEN 
     tmp  = 1 + (sspgrid%tip/sspgrid%solar-1)*pos%tife/0.3
     spec = spec * tmp
  ELSE IF (pos%tife.LT.0.0) THEN     
     tmp  = 1 + (sspgrid%tim/sspgrid%solar-1)*ABS(pos%tife)/0.3
     spec = spec * tmp
  ENDIF

  !only include these parameters in the "full" model
  IF (fitsimple.EQ.0) THEN

     !vary teff
     IF (pos%teff.GT.0.0) THEN 
        tmp  = 1 + (sspgrid%teffp/sspgrid%solar-1)*pos%teff/50.
        spec = spec * tmp
     ELSE IF (pos%teff.LT.0.0) THEN     
        tmp  = 1 + (sspgrid%teffm/sspgrid%solar-1)*ABS(pos%teff)/50.
        spec = spec * tmp
     ENDIF
     
     !vary young (3 Gyr) population
     fy = MAX(MIN(10**pos%logfy,1.0),0.0)
     spec = (1-fy)*spec + fy*10**sspgrid%logfkrpa(1,:)
     
     !add a hot star
     vt = MAX(MIN(locate(sspgrid%teffarrhot,pos%hotteff),3),1)
     dt = (pos%hotteff-sspgrid%teffarrhot(vt))/&
          (sspgrid%teffarrhot(vt+1)-sspgrid%teffarrhot(vt))
     fy = MAX(MIN(10**pos%loghot,1.0),0.0)
     tmp  = dt*sspgrid%hotspec(vt+1,:) + (1-dt)*sspgrid%hotspec(vt,:)
     spec = (1-fy)*spec + fy*tmp
     
     !add in an M7 giant
     fy = MAX(MIN(10**pos%logm7g,1.0),0.0)
     spec = (1-fy)*spec + fy*sspgrid%m7g

     !vary [K/Fe]
     tmp  = 1 + (sspgrid%kp/sspgrid%solar-1)*pos%kfe/0.3
     spec = spec * tmp

     !vary [V/Fe]
     tmp  = 1 + (sspgrid%vp/sspgrid%solar-1)*pos%vfe/0.3
     spec = spec * tmp

     !vary [Cr/Fe]
     tmp  = 1 + (sspgrid%crp/sspgrid%solar-1)*pos%crfe/0.3
     spec = spec * tmp

     !vary [Mn/Fe]
     tmp  = 1 + (sspgrid%mnp/sspgrid%solar-1)*pos%mnfe/0.3
     spec = spec * tmp
     
     !vary [Co/Fe]
     tmp  = 1 + (sspgrid%cop/sspgrid%solar-1)*pos%cofe/0.3
     spec = spec * tmp
     
     !vary [Ni/Fe]
     tmp  = 1 + (sspgrid%nip/sspgrid%solar-1)*pos%nife/0.3
     spec = spec * tmp

     !vary [Cu/Fe]
  !   tmp  = 1 + (sspgrid%cup/sspgrid%solar-1)*pos%cufe/0.3
  !   spec = spec * tmp
     
     !vary [Rb/Fe]
     tmp  = 1 + (sspgrid%rbp/sspgrid%solar-1)*pos%rbfe/0.3
     spec = spec * tmp
     
     !vary [Sr/Fe]
     tmp  = 1 + (sspgrid%srp/sspgrid%solar-1)*pos%srfe/0.3
     spec = spec * tmp
     
     !vary [Y/Fe]
     tmp  = 1 + (sspgrid%yp/sspgrid%solar-1)*pos%yfe/0.3
     spec = spec * tmp
     
     !vary [Zr/Fe]
     tmp  = 1 + (sspgrid%zrp/sspgrid%solar-1)*pos%zrfe/0.3
     spec = spec * tmp
     
     !vary [Ba/Fe]
     IF (pos%bafe.GT.0.0) THEN 
        tmp  = 1 + (sspgrid%bap/sspgrid%solar-1)*pos%bafe/0.3
        spec = spec * tmp
     ELSE IF (pos%bafe.LT.0.0) THEN     
        tmp  = 1 + (sspgrid%bam/sspgrid%solar-1)*ABS(pos%bafe)/0.3
        spec = spec * tmp
     ENDIF

     !vary [Eu/Fe]
     tmp  = 1 + (sspgrid%eup/sspgrid%solar-1)*pos%eufe/0.3
     spec = spec * tmp
     
     IF (.NOT.PRESENT(mw)) THEN
        
        !vary IMF
        vv1 = MAX(MIN(locate(sspgrid%imfx,pos%imf1),nimf-1),1)
        dx1 = (pos%imf1-sspgrid%imfx(vv1))/(sspgrid%imfx(vv1+1)-sspgrid%imfx(vv1))
        dx1 = MAX(MIN(dx1,1.0),-1.0)
        vv2 = MAX(MIN(locate(sspgrid%imfx,pos%imf2),nimf-1),1)
        dx2 = (pos%imf2-sspgrid%imfx(vv2))/(sspgrid%imfx(vv2+1)-sspgrid%imfx(vv2))
        dx2 = MAX(MIN(dx2,1.0),-1.0)
        tmp = (1-dx1)*(1-dx2)*sspgrid%imf(vv1,vv2,:)+&
             dx1*(1-dx2)*sspgrid%imf(vv1+1,vv2,:)+&
             (1-dx1)*dx2*sspgrid%imf(vv1,vv2+1,:)+&
             dx1*dx2*sspgrid%imf(vv1+1,vv2+1,:)
        tmp = tmp/sspgrid%imf(i13,i23,:)
        !turn off IMF sensitivity at lambda<7000A
        !wh = where(la LT 7E3)
        !tmp[wh] = 1.0
        spec = spec * tmp

     ENDIF

  ENDIF

  !add emission lines
  DO i=1,neml
     !allow the em lines to be offset in velocity from the continuum
     !NB: velz2 is a *relative* shift between continuum and lines
     vz   = emlines(i) / (1+pos%velz2/clight*1E5)
     lsig = MAX(vz*pos%sigma2/clight*1E5,1.0)
     spec = spec + 10**pos%logemnorm(i) * &
          EXP(-(sspgrid%lam-vz)**2/lsig**2/2.0)
  ENDDO

  !velocity broaden the model
  IF (pos%sigma.GT.20.0) &
       CALL VELBROAD(sspgrid%lam,spec,pos%sigma)

  IF (apply_temperrfcn.EQ.1) THEN
     spec = spec / temperrfcn
  ENDIF


END SUBROUTINE GETMODEL
