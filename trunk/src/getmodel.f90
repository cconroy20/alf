SUBROUTINE GETMODEL(pos,spec,mw)

  !routine to produce a model spectrum (spec) for an input 
  !set of parameters (pos).  The optional flag 'mw' is used
  !to force the IMF to be of the MW (Kroupa 2001) form

  USE sfvars; USE nr, ONLY : locate
  USE sfutils, ONLY : velbroad, add_response
  IMPLICIT NONE

  TYPE(PARAMS), INTENT(in) :: pos
  REAL(DP), DIMENSION(nl), INTENT(out) :: spec
  INTEGER, OPTIONAL :: mw
  REAL(DP), DIMENSION(nl) :: tmp,tmpr
  INTEGER :: vt,vv1,vv2,i,vr
  REAL(DP) :: dt,fy,dx1,dx2,lsig,vz,dr

  !---------------------------------------------------------------!
  !---------------------------------------------------------------!

  !vary age of the empirical SSPs
  vt   = MAX(MIN(locate(sspgrid%logagegrid,pos%logage),nage-1),1)
  dt   = (pos%logage-sspgrid%logagegrid(vt))/&
       (sspgrid%logagegrid(vt+1)-sspgrid%logagegrid(vt))
  dt   = MAX(dt,-0.25) !no extrapolation younger than 0.5 Gyr
  spec = dt*sspgrid%logfkrpa(vt+1,:) + (1-dt)*sspgrid%logfkrpa(vt,:)
  spec = 10**spec

  !vary age in the response functions
  IF (use_age_dep_resp_fcns.EQ.0) THEN
     !force the use of the 13 Gyr response fcn
     vr = nage_rfcn-1
     dr = 1.0
  ELSE
     vr = MAX(MIN(locate(sspgrid%logagegrid_rfcn,pos%logage),nage_rfcn-1),1)
     dr = (pos%logage-sspgrid%logagegrid_rfcn(vr))/&
          (sspgrid%logagegrid_rfcn(vr+1)-sspgrid%logagegrid_rfcn(vr))
     dr = MAX(MIN(dr,1.0),0.0)
  ENDIF

  !vary [Fe/H]
  CALL ADD_RESPONSE(spec,pos%feh,0.3,dr,vr,sspgrid%solar,sspgrid%fep,sspgrid%fem)
  !vary [O/H]
  CALL ADD_RESPONSE(spec,pos%ah,0.3,dr,vr,sspgrid%solar,sspgrid%ap)
  !vary [C/H]
  CALL ADD_RESPONSE(spec,pos%ch,0.15,dr,vr,sspgrid%solar,sspgrid%cp,sspgrid%cm)
  !vary [N/H]
  CALL ADD_RESPONSE(spec,pos%nh,0.3,dr,vr,sspgrid%solar,sspgrid%np,sspgrid%nm)
  !vary [Mg/H]
  CALL ADD_RESPONSE(spec,pos%mgh,0.3,dr,vr,sspgrid%solar,sspgrid%mgp,sspgrid%mgm)
  !vary [Si/H]
  CALL ADD_RESPONSE(spec,pos%sih,0.3,dr,vr,sspgrid%solar,sspgrid%sip,sspgrid%sim)
  !vary [Ca/H]
  CALL ADD_RESPONSE(spec,pos%cah,0.3,dr,vr,sspgrid%solar,sspgrid%cap,sspgrid%cam)
  !vary [Ti/H]
  CALL ADD_RESPONSE(spec,pos%tih,0.3,dr,vr,sspgrid%solar,sspgrid%tip,sspgrid%tim)

  !vary [Na/H] (special case)
  IF (pos%nah.GT.0.0.AND.pos%nah.LT.0.3) THEN
     tmpr = dr*sspgrid%nap(vr+1,:)/sspgrid%solar(vr+1,:) + &
          (1-dr)*sspgrid%nap(vr,:)/sspgrid%solar(vr,:)
     tmp  = 1 + (tmpr-1)*pos%nah/0.3
     spec = spec * tmp
  ELSE IF (pos%nah.GE.0.3.AND.pos%nah.LT.0.6) THEN
     tmpr = dr*sspgrid%nap(vr+1,:)/sspgrid%solar(vr+1,:) + &
          (1-dr)*sspgrid%nap(vr,:)/sspgrid%solar(vr,:)
     tmp = dr * (sspgrid%nap6(vr+1,:)-sspgrid%nap(vr+1,:))/sspgrid%solar(vr+1,:) + &
          (1-dr)*(sspgrid%nap6(vr,:)-sspgrid%nap(vr,:))/sspgrid%solar(vr,:)
     spec = spec * ( tmpr + tmp*(pos%nah-0.3)/0.3 )
  ELSE IF (pos%nah.GE.0.6) THEN
     tmpr = dr*sspgrid%nap6(vr+1,:)/sspgrid%solar(vr+1,:) + &
          (1-dr)*sspgrid%nap6(vr,:)/sspgrid%solar(vr,:)
     tmp = dr * (sspgrid%nap9(vr+1,:)-sspgrid%nap6(vr+1,:))/sspgrid%solar(vr+1,:) + &
          (1-dr)*(sspgrid%nap9(vr,:)-sspgrid%nap6(vr,:))/sspgrid%solar(vr,:)
     spec = spec * ( tmpr + tmp*(pos%nah-0.6)/0.6 )
  ELSE IF (pos%nah.LT.0.0) THEN
     tmpr = dr*sspgrid%nam(vr+1,:)/sspgrid%solar(vr+1,:) + &
          (1-dr)*sspgrid%nam(vr,:)/sspgrid%solar(vr,:)
     tmp  = 1 + (tmpr-1)*ABS(pos%nah)/0.3
     spec = spec * tmp
  ENDIF


  !only include these parameters in the "full" model
  IF (fitsimple.EQ.0) THEN

     !vary Teff (special case - force use of the 13 Gyr model)
     CALL ADD_RESPONSE(spec,pos%teff,50.,1.d0,nage_rfcn-1,sspgrid%solar,&
          sspgrid%teffp,sspgrid%teffm)
     
     !vary young (3 Gyr) population
     fy = MAX(MIN(10**pos%logfy,1.0),0.0)
     spec = (1-fy)*spec + fy*10**sspgrid%logfkrpa(2,:)
     
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

     !vary [K/H]
     CALL ADD_RESPONSE(spec,pos%kh,0.3,dr,vr,sspgrid%solar,sspgrid%kp)
     !vary [V/H]
     CALL ADD_RESPONSE(spec,pos%vh,0.3,dr,vr,sspgrid%solar,sspgrid%vp)
     !vary [Cr/H]
     CALL ADD_RESPONSE(spec,pos%crh,0.3,dr,vr,sspgrid%solar,sspgrid%crp)
     !vary [Mn/H]
     CALL ADD_RESPONSE(spec,pos%mnh,0.3,dr,vr,sspgrid%solar,sspgrid%mnp)
     !vary [Co/H]
     CALL ADD_RESPONSE(spec,pos%coh,0.3,dr,vr,sspgrid%solar,sspgrid%cop)
     !vary [Ni/H]
     CALL ADD_RESPONSE(spec,pos%nih,0.3,dr,vr,sspgrid%solar,sspgrid%nip)
     !vary [Cu/H]
     CALL ADD_RESPONSE(spec,pos%cuh,0.3,dr,vr,sspgrid%solar,sspgrid%cup)
     !vary [Rb/H]
     CALL ADD_RESPONSE(spec,pos%rbh,0.3,dr,vr,sspgrid%solar,sspgrid%rbp)
     !vary [Sr/H]
     CALL ADD_RESPONSE(spec,pos%srh,0.3,dr,vr,sspgrid%solar,sspgrid%srp)
     !vary [Y/H]
     CALL ADD_RESPONSE(spec,pos%yh,0.3,dr,vr,sspgrid%solar,sspgrid%yp)
     !vary [Zr/H]
     CALL ADD_RESPONSE(spec,pos%zrh,0.3,dr,vr,sspgrid%solar,sspgrid%zrp)
     !vary [Ba/H]
     CALL ADD_RESPONSE(spec,pos%bah,0.3,dr,vr,sspgrid%solar,sspgrid%bap,sspgrid%bam)
     !vary [Eu/H]
     CALL ADD_RESPONSE(spec,pos%euh,0.3,dr,vr,sspgrid%solar,sspgrid%eup)

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

  IF (maskem.EQ.0) THEN
     !add emission lines
     DO i=1,neml
        !allow the em lines to be offset in velocity from the continuum
        !NB: velz2 is a *relative* shift between continuum and lines
        vz   = emlines(i) / (1+pos%velz2/clight*1E5)
        lsig = MAX(vz*pos%sigma2/clight*1E5,1.0)
        spec = spec + 10**pos%logemnorm(i) * &
             EXP(-(sspgrid%lam-vz)**2/lsig**2/2.0)
     ENDDO
  ENDIF

  !velocity broaden the model
  IF (pos%sigma.GT.20.0) &
       CALL VELBROAD(sspgrid%lam,spec,pos%sigma,MINVAL(l1),MAXVAL(l2))

  IF (apply_temperrfcn.EQ.1) THEN
     spec = spec / temperrfcn
  ENDIF


END SUBROUTINE GETMODEL
