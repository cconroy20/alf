SUBROUTINE GETMODEL(pos,spec,mw)

  !routine to produce a model spectrum (spec) for an input 
  !set of parameters (pos).  The optional flag 'mw' is used
  !to force the IMF to be of the MW (Kroupa 2001) form

  USE alf_vars; USE nr, ONLY : locate
  USE alf_utils, ONLY : velbroad, add_response,linterp
  IMPLICIT NONE

  TYPE(PARAMS), INTENT(in) :: pos
  REAL(DP), DIMENSION(nl), INTENT(out) :: spec
  INTEGER, OPTIONAL :: mw
  REAL(DP), DIMENSION(nl) :: tmp,tmpr,yspec
  INTEGER :: vt,vv1,vv2,i,vr,vm
  REAL(DP) :: dt,fy,dx1,dx2,lsig,vz,dr,dm
  REAL(DP), DIMENSION(nl)   :: tmp_ltrans, tmp_ftrans
  REAL(DP), DIMENSION(neml) :: emnormall=1.0
  
  !---------------------------------------------------------------!
  !---------------------------------------------------------------!

  !vary age and metallicity of the empirical SSPs
  vt   = MAX(MIN(locate(sspgrid%logagegrid,pos%logage),nage-1),1)
  dt   = (pos%logage-sspgrid%logagegrid(vt))/&
       (sspgrid%logagegrid(vt+1)-sspgrid%logagegrid(vt))
  dt   = MAX(MIN(dt,1.5),-0.3)  !0.5<age<15 Gyr
  vm   = MAX(MIN(locate(sspgrid%logzgrid,pos%zh),nzmet-1),1)
  dm   = (pos%zh-sspgrid%logzgrid(vm)) / &
       (sspgrid%logzgrid(vm+1)-sspgrid%logzgrid(vm))
  dm   = MAX(MIN(dm,1.5),-1.0) ! -2.0<[Z/H]<0.45
  spec(1:nl_fit) = &
       dt*dm*sspgrid%logfkrpa(1:nl_fit,vt+1,vm+1) + &
       (1-dt)*dm*sspgrid%logfkrpa(1:nl_fit,vt,vm+1) + & 
       dt*(1-dm)*sspgrid%logfkrpa(1:nl_fit,vt+1,vm) + & 
       (1-dt)*(1-dm)*sspgrid%logfkrpa(1:nl_fit,vt,vm)
  spec(1:nl_fit) = 10**spec(1:nl_fit)

  !vary age in the response functions
  IF (use_age_dep_resp_fcns.EQ.0) THEN
     !force the use of the 13 Gyr response fcn
     !vr = nage_rfcn-1
     !dr = 1.0
     !force the use of the 5 Gyr response fcn
     vr = 3
     dr = 0.0
  ELSE
     !should be using mass-weighted age here
     vr = MAX(MIN(locate(sspgrid%logagegrid_rfcn,pos%logage),nage_rfcn-1),1)
     dr = (pos%logage-sspgrid%logagegrid_rfcn(vr))/&
          (sspgrid%logagegrid_rfcn(vr+1)-sspgrid%logagegrid_rfcn(vr))
     dr = MAX(MIN(dr,1.0),0.0)
  ENDIF

  !vary young population - both fraction and age
  !only include these parameters in the "full" model
  IF (fit_type.EQ.0.AND.powell_fitting.EQ.0) THEN
     fy    = MAX(MIN(10**pos%logfy,1.0),0.0)
     vt    = MAX(MIN(locate(sspgrid%logagegrid,pos%fy_logage),nage-1),1)
     dt    = (pos%fy_logage-sspgrid%logagegrid(vt))/&
          (sspgrid%logagegrid(vt+1)-sspgrid%logagegrid(vt))
     dt    = MAX(MIN(dt,1.0),-0.3) !0.5<age<13 Gyr
     yspec(1:nl_fit) = &
          dt*dm*sspgrid%logfkrpa(1:nl_fit,vt+1,vm+1) + &
          (1-dt)*dm*sspgrid%logfkrpa(1:nl_fit,vt,vm+1) + & 
          dt*(1-dm)*sspgrid%logfkrpa(1:nl_fit,vt+1,vm) + & 
          (1-dt)*(1-dm)*sspgrid%logfkrpa(1:nl_fit,vt,vm)
      spec(1:nl_fit)  = (1-fy)*spec(1:nl_fit) + fy*10**yspec(1:nl_fit)
  ENDIF
  
 
  !Only sigma, velz, logage, and [Z/H] are fit when either
  !fitting in Powell mode or "super simple" mode
  IF (powell_fitting.EQ.0.AND.fit_type.NE.2) THEN

     !vary [Fe/H]
     CALL ADD_RESPONSE(spec,pos%feh,0.3,dr,vr,dm,vm,sspgrid%solar,sspgrid%fep,sspgrid%fem)
     !vary [O/H]
     CALL ADD_RESPONSE(spec,pos%ah,0.3,dr,vr,dm,vm,sspgrid%solar,sspgrid%ap)
     !vary [C/H]
     CALL ADD_RESPONSE(spec,pos%ch,0.15,dr,vr,dm,vm,sspgrid%solar,sspgrid%cp,sspgrid%cm)
     !vary [N/H]
     CALL ADD_RESPONSE(spec,pos%nh,0.3,dr,vr,dm,vm,sspgrid%solar,sspgrid%np,sspgrid%nm)
     !vary [Mg/H]
     CALL ADD_RESPONSE(spec,pos%mgh,0.3,dr,vr,dm,vm,sspgrid%solar,sspgrid%mgp,sspgrid%mgm)
     !vary [Si/H]
     CALL ADD_RESPONSE(spec,pos%sih,0.3,dr,vr,dm,vm,sspgrid%solar,sspgrid%sip,sspgrid%sim)
     !vary [Ca/H]
     CALL ADD_RESPONSE(spec,pos%cah,0.3,dr,vr,dm,vm,sspgrid%solar,sspgrid%cap,sspgrid%cam)
     !vary [Ti/H]
     CALL ADD_RESPONSE(spec,pos%tih,0.3,dr,vr,dm,vm,sspgrid%solar,sspgrid%tip,sspgrid%tim)
     
     !vary [Na/H] (special case)
     IF (pos%nah.LT.0.3) THEN

        CALL ADD_RESPONSE(spec,pos%nah,0.3,dr,vr,dm,vm,sspgrid%solar,sspgrid%nap,sspgrid%nam)

     ELSE IF (pos%nah.GE.0.3.AND.pos%nah.LT.0.6) THEN

        tmpr(1:nl_fit) = &
             dr*dm*sspgrid%nap(1:nl_fit,vr+1,vm+1)/sspgrid%solar(1:nl_fit,vr+1,vm+1) + &
             (1-dr)*dm*sspgrid%nap(1:nl_fit,vr,vm+1)/sspgrid%solar(1:nl_fit,vr,vm+1) + &
             dr*(1-dm)*sspgrid%nap(1:nl_fit,vr+1,vm)/sspgrid%solar(1:nl_fit,vr+1,vm) + &
             (1-dr)*(1-dm)*sspgrid%nap(1:nl_fit,vr,vm)/sspgrid%solar(1:nl_fit,vr,vm)
        tmp(1:nl_fit) = &
             dr*dm*(sspgrid%nap6(1:nl_fit,vr+1,vm+1)-sspgrid%nap(1:nl_fit,vr+1,vm+1))/&
             sspgrid%solar(1:nl_fit,vr+1,vm+1)+&
             (1-dr)*dm*(sspgrid%nap6(1:nl_fit,vr,vm+1)-sspgrid%nap(1:nl_fit,vr,vm+1))/&
             sspgrid%solar(1:nl_fit,vr,vm+1)+&
             dr*(1-dm)*(sspgrid%nap6(1:nl_fit,vr+1,vm)-sspgrid%nap(1:nl_fit,vr+1,vm))/&
             sspgrid%solar(1:nl_fit,vr+1,vm)+&
             (1-dr)*(1-dm)*(sspgrid%nap6(1:nl_fit,vr,vm)-sspgrid%nap(1:nl_fit,vr,vm))/&
             sspgrid%solar(1:nl_fit,vr,vm)
 
       spec(1:nl_fit) = spec(1:nl_fit) * (tmpr(1:nl_fit)+tmp(1:nl_fit)*(pos%nah-0.3)/0.3 )

     ELSE IF (pos%nah.GE.0.6) THEN

        tmpr(1:nl_fit) = &
             dr*dm*sspgrid%nap6(1:nl_fit,vr+1,vm+1)/sspgrid%solar(1:nl_fit,vr+1,vm+1) + &
             (1-dr)*dm*sspgrid%nap6(1:nl_fit,vr,vm+1)/sspgrid%solar(1:nl_fit,vr,vm+1) + &
             dr*(1-dm)*sspgrid%nap6(1:nl_fit,vr+1,vm)/sspgrid%solar(1:nl_fit,vr+1,vm) + &
             (1-dr)*(1-dm)*sspgrid%nap6(1:nl_fit,vr,vm)/sspgrid%solar(1:nl_fit,vr,vm)

        tmp(1:nl_fit) = &
             dr*dm*(sspgrid%nap9(1:nl_fit,vr+1,vm+1)-sspgrid%nap6(1:nl_fit,vr+1,vm+1))/&
             sspgrid%solar(1:nl_fit,vr+1,vm+1)+&
             (1-dr)*dm*(sspgrid%nap9(1:nl_fit,vr,vm+1)-sspgrid%nap6(1:nl_fit,vr,vm+1))/&
             sspgrid%solar(1:nl_fit,vr,vm+1)+&
             dr*(1-dm)*(sspgrid%nap9(1:nl_fit,vr+1,vm)-sspgrid%nap6(1:nl_fit,vr+1,vm))/&
             sspgrid%solar(1:nl_fit,vr+1,vm)+&
             (1-dr)*(1-dm)*(sspgrid%nap9(1:nl_fit,vr,vm)-sspgrid%nap6(1:nl_fit,vr,vm))/&
             sspgrid%solar(1:nl_fit,vr,vm)

        spec(1:nl_fit) = spec(1:nl_fit) * (tmpr(1:nl_fit)+tmp*(pos%nah-0.6)/0.6 )

     ENDIF
     
  ENDIF

  !only include these parameters in the "full" model
  IF (fit_type.EQ.0.AND.powell_fitting.EQ.0) THEN

     !vary Teff (special case - force use of the 13 Gyr model)
     CALL ADD_RESPONSE(spec,pos%teff,50.,1.d0,nage_rfcn-1,dm,vm,sspgrid%solar,&
          sspgrid%teffp,sspgrid%teffm)

     !add a hot star
     vt = MAX(MIN(locate(sspgrid%teffarrhot,pos%hotteff),nhot-1),1)
     dt = (pos%hotteff-sspgrid%teffarrhot(vt))/&
          (sspgrid%teffarrhot(vt+1)-sspgrid%teffarrhot(vt))
     fy = MAX(MIN(10**pos%loghot,1.0),0.0)
     tmp(1:nl_fit)  = dt*sspgrid%hotspec(1:nl_fit,vt+1) + &
          (1-dt)*sspgrid%hotspec(1:nl_fit,vt)
     spec(1:nl_fit) = (1-fy)*spec(1:nl_fit) + fy*tmp(1:nl_fit)
     
     !add in an M7 giant
     fy = MAX(MIN(10**pos%logm7g,1.0),0.0)
     spec(1:nl_fit) = (1-fy)*spec(1:nl_fit) + fy*sspgrid%m7g(1:nl_fit)

     !vary [K/H]
     CALL ADD_RESPONSE(spec,pos%kh,0.3,dr,vr,dm,vm,sspgrid%solar,sspgrid%kp)
     !vary [V/H]
     CALL ADD_RESPONSE(spec,pos%vh,0.3,dr,vr,dm,vm,sspgrid%solar,sspgrid%vp)
     !vary [Cr/H]
     CALL ADD_RESPONSE(spec,pos%crh,0.3,dr,vr,dm,vm,sspgrid%solar,sspgrid%crp)
     !vary [Mn/H]
     CALL ADD_RESPONSE(spec,pos%mnh,0.3,dr,vr,dm,vm,sspgrid%solar,sspgrid%mnp)
     !vary [Co/H]
     CALL ADD_RESPONSE(spec,pos%coh,0.3,dr,vr,dm,vm,sspgrid%solar,sspgrid%cop)
     !vary [Ni/H]
     CALL ADD_RESPONSE(spec,pos%nih,0.3,dr,vr,dm,vm,sspgrid%solar,sspgrid%nip)
     !vary [Cu/H]
     CALL ADD_RESPONSE(spec,pos%cuh,0.3,dr,vr,dm,vm,sspgrid%solar,sspgrid%cup)
     !vary [Sr/H]
     CALL ADD_RESPONSE(spec,pos%srh,0.3,dr,vr,dm,vm,sspgrid%solar,sspgrid%srp)
     !vary [Ba/H]
     CALL ADD_RESPONSE(spec,pos%bah,0.3,dr,vr,dm,vm,sspgrid%solar,sspgrid%bap,sspgrid%bam)
     !vary [Eu/H]
     CALL ADD_RESPONSE(spec,pos%euh,0.3,dr,vr,dm,vm,sspgrid%solar,sspgrid%eup)

     !vary IMF
     !Note that this part of the model does not depend on age nor Z!
     IF (mwimf.EQ.0.AND..NOT.PRESENT(mw)) THEN
        vv1 = MAX(MIN(locate(sspgrid%imfx,pos%imf1),nimf-1),1)
        dx1 = (pos%imf1-sspgrid%imfx(vv1))/(sspgrid%imfx(vv1+1)-sspgrid%imfx(vv1))
        dx1 = MAX(MIN(dx1,1.0),-1.0)
        IF (fit_oneimf.EQ.0) THEN
           vv2 = MAX(MIN(locate(sspgrid%imfx,pos%imf2),nimf-1),1)
           dx2 = (pos%imf2-sspgrid%imfx(vv2))/(sspgrid%imfx(vv2+1)-sspgrid%imfx(vv2))
           dx2 = MAX(MIN(dx2,1.0),-1.0)
        ELSE
           vv2 = vv1
           dx2 = dx1
        ENDIF
        tmp(1:nl_fit) = (1-dx1)*(1-dx2)*sspgrid%imf(1:nl_fit,vv1,vv2)+&
             dx1*(1-dx2)*sspgrid%imf(1:nl_fit,vv1+1,vv2)+&
             (1-dx1)*dx2*sspgrid%imf(1:nl_fit,vv1,vv2+1)+&
             dx1*dx2*sspgrid%imf(1:nl_fit,vv1+1,vv2+1)
        tmp(1:nl_fit) = tmp(1:nl_fit)/sspgrid%imf(1:nl_fit,i13,i23)
        !turn off IMF sensitivity at lambda<7000A
        IF (blueimf_off.EQ.1) THEN
           tmp(1:lam7) = 1.0
        ENDIF
        spec(1:nl_fit) = spec(1:nl_fit) * tmp(1:nl_fit)
     ENDIF

     !add emission lines
     IF (maskem.EQ.0) THEN

        !these line ratios come from Cloudy or Osterbrock (Table 4.2)
        emnormall(1)  = 10**pos%logemline_h / 11.21  ! Hy
        emnormall(2)  = 10**pos%logemline_h / 6.16   ! Hd
        emnormall(3)  = 10**pos%logemline_h / 2.87   ! Hb
        emnormall(4)  = 10**pos%logemline_oiii / 3.0 ! [OIII]
        emnormall(5)  = 10**pos%logemline_oiii       ! [OIII]
        emnormall(6)  = 10**pos%logemline_ni         ! [NI]
        emnormall(7)  = 10**pos%logemline_nii / 2.95 ! [NII]
        emnormall(8)  = 10**pos%logemline_h          ! Ha
        emnormall(9)  = 10**pos%logemline_nii        ! [NII]
        emnormall(10) = 10**pos%logemline_sii        ! [SII]
        emnormall(11) = 10**pos%logemline_sii * 0.7  ! [SII]

        DO i=1,neml
           !allow the em lines to be offset in velocity from the continuum
           !NB: velz2 is a *relative* shift between continuum and lines
           vz   = emlines(i) / (1+pos%velz2/clight*1E5)
           lsig = MAX(vz*pos%sigma2/clight*1E5,0.5)
           spec(1:nl_fit) = spec(1:nl_fit) + emnormall(i) * &
                EXP(-(sspgrid%lam(1:nl_fit)-vz)**2/lsig**2/2.0)
        ENDDO
     ENDIF

     !apply an atmosphereric transmission function
     IF (fit_trans.EQ.1) THEN
        !atm trans applied in the observed frame
        tmp_ltrans = sspgrid%lam / (1+pos%velz/clight*1E5)
        tmp_ftrans(1:nl_fit) = linterp(tmp_ltrans,&
             sspgrid%atm_trans,sspgrid%lam(1:nl_fit))
        spec(1:nl_fit) = spec(1:nl_fit) * &
             (1+(tmp_ftrans(1:nl_fit)-1)*10**pos%logtrans)
     ENDIF


  ENDIF

  !velocity broaden the model
  IF (pos%sigma.GT.5.0) &
       CALL VELBROAD(sspgrid%lam,spec,pos%sigma,l1(1),l2(nlint))

  !apply a template error function
  IF (apply_temperrfcn.EQ.1) THEN
     spec(1:nl_fit) = spec(1:nl_fit) / temperrfcn(1:nl_fit)
  ENDIF


END SUBROUTINE GETMODEL
