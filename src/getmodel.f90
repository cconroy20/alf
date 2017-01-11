SUBROUTINE GETMODEL(pos,spec,mw)

  !routine to produce a model spectrum (spec) for an input 
  !set of parameters (pos).  The optional flag 'mw' is used
  !to force the IMF to be of the MW (Kroupa 2001) form

  USE alf_vars; USE nr, ONLY : locate
  USE alf_utils, ONLY : velbroad, add_response,linterp,getmass
  IMPLICIT NONE

  TYPE(PARAMS), INTENT(in) :: pos
  REAL(DP), DIMENSION(nl), INTENT(out) :: spec
  INTEGER, OPTIONAL :: mw
  REAL(DP), DIMENSION(nl) :: tmp,tmpr,yspec,tmp1,tmp2,tmp3,tmp4
  INTEGER  :: vt,vy,vv1,vv2,vv3,i,vr,vm,vh,vm2,vm3
  REAL(DP) :: dt,fy,dx1,dx2,dx3,lsig,ve,dr,dm,dm2,dm3,dh,dy,tmps,inorm,mass,msto
  REAL(DP), DIMENSION(nl)   :: tmp_ltrans,tmp_ftrans_h2o,tmp_ftrans_o2
  REAL(DP), DIMENSION(neml) :: emnormall=1.0
  REAL(DP), DIMENSION(nimfnp) :: imfw=0.0
  REAL(DP), DIMENSION(2) :: hermite=0.0

  !---------------------------------------------------------------!
  !---------------------------------------------------------------!

  !set up interpolants for age
  vt = MAX(MIN(locate(sspgrid%logagegrid,pos%logage),nage-1),1)
  dt = (pos%logage-sspgrid%logagegrid(vt))/&
       (sspgrid%logagegrid(vt+1)-sspgrid%logagegrid(vt))
  dt = MAX(MIN(dt,1.2),-0.3)  !0.5<age<14 Gyr

  !set up interpolants for metallicity
  vm = MAX(MIN(locate(sspgrid%logzgrid,pos%zh),nzmet-1),1)
  dm = (pos%zh-sspgrid%logzgrid(vm)) / &
       (sspgrid%logzgrid(vm+1)-sspgrid%logzgrid(vm))
  dm = MAX(MIN(dm,1.0),-1.0) ! -2.0<[Z/H]<0.25

  !compute the IMF-variable SSP
  IF (mwimf.EQ.0.AND..NOT.PRESENT(mw).AND.&
       fit_type.EQ.0.AND.powell_fitting.EQ.0) THEN
     
     vv1 = MAX(MIN(locate(sspgrid%imfx1,pos%imf1),nimf-1),1)
     dx1 = (pos%imf1-sspgrid%imfx1(vv1))/&
          (sspgrid%imfx1(vv1+1)-sspgrid%imfx1(vv1))
     dx1 = MAX(MIN(dx1,1.0),0.0)

     IF (imf_type.EQ.0.OR.imf_type.EQ.2) THEN
        !single power-law slope for IMF=0,2
        vv2 = vv1
        dx2 = dx1
     ELSE
        !two-part power-law for IMF=1,3
        vv2 = MAX(MIN(locate(sspgrid%imfx2,pos%imf2),nimf-1),1)
        dx2 = (pos%imf2-sspgrid%imfx2(vv2))/&
             (sspgrid%imfx2(vv2+1)-sspgrid%imfx2(vv2))
        dx2 = MAX(MIN(dx2,1.0),0.0)        
     ENDIF
     IF (imf_type.EQ.2.OR.imf_type.EQ.3) THEN
        !variable low-mass cutoff for IMF=2,3
        vv3 = MAX(MIN(locate(sspgrid%imfx3,pos%imf3),nmcut-1),1)
        dx3 = (pos%imf3-sspgrid%imfx3(vv3))/&
             (sspgrid%imfx3(vv3+1)-sspgrid%imfx3(vv3))
        dx3 = MAX(MIN(dx3,1.0),0.0)
     ENDIF

     IF (imf_type.EQ.2.OR.imf_type.EQ.3) THEN

        vm3 = MAX(MIN(locate(sspgrid%logzgrid2,pos%zh),nzmet3-1),1)
        dm3 = (pos%zh-sspgrid%logzgrid2(vm3)) / &
             (sspgrid%logzgrid2(vm3+1)-sspgrid%logzgrid2(vm3))
        dm3 = MAX(MIN(dm3,1.5),-1.0) 

        tmp1 = (1-dx1)*(1-dx2)*(1-dx3)*sspgrid%logsspm(:,vv1,vv2,vt+1,vv3,vm3+1) + &
               (dx1)*(1-dx2)*(1-dx3)*sspgrid%logsspm(:,vv1+1,vv2,vt+1,vv3,vm3+1) + &
               (1-dx1)*(dx2)*(1-dx3)*sspgrid%logsspm(:,vv1,vv2+1,vt+1,vv3,vm3+1) + &
               (1-dx1)*(1-dx2)*(dx3)*sspgrid%logsspm(:,vv1,vv2,vt+1,vv3+1,vm3+1) + &
               (dx1)*(1-dx2)*(dx3)*sspgrid%logsspm(:,vv1+1,vv2,vt+1,vv3+1,vm3+1) + &
               (1-dx1)*(dx2)*(dx3)*sspgrid%logsspm(:,vv1,vv2+1,vt+1,vv3+1,vm3+1) + &
               (dx1)*(dx2)*(1-dx3)*sspgrid%logsspm(:,vv1+1,vv2+1,vt+1,vv3,vm3+1) + &
                     dx1*dx2*dx3*sspgrid%logsspm(:,vv1+1,vv2+1,vt+1,vv3+1,vm3+1) 

        tmp2 = (1-dx1)*(1-dx2)*(1-dx3)*sspgrid%logsspm(:,vv1,vv2,vt,vv3,vm3+1) + &
               (dx1)*(1-dx2)*(1-dx3)*sspgrid%logsspm(:,vv1+1,vv2,vt,vv3,vm3+1) + &
               (1-dx1)*(dx2)*(1-dx3)*sspgrid%logsspm(:,vv1,vv2+1,vt,vv3,vm3+1) + &
               (1-dx1)*(1-dx2)*(dx3)*sspgrid%logsspm(:,vv1,vv2,vt,vv3+1,vm3+1) + &
               (dx1)*(1-dx2)*(dx3)*sspgrid%logsspm(:,vv1+1,vv2,vt,vv3+1,vm3+1) + &
               (1-dx1)*(dx2)*(dx3)*sspgrid%logsspm(:,vv1,vv2+1,vt,vv3+1,vm3+1) + &
               (dx1)*(dx2)*(1-dx3)*sspgrid%logsspm(:,vv1+1,vv2+1,vt,vv3,vm3+1) + &
                     dx1*dx2*dx3*sspgrid%logsspm(:,vv1+1,vv2+1,vt,vv3+1,vm3+1) 
        
        tmp3 = (1-dx1)*(1-dx2)*(1-dx3)*sspgrid%logsspm(:,vv1,vv2,vt+1,vv3,vm3) + &
               (dx1)*(1-dx2)*(1-dx3)*sspgrid%logsspm(:,vv1+1,vv2,vt+1,vv3,vm3) + &
               (1-dx1)*(dx2)*(1-dx3)*sspgrid%logsspm(:,vv1,vv2+1,vt+1,vv3,vm3) + &
               (1-dx1)*(1-dx2)*(dx3)*sspgrid%logsspm(:,vv1,vv2,vt+1,vv3+1,vm3) + &
               (dx1)*(1-dx2)*(dx3)*sspgrid%logsspm(:,vv1+1,vv2,vt+1,vv3+1,vm3) + &
               (1-dx1)*(dx2)*(dx3)*sspgrid%logsspm(:,vv1,vv2+1,vt+1,vv3+1,vm3) + &
               (dx1)*(dx2)*(1-dx3)*sspgrid%logsspm(:,vv1+1,vv2+1,vt+1,vv3,vm3) + &
                     dx1*dx2*dx3*sspgrid%logsspm(:,vv1+1,vv2+1,vt+1,vv3+1,vm3) 

        tmp4 = (1-dx1)*(1-dx2)*(1-dx3)*sspgrid%logsspm(:,vv1,vv2,vt,vv3,vm3) + &
               (dx1)*(1-dx2)*(1-dx3)*sspgrid%logsspm(:,vv1+1,vv2,vt,vv3,vm3) + &
               (1-dx1)*(dx2)*(1-dx3)*sspgrid%logsspm(:,vv1,vv2+1,vt,vv3,vm3) + &
               (1-dx1)*(1-dx2)*(dx3)*sspgrid%logsspm(:,vv1,vv2,vt,vv3+1,vm3) + &
               (dx1)*(1-dx2)*(dx3)*sspgrid%logsspm(:,vv1+1,vv2,vt,vv3+1,vm3) + &
               (1-dx1)*(dx2)*(dx3)*sspgrid%logsspm(:,vv1,vv2+1,vt,vv3+1,vm3) + &
               (dx1)*(dx2)*(1-dx3)*sspgrid%logsspm(:,vv1+1,vv2+1,vt,vv3,vm3) + &
                     dx1*dx2*dx3*sspgrid%logsspm(:,vv1+1,vv2+1,vt,vv3+1,vm3) 
  
        spec = 10**( dt*dm3*tmp1 + (1-dt)*dm3*tmp2 + &
             dt*(1-dm3)*tmp3 + (1-dt)*(1-dm3)*tmp4 )

     ELSE IF (imf_type.EQ.0.OR.imf_type.EQ.1) THEN

        tmp1 = (1-dx1)*(1-dx2)*sspgrid%logssp(:,vv1,vv2,vt+1,vm+1)+&
             dx1*(1-dx2)*sspgrid%logssp(:,vv1+1,vv2,vt+1,vm+1)+&
             (1-dx1)*dx2*sspgrid%logssp(:,vv1,vv2+1,vt+1,vm+1)+&
             dx1*dx2*sspgrid%logssp(:,vv1+1,vv2+1,vt+1,vm+1)
      
        tmp2 = (1-dx1)*(1-dx2)*sspgrid%logssp(:,vv1,vv2,vt,vm+1)+&
             dx1*(1-dx2)*sspgrid%logssp(:,vv1+1,vv2,vt,vm+1)+&
             (1-dx1)*dx2*sspgrid%logssp(:,vv1,vv2+1,vt,vm+1)+&
             dx1*dx2*sspgrid%logssp(:,vv1+1,vv2+1,vt,vm+1)
 
        tmp3 = (1-dx1)*(1-dx2)*sspgrid%logssp(:,vv1,vv2,vt+1,vm)+&
             dx1*(1-dx2)*sspgrid%logssp(:,vv1+1,vv2,vt+1,vm)+&
             (1-dx1)*dx2*sspgrid%logssp(:,vv1,vv2+1,vt+1,vm)+&
             dx1*dx2*sspgrid%logssp(:,vv1+1,vv2+1,vt+1,vm)
      
        tmp4 = (1-dx1)*(1-dx2)*sspgrid%logssp(:,vv1,vv2,vt,vm)+&
             dx1*(1-dx2)*sspgrid%logssp(:,vv1+1,vv2,vt,vm)+&
             (1-dx1)*dx2*sspgrid%logssp(:,vv1,vv2+1,vt,vm)+&
             dx1*dx2*sspgrid%logssp(:,vv1+1,vv2+1,vt,vm)
     
        spec = 10**( dt*dm*tmp1 + (1-dt)*dm*tmp2 + &
             dt*(1-dm)*tmp3 + (1-dt)*(1-dm)*tmp4 )

     ELSE IF (imf_type.EQ.4) THEN
        
        !non-parametric IMF

        imfw(1) = 10**pos%imf1
        imfw(2) = 10**((pos%imf2+pos%imf1)/2.)
        imfw(3) = 10**pos%imf2
        imfw(4) = 10**((pos%imf3+pos%imf2)/2.)
        imfw(5) = 10**pos%imf3
        imfw(6) = 10**((pos%imf4+pos%imf3)/2.)
        imfw(7) = 10**pos%imf4
        imfw(8) = 10**((imf5+pos%imf4)/2.)
        imfw(9) = 10**imf5

        tmp1 = 0.0
        tmp2 = 0.0
        tmp3 = 0.0
        tmp4 = 0.0
        DO i=1,nimfnp
           tmp1 = tmp1 + imfw(i)*sspgrid%sspnp(:,i,vt+1,vm+1)
           tmp2 = tmp2 + imfw(i)*sspgrid%sspnp(:,i,vt,vm+1)
           tmp3 = tmp3 + imfw(i)*sspgrid%sspnp(:,i,vt+1,vm)
           tmp4 = tmp4 + imfw(i)*sspgrid%sspnp(:,i,vt,vm)
        ENDDO

        msto = MAX(MIN(10**(msto_t0+msto_t1*sspgrid%logagegrid(vt+1)) * &
             (msto_z0+msto_z1*sspgrid%logzgrid(vm+1)+&
             msto_z2*sspgrid%logzgrid(vm+1)**2),3.0),0.75)
        mass = getmass(imflo,msto,pos%imf1,pos%imf2,krpa_imf3,&
             pos%imf3,pos%imf4,inorm)
        tmp1 = tmp1/inorm

        msto = MAX(MIN(10**(msto_t0+msto_t1*sspgrid%logagegrid(vt)) * &
             (msto_z0+msto_z1*sspgrid%logzgrid(vm+1)+&
             msto_z2*sspgrid%logzgrid(vm+1)**2),3.0),0.75)
        mass = getmass(imflo,msto,pos%imf1,pos%imf2,krpa_imf3,&
             pos%imf3,pos%imf4,inorm)
        tmp2 = tmp2/inorm

        msto = MAX(MIN(10**(msto_t0+msto_t1*sspgrid%logagegrid(vt+1)) * &
             (msto_z0+msto_z1*sspgrid%logzgrid(vm)+&
             msto_z2*sspgrid%logzgrid(vm)**2),3.0),0.75)
        mass = getmass(imflo,msto,pos%imf1,pos%imf2,krpa_imf3,&
             pos%imf3,pos%imf4,inorm)
        tmp3 = tmp3/inorm

        msto = MAX(MIN(10**(msto_t0+msto_t1*sspgrid%logagegrid(vt)) * &
             (msto_z0+msto_z1*sspgrid%logzgrid(vm)+&
             msto_z2*sspgrid%logzgrid(vm)**2),3.0),0.75)
        mass = getmass(imflo,msto,pos%imf1,pos%imf2,krpa_imf3,&
             pos%imf3,pos%imf4,inorm)
        tmp4 = tmp4/inorm

        spec = 10**( dt*dm*LOG10(tmp1) + (1-dt)*dm*LOG10(tmp2) + &
             dt*(1-dm)*LOG10(tmp3) + (1-dt)*(1-dm)*LOG10(tmp4) )

     ENDIF

  ELSE
  
     !compute a Kroupa IMF
     spec = 10**( &
          dt*dm*sspgrid%logssp(:,imfr1,imfr2,vt+1,vm+1) + &
          (1-dt)*dm*sspgrid%logssp(:,imfr1,imfr2,vt,vm+1) + & 
          dt*(1-dm)*sspgrid%logssp(:,imfr1,imfr2,vt+1,vm) + & 
          (1-dt)*(1-dm)*sspgrid%logssp(:,imfr1,imfr2,vt,vm) )

  ENDIF
  
  !vary young population - both fraction and age
  !only include these parameters in the "full" model
  IF (fit_type.EQ.0.AND.powell_fitting.EQ.0.AND.fit_two_ages.EQ.1) THEN
     fy    = MAX(MIN(10**pos%logfy,1.0),0.0)
     vy    = MAX(MIN(locate(sspgrid%logagegrid,pos%fy_logage),nage-1),1)
     dy    = (pos%fy_logage-sspgrid%logagegrid(vy))/&
          (sspgrid%logagegrid(vy+1)-sspgrid%logagegrid(vy))
     dy    = MAX(MIN(dy,1.0),-0.3) !0.5<age<13.5 Gyr
     yspec = &
          dy*dm*sspgrid%logssp(:,imfr1,imfr2,vy+1,vm+1) + &
          (1-dy)*dm*sspgrid%logssp(:,imfr1,imfr2,vy,vm+1) + & 
          dy*(1-dm)*sspgrid%logssp(:,imfr1,imfr2,vy+1,vm) + & 
          (1-dy)*(1-dm)*sspgrid%logssp(:,imfr1,imfr2,vy,vm)
      spec  = (1-fy)*spec + fy*10**yspec
  ENDIF

  !vary age in the response functions
  IF (use_age_dep_resp_fcns.EQ.0) THEN
     !force the use of the response fcn at age=fix_age_dep_resp_fcns
     vr = MAX(MIN(locate(sspgrid%logagegrid_rfcn,LOG10(fix_age_dep_resp_fcns)),&
          nage_rfcn-1),1)
     dr = (LOG10(fix_age_dep_resp_fcns)-sspgrid%logagegrid_rfcn(vr))/&
          (sspgrid%logagegrid_rfcn(vr+1)-sspgrid%logagegrid_rfcn(vr))
     dr = MAX(MIN(dr,1.0),0.0)
   ELSE
     !should be using mass-weighted age here
     vr = MAX(MIN(locate(sspgrid%logagegrid_rfcn,pos%logage),nage_rfcn-1),1)
     dr = (pos%logage-sspgrid%logagegrid_rfcn(vr))/&
          (sspgrid%logagegrid_rfcn(vr+1)-sspgrid%logagegrid_rfcn(vr))
     dr = MAX(MIN(dr,1.0),0.0)
  ENDIF

  !vary metallicity in the response functions
  IF (use_z_dep_resp_fcns.EQ.0) THEN
     vm2   = MAX(MIN(locate(sspgrid%logzgrid,fix_z_dep_resp_fcns),nzmet-1),1)
     dm2   = (fix_z_dep_resp_fcns-sspgrid%logzgrid(vm2)) / &
          (sspgrid%logzgrid(vm2+1)-sspgrid%logzgrid(vm2))
     dm2 = MAX(MIN(dm2,1.0),0.0)
  ELSE
     vm2 = vm
     dm2 = dm
  ENDIF

  !Only sigma, velz, logage, and [Z/H] are fit when either
  !fitting in Powell mode or "super simple" mode
  IF (powell_fitting.EQ.0.AND.fit_type.NE.2) THEN

     !vary [Fe/H]
     CALL ADD_RESPONSE(spec,pos%feh,0.3,dr,vr,dm2,vm2,sspgrid%solar,&
          sspgrid%fep,sspgrid%fem)
     !vary [O/H]
     CALL ADD_RESPONSE(spec,pos%ah,0.3,dr,vr,dm2,vm2,sspgrid%solar,sspgrid%ap)
     !vary [C/H]
     CALL ADD_RESPONSE(spec,pos%ch,0.15,dr,vr,dm2,vm2,sspgrid%solar,&
          sspgrid%cp,sspgrid%cm)
     !vary [N/H]
     CALL ADD_RESPONSE(spec,pos%nh,0.3,dr,vr,dm2,vm2,sspgrid%solar,&
          sspgrid%np,sspgrid%nm)
     !vary [Mg/H]
     CALL ADD_RESPONSE(spec,pos%mgh,0.3,dr,vr,dm2,vm2,sspgrid%solar,&
          sspgrid%mgp,sspgrid%mgm)
     !vary [Si/H]
     CALL ADD_RESPONSE(spec,pos%sih,0.3,dr,vr,dm2,vm2,sspgrid%solar,&
          sspgrid%sip,sspgrid%sim)
     !vary [Ca/H]
     CALL ADD_RESPONSE(spec,pos%cah,0.3,dr,vr,dm2,vm2,sspgrid%solar,&
          sspgrid%cap,sspgrid%cam)
     !vary [Ti/H]
     CALL ADD_RESPONSE(spec,pos%tih,0.3,dr,vr,dm2,vm2,sspgrid%solar,&
          sspgrid%tip,sspgrid%tim)
     
     !vary [Na/H] (special case)
     IF (pos%nah.LT.0.3) THEN

        CALL ADD_RESPONSE(spec,pos%nah,0.3,dr,vr,dm2,vm2,sspgrid%solar,&
             sspgrid%nap,sspgrid%nam)

     ELSE IF (pos%nah.GE.0.3.AND.pos%nah.LT.0.6) THEN

        tmpr = &
             dr*dm2*sspgrid%nap(:,vr+1,vm2+1)/sspgrid%solar(:,vr+1,vm2+1) + &
             (1-dr)*dm2*sspgrid%nap(:,vr,vm2+1)/sspgrid%solar(:,vr,vm2+1) + &
             dr*(1-dm2)*sspgrid%nap(:,vr+1,vm2)/sspgrid%solar(:,vr+1,vm2) + &
             (1-dr)*(1-dm2)*sspgrid%nap(:,vr,vm2)/sspgrid%solar(:,vr,vm2)
        tmp = &
             dr*dm2*(sspgrid%nap6(:,vr+1,vm2+1)-sspgrid%nap(:,vr+1,vm2+1))/&
             sspgrid%solar(:,vr+1,vm2+1)+&
             (1-dr)*dm2*(sspgrid%nap6(:,vr,vm2+1)-sspgrid%nap(:,vr,vm2+1))/&
             sspgrid%solar(:,vr,vm2+1)+&
             dr*(1-dm2)*(sspgrid%nap6(:,vr+1,vm2)-sspgrid%nap(:,vr+1,vm2))/&
             sspgrid%solar(:,vr+1,vm2)+&
             (1-dr)*(1-dm2)*(sspgrid%nap6(:,vr,vm2)-sspgrid%nap(:,vr,vm2))/&
             sspgrid%solar(:,vr,vm2)
 
       spec = spec * (tmpr+tmp*(pos%nah-0.3)/0.3 )

     ELSE IF (pos%nah.GE.0.6) THEN

        tmpr = &
             dr*dm2*sspgrid%nap6(:,vr+1,vm2+1)/sspgrid%solar(:,vr+1,vm2+1) + &
             (1-dr)*dm2*sspgrid%nap6(:,vr,vm2+1)/sspgrid%solar(:,vr,vm2+1) + &
             dr*(1-dm2)*sspgrid%nap6(:,vr+1,vm2)/sspgrid%solar(:,vr+1,vm2) + &
             (1-dr)*(1-dm2)*sspgrid%nap6(:,vr,vm2)/sspgrid%solar(:,vr,vm2)
        
        tmp = &
             dr*dm2*(sspgrid%nap9(:,vr+1,vm2+1)-sspgrid%nap6(:,vr+1,vm2+1))/&
             sspgrid%solar(:,vr+1,vm2+1)+&
             (1-dr)*dm2*(sspgrid%nap9(:,vr,vm2+1)-sspgrid%nap6(:,vr,vm2+1))/&
             sspgrid%solar(:,vr,vm2+1)+&
             dr*(1-dm2)*(sspgrid%nap9(:,vr+1,vm2)-sspgrid%nap6(:,vr+1,vm2))/&
             sspgrid%solar(:,vr+1,vm2)+&
             (1-dr)*(1-dm2)*(sspgrid%nap9(:,vr,vm2)-sspgrid%nap6(:,vr,vm2))/&
             sspgrid%solar(:,vr,vm2)

        spec = spec * (tmpr+tmp*(pos%nah-0.6)/0.6 )

     ENDIF
     
  ENDIF

  !only include these parameters in the "full" model
  IF (fit_type.EQ.0.AND.powell_fitting.EQ.0) THEN

     !vary Teff (special case - force use of the 13 Gyr model)
     CALL ADD_RESPONSE(spec,pos%teff,50.,1.d0,nage_rfcn-1,dm2,vm2,sspgrid%solar,&
          sspgrid%teffp,sspgrid%teffm)

     !add a hot star
     vh   = MAX(MIN(locate(sspgrid%teffarrhot,pos%hotteff),nhot-1),1)
     dh   = (pos%hotteff-sspgrid%teffarrhot(vh))/&
          (sspgrid%teffarrhot(vh+1)-sspgrid%teffarrhot(vh))
     fy   = MAX(MIN(10**pos%loghot,1.0),0.0)
     tmp  = dh*sspgrid%hotspec(:,vh+1) + (1-dh)*sspgrid%hotspec(:,vh)
     spec = (1-fy)*spec + fy*tmp
     
     !add in an M7 giant
     fy   = MAX(MIN(10**pos%logm7g,1.0),0.0)
     spec = (1-fy)*spec + fy*sspgrid%m7g

     !vary [K/H]
     CALL ADD_RESPONSE(spec,pos%kh,0.3,dr,vr,dm2,vm2,sspgrid%solar,sspgrid%kp)
     !vary [V/H]
     CALL ADD_RESPONSE(spec,pos%vh,0.3,dr,vr,dm2,vm2,sspgrid%solar,sspgrid%vp)
     !vary [Cr/H]
     CALL ADD_RESPONSE(spec,pos%crh,0.3,dr,vr,dm2,vm2,sspgrid%solar,sspgrid%crp)
     !vary [Mn/H]
     CALL ADD_RESPONSE(spec,pos%mnh,0.3,dr,vr,dm2,vm2,sspgrid%solar,sspgrid%mnp)
     !vary [Co/H]
     CALL ADD_RESPONSE(spec,pos%coh,0.3,dr,vr,dm2,vm2,sspgrid%solar,sspgrid%cop)
     !vary [Ni/H]
     CALL ADD_RESPONSE(spec,pos%nih,0.3,dr,vr,dm2,vm2,sspgrid%solar,sspgrid%nip)
     !vary [Cu/H]
     CALL ADD_RESPONSE(spec,pos%cuh,0.3,dr,vr,dm2,vm2,sspgrid%solar,sspgrid%cup)
     !vary [Sr/H]
     CALL ADD_RESPONSE(spec,pos%srh,0.3,dr,vr,dm2,vm2,sspgrid%solar,sspgrid%srp)
     !vary [Ba/H]
     CALL ADD_RESPONSE(spec,pos%bah,0.3,dr,vr,dm2,vm2,sspgrid%solar,&
          sspgrid%bap,sspgrid%bam)
     !vary [Eu/H]
     CALL ADD_RESPONSE(spec,pos%euh,0.3,dr,vr,dm2,vm2,sspgrid%solar,sspgrid%eup)

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
           ve   = emlines(i) / (1+pos%velz2/clight*1E5)
           lsig = MAX(ve*pos%sigma2/clight*1E5,1.0) !min dlam=1.0A
           !maybe this should be spec*(1+emission lines)
           !otherwise the em lines change when the continuum changes...
           spec = spec + emnormall(i) * EXP(-(sspgrid%lam-ve)**2/lsig**2/2.0)
        ENDDO
     ENDIF

  ENDIF

  !velocity broaden the model
  IF (pos%sigma.GT.5.0.AND.fit_indices.EQ.0) THEN

     IF (fit_hermite.EQ.1) THEN
        hermite(1) = pos%h3
        hermite(2) = pos%h4
        CALL VELBROAD(sspgrid%lam,spec,pos%sigma,l1(1),l2(nlint),hermite)
     ELSE
        CALL VELBROAD(sspgrid%lam,spec,pos%sigma,l1(1),l2(nlint))
     ENDIF
     
  ENDIF

  !apply an atmospheric transmission function only in full mode
  !note that this is done *after* velocity broadening
  IF (fit_type.EQ.0.AND.powell_fitting.EQ.0.AND.fit_trans.EQ.1) THEN
     !applied in the observed frame
     tmp_ltrans     = sspgrid%lam / (1+pos%velz/clight*1E5)
     tmp_ftrans_h2o = linterp(tmp_ltrans,sspgrid%atm_trans_h2o,sspgrid%lam)
     tmp_ftrans_o2  = linterp(tmp_ltrans,sspgrid%atm_trans_o2,sspgrid%lam)
     spec = spec * (1+(tmp_ftrans_h2o-1)*10**pos%logtrans)
     spec = spec * (1+(tmp_ftrans_o2-1)*10**pos%logtrans)
  ENDIF

  !apply a template error function
  IF (apply_temperrfcn.EQ.1) THEN
     spec = spec / temperrfcn
  ENDIF


END SUBROUTINE GETMODEL
