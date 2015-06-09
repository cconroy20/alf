FUNCTION READ_ALF_ONE, file, old=old, nwalker=nwalker

  sdir = getenv('SPECFIT_HOME')
  dir  = sdir+'/results/'
  openr,lun,dir+file,/get_lun
  ss = ''
  readf,lun,ss
  close,lun
  free_lun,lun
  ts = strsplit(ss,' ',/extr)

  IF n_elements(ts) EQ 56 THEN BEGIN
     readcol,dir+file,chi2,velz,sigma,logage,feh,afe,nhe,cfe,$
             nfe,nafe,mgfe,sife,kfe,cafe,tife,vfe,crfe,mnfe,cofe,nife,$
             cufe,rbfe,srfe,yfe,zrfe,bafe,eufe,teff,imf1,imf2,logfy,sigma2,$
             velz2,logm7g,hotteff,loghot,fy_logage,d1,d2,d3,d4,d5,d6,d7,d8,$
             d9,d10,d11,d12,d13,m2lr,m2li,m2lk,m2lmwr,m2lmwi,m2lmwk,/sil
  ENDIF ELSE IF n_elements(ts) EQ 51 THEN BEGIN
     readcol,dir+file,chi2,velz,sigma,logage,feh,afe,nhe,cfe,$
             nfe,nafe,mgfe,sife,kfe,cafe,tife,vfe,crfe,mnfe,cofe,nife,$
             cufe,srfe,bafe,eufe,teff,imf1,imf2,logfy,sigma2,$
             velz2,logm7g,hotteff,loghot,fy_logage,d1,d2,d3,d4,d5,d6,d7,d8,$
             d9,d10,d11,m2lr,m2li,m2lk,m2lmwr,m2lmwi,m2lmwk,/sil
     yfe  = fltarr(n_elements(chi2))
     rbfe = fltarr(n_elements(chi2))
     zrfe = fltarr(n_elements(chi2))
     d12  = fltarr(n_elements(chi2))
     d13  = fltarr(n_elements(chi2))
  ENDIF ELSE BEGIN
     print,'READ_SPECFIT ERROR: file format not recognized, returning...'
     return,0
  ENDELSE 

  str = {name:'',ngc:0L,feh:0.0,afe:0.0,cfe:0.0,nfe:0.0,nafe:0.0,mgfe:0.0,$
         sife:0.0,cafe:0.0,tife:0.0,crfe:0.0,mnfe:0.0,bafe:0.0,$
         srfe:0.0,nife:0.0,cufe:0.0,cofe:0.0,eufe:0.0,kfe:0.0,vfe:0.0,$
         yfe:0.0,zrfe:0.0,rbfe:0.0,nhe:0.0,teff:0.0,imf1:0.0,logage:0.0,$
         logfy:0.0,fy_logage:0.0,sigma:0.0,sigma2:0.0,sigma3:0.0,sigma4:0.0,$
         sigma5:0.0,velz:0.0,velz2:0.0,logm7g:0.0,hotteff:0.0,loghot:0.0,$
         imf2:0.0,chi2:0.0,mlk:0.0,mli:0.0,mlr:0.0,mlk_mw:0.0,$
         mli_mw:0.0,mlr_mw:0.0,indgb:0.0,emline:fltarr(13),lsig:0.0,ml:0.0,$
         lage:99.0,vmag:99.,fuv:99.,nuv:99.}

  ss = strpos(file,'errp')
  IF ss EQ -1 THEN BEGIN
     tfeh = feh
  ENDIF ELSE BEGIN
     tfeh = fltarr(n_elements(feh))
  ENDELSE

  res = replicate(str,n_elements(logage))

  res.chi2 = chi2
  res.feh = feh
  res.afe = afe-tfeh
  res.cfe = cfe-tfeh
  res.nfe = nfe-tfeh
  res.nafe = nafe-tfeh
  res.mgfe = mgfe-tfeh
  res.sife = sife-tfeh
  res.cafe = cafe-tfeh
  res.tife = tife-tfeh
  res.crfe = crfe-tfeh
  res.mnfe = mnfe-tfeh
  res.bafe = bafe-tfeh
  res.eufe = eufe-tfeh
  res.srfe = srfe-tfeh
  res.nife = nife-tfeh
  res.cufe = cufe-tfeh
  res.cofe = cofe-tfeh
  res.kfe  = kfe-tfeh
  res.vfe  = vfe-tfeh
  res.nhe    = nhe
  res.teff   = teff
  res.imf1   = imf1
  res.imf2   = imf2
  res.logage = logage
  res.logfy  = logfy
  res.fy_logage = fy_logage
  res.sigma  = sigma
  res.velz   = velz
  res.sigma2 = sigma2
  res.velz2  = velz2
  res.logm7g = logm7g
  res.hotteff = hotteff
  res.loghot = loghot
  res.mlr    = m2lr
  res.mli    = m2li
  res.mlk    = m2lk
  res.mlr_mw = m2lmwr
  res.mli_mw = m2lmwi
  res.mlk_mw = m2lmwk
  FOR i=0,n_elements(logage)-1 DO $
     res[i].emline = [d1[i],d2[i],d3[i],d4[i],d5[i],d6[i],d7[i],d8[i],$
                      d9[i],d10[i],d11[i],d12[i],d13[i]]
  
  IF keyword_set(old) THEN BEGIN
     res.yfe = yfe-tfeh
     res.zrfe = zrfe-tfeh
     res.rbfe = rbfe-tfeh
  ENDIF

  res.lsig = alog10(res.sigma)
  res.ml   = res.mli/res.mli_mw

  IF strpos(file,'errp') EQ -1 THEN BEGIN
     ;light-weighted age - actually this is a mass-weighted age!
     res.lage = (1-10^res.logfy)*10^res.logage + 10^res.fy_logage*10^res.logfy
  ENDIF ELSE res.lage = -99.

  IF strpos(file,'simple') NE -1 THEN BEGIN
     ;in simple mode, we only have one component
     res.lage = 10^res.logage
  ENDIF

  ;split the chain into separate walkers
  IF keyword_set(nwalker) THEN BEGIN
     tres = replicate(res[0],nwalker,n_elements(res)/nwalker)
     FOR i=0,nwalker-1 DO BEGIN
        FOR j=0,n_elements(res)/nwalker-1 DO $
           tres[i,j] = res[j*nwalker+i]
     ENDFOR
     res = tres
  ENDIF

  RETURN,res

END

;----------------------------------------------------------------------;
;----------------------------------------------------------------------;

FUNCTION READ_ALF, file, old=old, nwalker=nwalker

  sdir = getenv('SPECFIT_HOME')
  IF sdir EQ '' THEN BEGIN
     print,'READ_ALF ERROR: SPECFIT_HOME environment '+$
           'variable not set, returning...'
     return,0
  ENDIF

  ff = findfile(sdir+'/results/'+file,count=ct)
  IF ct EQ 0 THEN print,'no files found'

  FOR i=0,ct-1 DO BEGIN
     tf  = strsplit(ff[i],'/',/extr)
     tf  = tf[n_elements(tf)-1]
     one = read_specfit_one(tf,old=old, nwalker=nwalker)
     all = (i EQ 0) ? one : [all,one]
  ENDFOR

  RETURN,all

END

