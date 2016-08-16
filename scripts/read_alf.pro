FUNCTION READ_ALF_ONE, file, nwalker=nwalker

  ;library correction factors from Schiavon 2007 Table 6
  ;note that I have forced the library corrections factors to
  ;be zero for [Fe/H]=0.0,0.2
  libfeh  = [-1.6,-1.4,-1.2,-1.0,-0.8,-0.6,-0.4,-0.2,0.0,0.2]
  libofe  = [0.6,0.5,0.5,0.4,0.3,0.2,0.2,0.1,0.0,0.0]
  libmgfe = [0.4,0.4,0.4,0.4,0.29,0.20,0.13,0.08,0.0,0.0]
  libcafe = [0.32,0.3,0.28,0.26,0.20,0.12,0.06,0.02,0.0,0.0]

  sdir = getenv('ALF_HOME')
  dir  = sdir+'/results/'
  openr,lun,dir+file,/get_lun
  ss = '#'
  WHILE strmid(ss,0,1) EQ '#' DO readf,lun,ss
  ts = strsplit(ss,' ',/extr)
  close,lun
  free_lun,lun

  ;what kind of file are we dealing with?
  errp   = strpos(file,'errp')
  sum    = strpos(file,'sum')
  simple = strpos(file,'simple')

  IF n_elements(ts) EQ 51 THEN BEGIN
     readcol,dir+file,chi2,velz,sigma,logage,zh,feh,afe,cfe,$
             nfe,nafe,mgfe,sife,kfe,cafe,tife,vfe,crfe,mnfe,cofe,nife,$
             cufe,srfe,bafe,eufe,teff,imf1,imf2,logfy,sigma2,$
             velz2,logm7g,hotteff,loghot,fy_logage,logtrans,d1,d2,d3,d4,d5,$
             jitter,imf3,logsky,imf4,imf5,m2lr,m2li,m2lk,m2lmwr,m2lmwi,m2lmwk,/sil
  ENDIF ELSE IF n_elements(ts) EQ 50 THEN BEGIN
     readcol,dir+file,chi2,velz,sigma,logage,zh,feh,afe,cfe,$
             nfe,nafe,mgfe,sife,kfe,cafe,tife,vfe,crfe,mnfe,cofe,nife,$
             cufe,srfe,bafe,eufe,teff,imf1,imf2,logfy,sigma2,$
             velz2,logm7g,hotteff,loghot,fy_logage,logtrans,d1,d2,d3,d4,d5,$
             jitter,imf3,logsky,imf4,m2lr,m2li,m2lk,m2lmwr,m2lmwi,m2lmwk,/sil
  ENDIF ELSE IF n_elements(ts) EQ 49 THEN BEGIN
     readcol,dir+file,chi2,velz,sigma,logage,zh,feh,afe,cfe,$
             nfe,nafe,mgfe,sife,kfe,cafe,tife,vfe,crfe,mnfe,cofe,nife,$
             cufe,srfe,bafe,eufe,teff,imf1,imf2,logfy,sigma2,$
             velz2,logm7g,hotteff,loghot,fy_logage,logtrans,d1,d2,d3,d4,d5,$
             jitter,imf3,logsky,m2lr,m2li,m2lk,m2lmwr,m2lmwi,m2lmwk,/sil
     imf4 = findgen(n_elements(chi2))
  ENDIF ELSE IF n_elements(ts) EQ 48 THEN BEGIN
     readcol,dir+file,chi2,velz,sigma,logage,zh,feh,afe,cfe,$
             nfe,nafe,mgfe,sife,kfe,cafe,tife,vfe,crfe,mnfe,cofe,nife,$
             cufe,srfe,bafe,eufe,teff,imf1,imf2,logfy,sigma2,$
             velz2,logm7g,hotteff,loghot,fy_logage,logtrans,d1,d2,d3,d4,d5,$
             jitter,imf3,m2lr,m2li,m2lk,m2lmwr,m2lmwi,m2lmwk,/sil
     imf4   = findgen(n_elements(chi2))
     logsky = findgen(n_elements(chi2))
  ENDIF ELSE IF n_elements(ts) EQ 47 THEN BEGIN
     readcol,dir+file,chi2,velz,sigma,logage,zh,feh,afe,cfe,$
             nfe,nafe,mgfe,sife,kfe,cafe,tife,vfe,crfe,mnfe,cofe,nife,$
             cufe,srfe,bafe,eufe,teff,imf1,imf2,logfy,sigma2,$
             velz2,logm7g,hotteff,loghot,fy_logage,logtrans,d1,d2,d3,d4,d5,$
             jitter,m2lr,m2li,m2lk,m2lmwr,m2lmwi,m2lmwk,/sil
     imf4   = findgen(n_elements(chi2))
     imf3   = findgen(n_elements(chi2))
     logsky = findgen(n_elements(chi2))
  ENDIF ELSE BEGIN
     print,'READ_ALF ERROR: file format not recognized, returning...'
     RETURN,0
  ENDELSE 

  str = {name:'',ngc:0L,logage:0.0,zh:0.0,feh:0.0,afe:0.0,cfe:0.0,nfe:0.0,$
         nafe:0.0,mgfe:0.0,sife:0.0,cafe:0.0,tife:0.0,crfe:0.0,mnfe:0.0,bafe:0.0,$
         srfe:0.0,nife:0.0,cufe:0.0,cofe:0.0,eufe:0.0,kfe:0.0,vfe:0.0,teff:0.0,$
         imf1:0.0,imf2:0.0,imf3:0.0,imf4:0.0,imf5:0.0,logfy:0.0,fy_logage:0.0,$
         sigma:0.0,sigma2:0.0,tfeh:0.0,velz:0.0,velz2:0.0,logm7g:0.0,hotteff:0.0,$
         loghot:0.0,logtrans:0.0,chi2:0.0,mlk:0.0,mli:0.0,mlr:0.0,mlk_mw:0.0,$
         mli_mw:0.0,mlr_mw:0.0,indgb:0.0,emline:fltarr(5),lsig:0.0,ml:0.0,$
         lage:99.0,vmag:99.,fuv:99.,nuv:99.,logemline_h:0.0,logemline_oiii:0.0,$
         logemline_sii:0.0,logemline_ni:0.0,logemline_nii:0.0,delafe:0.0,$
         delmgfe:0.0,delcafe:0.0,jitter:0.0,logsky:0.0,logm:0.0,imf:fltarr(5)}

  IF errp EQ -1 THEN BEGIN
     tfeh = feh
  ENDIF ELSE BEGIN
     tfeh = fltarr(n_elements(feh))
  ENDELSE

  res = replicate(str,n_elements(logage))

  ;these are the indices that are in the 'natural' units of the parameter
  ;as opposed to the 1 sigma error
  IF sum GT -1 THEN BEGIN
     sind = [0,1,3,4,5,6,7,8,9]
     eind = 2
  ENDIF ELSE BEGIN
     sind = indgen(n_elements(logage))
  ENDELSE


  res.chi2 = chi2
  res.zh   = zh
  res.tfeh = feh
  
  ;fold in the [Z/H] result into [Fe/H]
  IF errp EQ -1 THEN BEGIN
     res[sind].feh = res[sind].zh + feh[sind]
     IF sum GT -1 THEN $
        res[eind].feh = sqrt(res[eind].zh^2+feh[eind]^2)
  ENDIF ELSE BEGIN
     res.feh = sqrt(res.zh^2+feh^2)
  ENDELSE

  ;compute the library enhancement factors
  IF errp EQ -1 THEN BEGIN
     delafe  = fltarr(n_elements(feh))
     delmgfe = delafe
     delcafe = delafe
     FOR i=0,n_elements(feh)-1 DO BEGIN
        delafe[i]  = mylinterp(libfeh,libofe,res[i].zh)
        delmgfe[i] = mylinterp(libfeh,libmgfe,res[i].zh)
        delcafe[i] = mylinterp(libfeh,libcafe,res[i].zh)
     ENDFOR
  ENDIF ELSE BEGIN
     delafe  = 0.0
     delmgfe = 0.0
     delcafe = 0.0
  ENDELSE

  res[sind].delafe  = delafe[sind]
  res[sind].delcafe = delcafe[sind]
  res[sind].delmgfe = delmgfe[sind]

  res[sind].afe  = afe[sind] -tfeh[sind] + delafe[sind]
  res[sind].mgfe = mgfe[sind]-tfeh[sind] + delmgfe[sind]

  ;assume that Ca~Ti~Si
  res[sind].sife = sife[sind]-tfeh[sind] + delcafe[sind]
  res[sind].cafe = cafe[sind]-tfeh[sind] + delcafe[sind]
  res[sind].tife = tife[sind]-tfeh[sind] + delcafe[sind]

  ;these elements seem to show no net enhancement at low Z
  res[sind].cfe  = cfe[sind] -tfeh[sind]
  res[sind].nfe  = nfe[sind] -tfeh[sind]
  res[sind].crfe = crfe[sind]-tfeh[sind]
  res[sind].nife = nife[sind]-tfeh[sind]
  res[sind].nafe = nafe[sind]-tfeh[sind]

  ;these have enhancements that I have not yet quantified
  res[sind].bafe = bafe[sind]-tfeh[sind]
  res[sind].eufe = eufe[sind]-tfeh[sind]
  res[sind].srfe = srfe[sind]-tfeh[sind]
  res[sind].cufe = cufe[sind]-tfeh[sind]
  res[sind].cofe = cofe[sind]-tfeh[sind]
  res[sind].kfe  = kfe[sind] -tfeh[sind]
  res[sind].vfe  = vfe[sind] -tfeh[sind]
  res[sind].mnfe = mnfe[sind]-tfeh[sind]

  res.imf1   = imf1
  res.imf2   = imf2
  res.imf3   = imf3
  res.imf4   = imf4
  IF n_elements(ts) EQ 51 THEN BEGIN
     res.imf5 = imf5
  ENDIF ELSE BEGIN
     res.imf5   = 0.0
  ENDELSE
  FOR i=0,n_elements(logage)-1 DO $
     res[i].imf = [res[i].imf1,res[i].imf2,res[i].imf3,$
                   res[i].imf4,res[i].imf5]

  res.teff      = teff
  res.logage    = logage
  res.logfy     = logfy
  res.fy_logage = fy_logage
  res.sigma     = sigma
  res.velz      = velz
  res.sigma2    = sigma2
  res.velz2     = velz2
  res.logm7g    = logm7g
  res.hotteff   = hotteff
  res.loghot    = loghot
  res.logtrans  = logtrans
  res.jitter    = jitter
  res.logsky    = logsky

  res.mlr    = m2lr
  res.mli    = m2li
  res.mlk    = m2lk
  res.mlr_mw = m2lmwr
  res.mli_mw = m2lmwi
  res.mlk_mw = m2lmwk
  
  FOR i=0,n_elements(logage)-1 DO $
     res[i].emline = [d1[i],d2[i],d3[i],d4[i],d5[i]]
  res.logemline_h    = res.emline[0]
  res.logemline_oiii = res.emline[1]
  res.logemline_sii  = res.emline[2]
  res.logemline_ni   = res.emline[3]
  res.logemline_nii  = res.emline[4]

  IF errp EQ -1 THEN BEGIN
     res[sind].lsig = alog10(res[sind].sigma)
     res[sind].ml   = res[sind].mli/res[sind].mli_mw
     IF sum GT -1 THEN BEGIN
        res[eind].ml = sqrt(res[eind].mli^2+res[eind].mli_mw^2)/res[0].mli_mw
     ENDIF
  ENDIF

  ;mass-weighted age
  IF errp EQ -1 THEN BEGIN
     IF simple EQ -1 THEN BEGIN
        res[sind].lage = (1-10^res[sind].logfy)*10^res[sind].logage + $
                         10^res[sind].fy_logage*10^res[sind].logfy
     ENDIF ELSE BEGIN
        res[sind].lage = 10^res[sind].logage
     ENDELSE 
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

FUNCTION READ_ALF, file,nwalker=nwalker,bestp=bestp,errp=errp,$
                   cl2=cl2,cl16=cl16,cl50=cl50,cl84=cl84,cl97=cl97,$
                   minp=minp,prlo=prlo,prhi=prhi

  sdir = getenv('ALF_HOME')
  IF sdir EQ '' THEN BEGIN
     print,'READ_ALF ERROR: ALF_HOME environment '+$
           'variable not set, returning...'
     RETURN,0
  ENDIF

  IF n_elements(file) EQ 1 THEN BEGIN
     ff = findfile(sdir+'/results/'+file,count=ct)
     IF ct EQ 0 THEN print,'READ_ALF ERROR: no files found'
  ENDIF ELSE BEGIN
     ct = n_elements(file)
     ff = file
  ENDELSE

  FOR i=0,ct-1 DO BEGIN
     tf  = strsplit(ff[i],'/',/extr)
     tf  = tf[n_elements(tf)-1]
     one = read_alf_one(tf, nwalker=nwalker)
     none = n_elements(one)
     IF i EQ 0 THEN $
        all = replicate(one[0],n_elements(one),ct)
     all[*,i] = one
  ENDFOR

  IF none EQ 10 THEN BEGIN
     IF ct EQ 0 THEN bestp = all[0] ELSE bestp = reform(all[0,*])
     IF ct EQ 0 THEN minp  = all[1] ELSE minp  = reform(all[1,*])
     IF ct EQ 0 THEN errp  = all[2] ELSE errp  = reform(all[2,*])
     IF ct EQ 0 THEN cl2   = all[3] ELSE cl2   = reform(all[3,*])
     IF ct EQ 0 THEN cl16  = all[4] ELSE cl16  = reform(all[4,*])
     IF ct EQ 0 THEN cl50  = all[5] ELSE cl50  = reform(all[5,*])
     IF ct EQ 0 THEN cl84  = all[6] ELSE cl84  = reform(all[6,*])
     IF ct EQ 0 THEN cl97  = all[7] ELSE cl97  = reform(all[7,*])
     IF ct EQ 0 THEN prlo  = all[8] ELSE prlo  = reform(all[8,*])
     IF ct EQ 0 THEN prhi  = all[8] ELSE prhi  = reform(all[9,*])
  ENDIF
     
  RETURN,all

END

