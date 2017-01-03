FUNCTION READ_SPECBIN, file, lam=lam

  sdir = getenv('ALF_HOME')

  readcol,sdir+'/results/'+file+'.sum',s1,name,d2,val,/sil,format='(A,A,A,I)'

  wh = where(name EQ 'Nwalkers',ct1)
  nw = long(val[wh])
  wh = where(name EQ 'Nchain',ct2)
  nc = long(val[wh])
  wh = where(name EQ 'Nsample',ct3)
  ns = long(val[wh])
  wh = where(name EQ 'Nwave',ct4)
  nl = long(val[wh])

  IF ct1 EQ 0 OR ct2 EQ 0 OR ct3 EQ 0 OR ct4 EQ 0 THEN BEGIN
     print,'ERROR: missing one or more header entries in *sum file'
     RETURN,0
  ENDIF

  arr = read_binary(sdir+'/results/'+file+'.spec',data_type=4,$
                    data_dims=[nl,nw*nc/ns+1])

  lam = reform(arr[*,0])
  arr = arr[*,1:*]

  RETURN, arr

END
