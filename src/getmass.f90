FUNCTION GETMASS(mlo,mto,imf1,imf2,imfup,imf3,imf4,timfnorm)

  !compute mass in stars and remnants (normalized to 1 Msun at t=0)
  !assume an IMF that runs from 0.08 to 100 Msun.

  USE alf_vars
  IMPLICIT NONE

  !turnoff mass
  REAL(DP), INTENT(in) :: mlo,mto,imf1,imf2,imfup
  REAL(DP), INTENT(in), OPTIONAL :: imf3,imf4
  REAL(DP), INTENT(inout), OPTIONAL :: timfnorm
  INTEGER :: i
  REAL(DP) :: imfnorm, getmass
  REAL(DP), PARAMETER :: bhlim=40.0,nslim=8.5
  REAL(DP) :: m2=0.5,m3=1.0,alpha
  REAL(DP), DIMENSION(nimfnp) :: imfw=0.0, alpha2

  !---------------------------------------------------------------!
  !---------------------------------------------------------------!

  IF (mlo.GT.m2) THEN
     WRITE(*,*) 'GETMASS ERROR: mlo>m2'
     STOP
  ENDIF

  getmass  = 0.d0

  IF (.NOT.PRESENT(imf4)) THEN

     !normalize the weights so that 1 Msun formed at t=0
     imfnorm = (m2**(-imf1+2)-mlo**(-imf1+2))/(-imf1+2) + &
          m2**(-imf1+imf2)*(m3**(-imf2+2)-m2**(-imf2+2))/(-imf2+2) + &
          m2**(-imf1+imf2)*(imfhi**(-imfup+2)-m3**(-imfup+2))/(-imfup+2)

     !stars still alive
     getmass = (m2**(-imf1+2)-mlo**(-imf1+2))/(-imf1+2)
     IF (mto.LT.m3) THEN
        getmass = getmass + m2**(-imf1+imf2)*(mto**(-imf2+2)-m2**(-imf2+2))/(-imf2+2)
     ELSE
        getmass = getmass + m2**(-imf1+imf2)*(m3**(-imf2+2)-m2**(-imf2+2))/(-imf2+2) + &
             m2**(-imf1+imf2)*(mto**(-imfup+2)-m3**(-imfup+2))/(-imfup+2)
     ENDIF
     getmass = getmass/imfnorm

     !BH remnants
     !40<M<imf_up leave behind a 0.5*M BH
     getmass = getmass + &
          0.5*m2**(-imf1+imf2)*(imfhi**(-imfup+2)-bhlim**(-imfup+2))/(-imfup+2)/imfnorm

     !NS remnants
     !8.5<M<40 leave behind 1.4 Msun NS
     getmass = getmass + &
          1.4*m2**(-imf1+imf2)*(bhlim**(-imfup+1)-nslim**(-imfup+1))/(-imfup+1)/imfnorm

     !WD remnants
     !M<8.5 leave behind 0.077*M+0.48 WD
     IF (mto.LT.m3) THEN
        getmass = getmass + &
             0.48*m2**(-imf1+imf2)*(nslim**(-imfup+1)-m3**(-imfup+1))/(-imfup+1)/imfnorm
        getmass = getmass + &
             0.48*m2**(-imf1+imf2)*(m3**(-imf2+1)-mto**(-imf2+1))/(-imf2+1)/imfnorm
        getmass = getmass + &
             0.077*m2**(-imf1+imf2)*(nslim**(-imfup+2)-m3**(-imfup+2))/(-imfup+2)/imfnorm
        getmass = getmass + &
             0.077*m2**(-imf1+imf2)*(m3**(-imf2+2)-mto**(-imf2+2))/(-imf2+2)/imfnorm
     ELSE
        getmass = getmass + &
             0.48*m2**(-imf1+imf2)*(nslim**(-imfup+1)-mto**(-imfup+1))/(-imfup+1)/imfnorm
        getmass = getmass + &
             0.077*m2**(-imf1+imf2)*(nslim**(-imfup+2)-mto**(-imfup+2))/(-imfup+2)/imfnorm
     ENDIF

  ELSE

     !non-parametric IMF
     
     !mbin_nimf = (/0.2,0.4,0.6,0.8,1.0/)

     alpha2 = 2.0-npi_alphav

     imfw(1) = 10**imf1 
     imfw(2) = 10**((imf2+imf1)/2.)
     imfw(3) = 10**imf2
     imfw(4) = 10**((imf3+imf2)/2.)
     imfw(5) = 10**imf3
     imfw(6) = 10**((imf4+imf3)/2.)
     imfw(7) = 10**imf4
     imfw(8) = 10**((imf5+imf4)/2.)
     imfw(9) = 10**imf5

     DO i=1,nimfnp
        imfw(i) = imfw(i) * npi_renorm(i)
     ENDDO
     
     imfnorm = 0.0
     DO i=1,nimfnp
        imfnorm = imfnorm + imfw(i)/alpha2(i) * &
             (mbin_nimf9(i+1)**alpha2(i)-mbin_nimf9(i)**alpha2(i))
     ENDDO
     imfnorm = imfnorm + imfw(9)/(-imfup+2)/(mbin_nimf9(10)**(-imfup)) * &
          (imfhi**(-imfup+2)-mbin_nimf9(10)**(-imfup+2)) 

     IF (mto.GT.mbin_nimf9(10)) THEN

        ! MSTO > 1.0
        
        DO i=1,nimfnp
           getmass = getmass + imfw(i)/alpha2(i)*(mbin_nimf9(i+1)**alpha2(i) - &
                mbin_nimf9(i)**alpha2(i))
        ENDDO
        getmass = getmass + imfw(9)/(-imfup+2)/(mbin_nimf9(10)**(-imfup)) * &
             (mto**(-imfup+2)-mbin_nimf9(10)**(-imfup+2)) 
        
        !WD remnants from MTO-8.5
        getmass = getmass + 0.48*imfw(9)/(mbin_nimf9(10)**(-imfup))*&
             (nslim**(-imfup+1)-mto**(-imfup+1))/(-imfup+1)
        getmass = getmass + 0.077*imfw(9)/(mbin_nimf9(10)**(-imfup))*&
             (nslim**(-imfup+2)-mto**(-imfup+2))/(-imfup+2)

     ELSE IF (mto.GT.mbin_nimf9(9).AND.mto.LE.mbin_nimf9(10)) THEN

        ! 0.9 < MSTO < 1.0

        DO i=1,nimfnp-1
           getmass = getmass + imfw(i)/alpha2(i)*(mbin_nimf9(i+1)**alpha2(i) - &
                mbin_nimf9(i)**alpha2(i))
        ENDDO
        getmass = getmass + imfw(9)/alpha2(i)*(mto**alpha2(i)-mbin_nimf9(9)**alpha2(i)) 
              
        !WD remnants from 1.0-8.5
        getmass = getmass + 0.48*imfw(9)/(mbin_nimf9(10)**(-imfup))*&
             (nslim**(-imfup+1)-mbin_nimf9(10)**(-imfup+1))/(-imfup+1)
        getmass = getmass + 0.077*imfw(9)/(mbin_nimf9(10)**(-imfup))*&
             (nslim**(-imfup+2)-mbin_nimf9(10)**(-imfup+2))/(-imfup+2)
        
        !WD remnants from MSTO-1.0
        getmass = getmass + 0.48*imfw(9)/(1-npi_alphav(i)) * &
             (mbin_nimf9(10)**(1-npi_alphav(i))-mto**(1-npi_alphav(i)))
        getmass = getmass + 0.077*imfw(9)/alpha2(i) * (mbin_nimf9(10)**alpha2(i)-mto**alpha2(i))

     ELSE IF (mto.GT.mbin_nimf9(8).AND.mto.LE.mbin_nimf9(9)) THEN

        ! 0.8 < MSTO < 0.9

        DO i=1,nimfnp-2
           getmass = getmass + imfw(i)/alpha2(i)*(mbin_nimf9(i+1)**alpha2(i) - &
                mbin_nimf9(i)**alpha2(i))
        ENDDO
        getmass = getmass + imfw(8)/alpha2(i)*(mto**alpha2(i)-mbin_nimf9(8)**alpha2(i)) 
        
        !WD remnants from 1.0-8.5
        getmass = getmass + 0.48*imfw(9)/(mbin_nimf9(10)**(-imfup))*&
             (nslim**(-imfup+1)-mbin_nimf9(10)**(-imfup+1))/(-imfup+1)
        getmass = getmass + 0.077*imfw(9)/(mbin_nimf9(10)**(-imfup))*&
             (nslim**(-imfup+2)-mbin_nimf9(10)**(-imfup+2))/(-imfup+2)

        !WD remnants from 0.9-1.0
        getmass = getmass + 0.48*imfw(9)/(1-npi_alphav(i)) * &
             (mbin_nimf9(10)**(1-npi_alphav(i))-mbin_nimf9(9)**(1-npi_alphav(i)))
        getmass = getmass + 0.077*imfw(9)/alpha2(i) * &
             (mbin_nimf9(10)**alpha2(i)-mbin_nimf9(9)**alpha2(i))
        
        !WD remnants from MSTO-0.9
        getmass = getmass + 0.48*imfw(8)/(1-npi_alphav(i)) * &
             (mbin_nimf9(9)**(1-npi_alphav(i))-mto**(1-npi_alphav(i)))
        getmass = getmass + 0.077*imfw(8)/alpha2(i) * (mbin_nimf9(9)**alpha2(i)-mto**alpha2(i))

     ELSE IF (mto.GT.mbin_nimf9(7).AND.mto.LE.mbin_nimf9(8)) THEN

        ! 0.7 < MSTO < 0.8

        DO i=1,nimfnp-3
           getmass = getmass + imfw(i)/alpha2(i)*(mbin_nimf9(i+1)**alpha2(i) - &
                mbin_nimf9(i)**alpha2(i))
        ENDDO
        getmass = getmass + imfw(7)/alpha2(i)*(mto**alpha2(i)-mbin_nimf9(7)**alpha2(i)) 
        
        !WD remnants from 1.0-8.5
        getmass = getmass + 0.48*imfw(9)/(mbin_nimf9(10)**(-imfup))*&
             (nslim**(-imfup+1)-mbin_nimf9(10)**(-imfup+1))/(-imfup+1)
        getmass = getmass + 0.077*imfw(9)/(mbin_nimf9(10)**(-imfup))*&
             (nslim**(-imfup+2)-mbin_nimf9(10)**(-imfup+2))/(-imfup+2)

        !WD remnants from 0.9-1.0
        getmass = getmass + 0.48*imfw(9)/(1-npi_alphav(i)) * &
             (mbin_nimf9(10)**(1-npi_alphav(i))-mbin_nimf9(9)**(1-npi_alphav(i)))
        getmass = getmass + 0.077*imfw(9)/alpha2(i) * &
             (mbin_nimf9(10)**alpha2(i)-mbin_nimf9(9)**alpha2(i))
 
       !WD remnants from 0.8-0.9
        getmass = getmass + 0.48*imfw(8)/(1-npi_alphav(i)) * &
             (mbin_nimf9(9)**(1-npi_alphav(i))-mbin_nimf9(8)**(1-npi_alphav(i)))
        getmass = getmass + 0.077*imfw(8)/alpha2(i) * &
             (mbin_nimf9(9)**alpha2(i)-mbin_nimf9(8)**alpha2(i))
        
        !WD remnants from MSTO-0.8
        getmass = getmass + 0.48*imfw(7)/(1-npi_alphav(i)) * &
             (mbin_nimf9(8)**(1-npi_alphav(i))-mto**(1-npi_alphav(i)))
        getmass = getmass + 0.077*imfw(7)/alpha2(i) * (mbin_nimf9(8)**alpha2(i)-mto**alpha2(i))

     ELSE

        WRITE(*,*) 'GETMASS ERROR, msto<0.7',mto
        STOP

     ENDIF

     !BH remnants
     getmass = getmass + 0.5*10**imf5/(mbin_nimf9(10)**(-imfup))*&
          (imfhi**(-imfup+2)-bhlim**(-imfup+2))/(-imfup+2)
     !NS remnants
     getmass = getmass + 1.4*10**imf5/(mbin_nimf9(10)**(-imfup))*&
          (bhlim**(-imfup+1)-nslim**(-imfup+1))/(-imfup+1)
 
     getmass = getmass / imfnorm

  ENDIF


  IF (PRESENT(timfnorm)) timfnorm = imfnorm

  RETURN

END FUNCTION GETMASS
