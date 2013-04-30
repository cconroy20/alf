FUNCTION IMF(mass)

  !define IMFs (dn/dM)

  !the following is just a fast-and-loose way to pass things around:
  !if the imf_type var is +10 then we calculate dn/dm*m
  !if the imf_type var is <10 then we calculate dn/dm
  
  USE sfvars
  IMPLICIT NONE

  REAL(DP), DIMENSION(:), INTENT(in) :: mass
  REAL(DP), DIMENSION(SIZE(mass))    :: imf
  INTEGER :: i,n
  REAL(DP) :: imfcu, m1=0.08,m2=0.5,m3=1.0

  !---------------------------------------------------------------!
  !---------------------------------------------------------------!
 
  n = SIZE(mass)
  imf(1:n) = 0.0
  
  !Kroupa (2001) IMF
  DO i=1,n
     IF (mass(i).GE.m1.AND.mass(i).LT.m2) &
          imf(i) = mass(i)**(-imf_alpha(1))
     IF (mass(i).GE.m2.AND.mass(i).LT.m3) &
          imf(i) = m2**(-imf_alpha(1)+imf_alpha(2))*&
          mass(i)**(-imf_alpha(2))
     IF (mass(i).GE.m3) &
          imf(i) = m2**(-imf_alpha(1)+imf_alpha(2))*&
          mass(i)**(-imf_alpha(3))
  ENDDO

  IF (imf_flag.EQ.10) imf = mass*imf
  

END FUNCTION IMF
