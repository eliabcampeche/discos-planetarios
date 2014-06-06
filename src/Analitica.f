      PROGRAM ANALITICA
      IMPLICIT NONE
c
      DOUBLE PRECISION DR,sca
      REAL gamma,RINIT,RMAX,del 
      INTEGER NG,I,J
      PARAMETER(NG=3000)
      PARAMETER(gamma=1.)
      PARAMETER(del=1.5)
      REAL,DIMENSION(NG):: sig, r
      REAL, DIMENSION(6)::T
c      
      RINIT=0.01
      RMAX=50.
      DR=(RMAX-RINIT)/(NG-1) 
c
	T(1)=1.
	T(2)=26.
	T(3)=51.
	T(4)=76.
	T(5)=101.
	T(6)=126.
c
      DO J=1,6
         DO I=1,NG
            r(I)=RINIT+ (I-1.)*DR
            sca=(1./((T(J)**del)*r(I)**gamma))
            sig(I)=sca*exp(1.-(R(I)**(2.-gamma))/T(J))
            write(11,*) r(I),sig(I)
         ENDDO
      ENDDO
c
      STOP
      END


