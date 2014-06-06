      PROGRAM DISCO 	
C
      implicit none
C
c      double precision DR,DT,T,TPRINT,DTPRINT,ex,sca,der1T,derNGT,DRM
      real DR,DT,T,TPRINT,DTPRINT,ex,sca,der1T,derNGT,DRM
      real gamma, RMAX, RINIT,TMAX
      integer I,NG
      PARAMETER(NG=2000)
      PARAMETER(gamma=1.0)
      real, DIMENSION(NG)::SIG,R,RM,RS
      real, DIMENSION(NG)::FP,FM,FP1,FM1,SAG
C
      
      CALL INITF(SIG,R,RM,RS,DR,TMAX,DTPRINT,DRM)
C
      T=1. !equivale a t=0.
      TPRINT=0.
C
      DO WHILE(T.LE.TMAX)
        IF(T.GE.TPRINT) THEN
          TPRINT=TPRINT+DTPRINT
          write(*,*) TPRINT 
          CALL OUTP(SIG,R)
        END IF
        CALL TIMESTEP(DR,DT)
        CALL TSTEP(SIG,SAG,R,T,RM,RS,DR,DT,FP,FM,FP1,FM1,DRM)
        T=T+DT
      END DO
C
	write(*,*), DR,DT,DRM
      STOP
      END
C
      SUBROUTINE OUTP(SIG,R)
C
      implicit none
c
      integer NG,I
      PARAMETER(NG=2000)
      real,DIMENSION(NG)::SIG,R
c      real,DIMENSION(NG)::SAG !se los quito
C
      DO I=1,NG
        WRITE(10,*) R(I),SIG(I)
c        SAG(I)=SIG(I)
      END DO
C
      I=0.0
      RETURN
      END
C
      SUBROUTINE INITF(SIG,R,RM,RS,DR,TMAX,DTPRINT,DRM)
      implicit none
c
c      double precision DR,DTPRINT,ex,sca,DRM
      integer NG, gamma, I
      real DR,DTPRINT,ex,sca,DRM
      real TMAX,RINIT,RMAX
      PARAMETER(NG=2000)
      real,DIMENSION(NG)::SIG,R,RM,RS
c      common/expon/gamma
C
      RMAX=50.   !RMAX/R1 ; R1=10
      RINIT=0.01  !RINIT/R1
      DR=(RMAX-RINIT)/(NG-1) !0.01
      gamma=1.0
c
      DO I=1, NG
        R(I)=RINIT+ (I-1.)*DR
        sca=1./(R(I)**(gamma))
        ex= exp(1.-(R(I)**(2.-gamma)))
        SIG(I)=sca*ex	!sigma prima
c
        DRM=1./DR
c        XNU(I)=1.*(R(I))**gamma
        RM(I)=1./R(I)	
        RS(I)=(SQRT(R(I)))**3.	!gamma+1/2
      END DO
      I=0.0
c
      TMAX=150 !tiempo máximo t/t_s +1, ts=40295 años, Tmax=6 millones de años
      DTPRINT=TMAX/6. !imprime cada millón de años
C
      RETURN
      END
C
      SUBROUTINE TIMESTEP(DR,DT)
C
      implicit none
      integer NG
c      double precision DR,DT
      real DR,DT
      PARAMETER(NG=2000)
C
c      DIMENSION SIG(NG),R(NG)
C
      DT=7.5e-3*DR**2. ! convergencia
C
      RETURN
      END
C
      SUBROUTINE TSTEP(SIG,SAG,R,T,RM,RS,DR,DT,FP,FM,FP1,FM1,DRM)
C
      implicit none
c      double precision T,DR,DT,der1T,derNGT,DRM
      integer I, IM, IP, NG
      real T,DR,DT,der1T,derNGT,DRM
      real gamma
      PARAMETER(NG=2000)
C
      real,DIMENSION(NG)::SIG,R,RM,RS
      real,DIMENSION(NG)::FP,FM,FP1,FM1,SAG
C
      DO I=2,NG-1
        IM=I-1
        IP=I+1
      FP(I)=SQRT(R(I)+0.5*DR)*DRM*(RS(IP)*SIG(IP)
     &      -RS(I)*SIG(I))
      FM(I)=SQRT(R(I)-0.5*DR)*DRM*(RS(I)*SIG(I)
     &      -RS(IM)*SIG(IM))

      END DO
C
      DO I=2,NG-1
        SAG(I)= (0.3e-3)*DR*RM(I) !DRM*DT*3.*RM(I)
        SIG(I)=SIG(I) + SAG(I)*(FP(I)-FM(I)) 
      END DO
c
      der1T=-((T+R(1))/(SQRT(T)**5.*R(1)**2))*exp(1.-(R(1)/T))
c
      derNGT=-((T+R(NG))/(SQRT(T)**5.*R(NG)**2))*exp(1.-(R(NG)/T))

      SIG(1)=SIG(2)-der1T*DR
      SIG(NG)=SIG(NG-1)+derNGT*DR	 
C
      RETURN
      END
