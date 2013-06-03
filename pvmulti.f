ccc      implicit real*8 (a-h,o-z)
ccc      parameter (nmax = 32,nd = 32,nsp=30000)
ccc      integer *4 iout(2),inform(10),ierr(10)
ccc      complex*16 ucmplx(nmax),phicmplx(nmax),eye,zero,
ccc     *           qa(2*nmax),u2(2*nd),cfield(2*nmax),zz
ccc      dimension x(nmax),y(nmax),wksp(nsp),x2(2*nmax),
ccc     *          y2(2*nmax)
ccc      dimension xp(nmax),yp(nmax)
cccc
cccc     real and complex constants
cccc
ccc      call prini(6,13)
ccc      pi = 4.d0*datan(1.d0)
ccc      eye = dcmplx(0.d0,1.d0)
ccc      zero = dcmplx(0.d0,0.d0)
ccc      k0 = 0 
ccc      k = 0 
ccc      istart = 0
ccc      x0 = 3.0d0
ccc      y0 = 2.0d0
ccc      do i = 1,nd
ccc         thet = 2*pi/nd
ccc         x(istart+i) = dcos((i-1)*thet)
ccc         y(istart+i) = dsin((i-1)*thet)
ccc         xp(istart+i) = -y(istart+i)
ccc         yp(istart+i) = x(istart+i)
ccc         zz = 1.0d0/dcmplx(x(istart+i)-x0,y(istart+i)-y0)
ccc         ucmplx(i) = zz
ccc         ucmplx(i) = dcmplx(x(istart+i),y(istart+i))
ccc         u2(i) = zz/2
ccc         ucmplx(i) = ucmplx(i)*dcmplx(xp(istart+i),yp(istart+i))
ccc      enddo
ccc      do i = 1,2*nd
ccc         thet = pi/nd
ccc         x2(istart+i) = dcos((i-1)*thet)
ccc         y2(istart+i) = dsin((i-1)*thet)
ccc      enddo
cccccc      call prin2(' x2 is *',x2,2*nd)
cccccc      call prin2(' y2 is *',y2,2*nd)
ccc      call prin2(' u2 is *',u2,2*nd)
cccc
ccc      call prinf(' calling PVINTEV *',nd,0)
ccc      CALL PVINTEV(k0,k,nd,nmax,x,y,x2,y2,
ccc     *                    ucmplx,u2,qa,cfield,nsp,wksp,
ccc     *                    phicmplx)
cccc
ccc      call prin2(' phicmplx = *',phicmplx,2*nd)
ccc      stop 
ccc      end
C
C********x*********x*********x*********x*********x*********x*********x**
C
      SUBROUTINE PVINTEV(K0,K,ND,NMAX,X,Y,X2,Y2,UCMPLX,
     *                   U2,QA,CFIELD,NSP,WKSP,PHICMPLX)
c
c     THIS HAS BEEN MODIFIED TO ALLOW A DIFFERENT NUMBER OF POINTS ON
c     EACH HOLE.
C
C     Calculates Cauchy Principle Value Integrals in bounded
C     and unbounded multiply connected domains.
C
C     The precise analytical form of the integral evaluated is
C
C                       1     /   UCMPLX(r)
C     PHICMPLX(Z) =  ------  |   ----------    dr ,
C                    2 PI i  /    zeta(r) - z    
C                           Gamma
C
C     where r is a (REAL) parametrization of the curve
C     and zeta(r) is the COMPLEX curve coordinate at r. 
C     Thus, PHICMPLX(Z) is not a Cauchy integral. UCMPLX(r)
C     in this case should be the product of the density and
C     d zeta / dr. 
C
C     INPUT PARAMETERS:
C
C     K0 = 0 for bounded domains, 1 for unbounded domains.
C     K  = number of holes.
C     Thus,
C          K0 = 0, K = 0 is a bounded simply connected domain
C          K0 = 0, K = 1 is a bounded doubly connected domain
C          K0 = 1, K = 1 is a domain exterior to a single closed
C                             curve,
C          etc.
C
C     ND(k)  = number of points on hole k
C     NMAX   = total number of points on all holes
C     X,Y    = physical coordinates of boundaries, listed at equispaced
C              points in r, counter clockwise on all boundaries
C     UCMPLX = complex density at equispaced points in r, listed
C              counter clockwise on all boundaries
C
C     WORKSPACE PARAMETERS:
C
C     U2     = a complex array of dimension at least 2*ND
C     QA     = a complex array of dimension at least 2*NMAX
C     CFIELD = a complex array of dimension at least 2*NMAX
C     WKSP   = a complex array of dimension greater than 
C              max(8*NMAX+15,20*NMAX), the first being the
C              length of the workspace needed for Fourier
C              interpolation and the second for the fast 
C              multipole routine.
C
C     OUTPUT PARAMETERS:
C
C     PHICMPLX = value of principal value integral at each point
C                on the boundary
C
C----------------------------------------------------------------------
C     Method:
C     We use the standard alternating point method for evaluating
C     principal value integrals. In order not to lose resolution,
C     however, we first interpolate the density and curve to
C     twice the original number of points.
C----------------------------------------------------------------------
C
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER *4 IOUT(2),INFORM(10),IERR(10),nd(k0:k)
      COMPLEX*16 UCMPLX(NMAX),PHICMPLX(NMAX),EYE,ZERO,
     *           QA(2*NMAX),U2(2*NMAX),CFIELD(2*NMAX)
      DIMENSION X(NMAX),Y(NMAX),WKSP(NSP),X2(2*NMAX),
     *          Y2(2*NMAX),POTEN(1)
       REAL *4  TIMEP(2),ETIME
c
c    real and complex constants
c
      PI = 4.D0*DATAN(1.D0)
      EYE = DCMPLX(0.D0,1.D0)
      ZERO = DCMPLX(0.D0,0.D0)
      if (nsp.lt.40*nmax) write (6,*) 'NOT ENOUGH WKSP IN PVINTEV'
C
C     set QA to zero.
C
ccc      call prinf(' setting qa to 0 *',j,0)
      DO J = 1,2*NMAX
         QA(J) = 0.0D0
      ENDDO
C
C     Interpolate density, x and y coordinates
C
      ISTART = 0
      ISTART2 = 0
      DO NBOD = K0,K
         ND2 = 2*ND(nbod)
         NDM1= ND(nbod)-1
         ND2M1= ND2-1
ccc         CALL PRINF(' CALLING FTRPINC *',J,0)
         CALL FTRPINC(WKSP,NSP,IP2,IPT,NDM1,ND2M1)
ccc         CALL PRINF(' CALLING FINTERC *',J,0)
         CALL FINTERC(UCMPLX(ISTART+1),U2,NDM1,ND2M1,WKSP,NSP,
     1                IP2,IPT)
ccc         CALL PRIN2(' UCMPLX = *',UCMPLX(ISTART+1),2*ND)
ccc         CALL PRIN2(' U2 = *',U2,2*ND2)
         DO I = 2,ND2,2
            QA(ISTART2+I) = 4.d0*PI*U2(I)/ND2
         ENDDO
ccc         CALL PRIN2(' QA = *',QA,2*ND2)
         ISTART = ISTART+ND(nbod)
         ISTART2 = ISTART2+ND2
      END DO 
C
C     set parameters for FMM routine DAPIF2
C
      IOUT(1) = 0
      IOUT(2) = 0
c      
      IFLAG7 = 3
      NAPB = 30
      NINIRE = 2
      TOL = 1.0d-14
      MEX = 300
      EPS7 = 1.0d-14
      NNN = 2*NMAX
      CALL PRINF(' NNN = *',NNN,1)
       T0 = ETIME(TIMEP)
      CALL DAPIF2 (IOUT,IFLAG7,NNN,NAPB,NINIRE,MEX,IERR,INFORM,
     *             TOL,EPS7,X2,Y2,QA,POTEN,CFIELD,WKSP,NSP,CLOSE)
      if (IERR(1).ne.0) then
         write(6,*) '  ERROR IN DAPIF2, IERR = *', (ierr(ii),ii=1,6)
         write(6,*) '  INFORM = *', (inform(ii),ii=1,6)
         stop
      end if
        T1 = ETIME(TIMEP)
         TSEC = T1-T0
       write(6,*)' TIME for DAPIF is ',TSEC
ccc      CALL PRIN2(' CFIELD = *',CFIELD,2*ND2)
C
C     extract imaginary part.
C
      ISTART = 0
      ISTART2 = 0
      DO NBOD = K0,K
         DO I = 1,ND(nbod)
            I2 = 2*I-1
            PHICMPLX(ISTART+I) = -CFIELD(ISTART2+I2)/(2*PI*EYE)
         ENDDO
         ISTART = ISTART+ND(nbod)
         ISTART2 = ISTART2+2*nd(nbod)
      ENDDO
C
      RETURN
      END
C
      SUBROUTINE FINTER(FCORSE,FFINE,NCORSE,NFINE,WORK,LW,IP2,IPT)
C
C     Fourier interpolation subroutine for periodic function.
C
C     INPUT:
C
C     fcorse   = array of function values on coarse grid
C     ncorse+1 = number of coarse grid points
C                IMPORTANT: ncorse is assumed to be odd
C     work     = workspace of length at least
C                4*ncorse + 6*nfine + 50
C     nfine    = number of fine grid points, with
C                (nfine+1) an integral multiple of (ncorse+1)
C
C     OUTPUT:
C
C     ffine = array of function values on fine grid with
C             nfine points. The zero entry positions of
C             the two arrays are chosen to coincide
C
      IMPLICIT REAL *8 (A-H,O-Z)
      INTEGER *4 NCORSE,NFINE,LW,IPT,IP2
      REAL *8 FCORSE(0:NCORSE),FFINE(0:NFINE),WORK(LW)
C
      NC = NCORSE+1
      NF = NFINE+1
c
c---- complexify FCORSE data
c
      DO 100 I = 0,NCORSE
	 WORK(IPT+2*I) = FCORSE(I)
	 WORK(IPT+2*I+1) = 0.0D0
100   CONTINUE
c
      CALL DCFFTF(NC,WORK(IPT),WORK(1))
ccc   call prin2(' work(ipt) = *',work(ipt),2*nc)
      DO 120 I = 0,NCORSE+2
	 WORK(IPT+I) = WORK(IPT+I)/NC
120   CONTINUE
c
      DO 140 I = 0,NCORSE-1
	 WORK(IPT+2*NFINE+1-I) = WORK(IPT+2*NCORSE+1-I)/NC
140   CONTINUE
      do 160 i = NCORSE+3,2*NFINE-NCORSE+1
	 WORK(IPT+I) = 0
160   CONTINUE
      CALL DCFFTB(NF,WORK(IPT),WORK(IP2))
      DO 200 I = 0,NFINE
	 FFINE(i) = WORK(IPT+2*I)
200   CONTINUE
      RETURN
      END
c
      SUBROUTINE FTRPIN(WORK,LW,IP2,IPT,NCORSE,NFINE)
      IMPLICIT REAL *8 (A-H,O-Z)
      INTEGER *4 NCORSE,NFINE,LW,IPT,IP2
      REAL *8 WORK(LW)
C
C---- initialize the two fft work arrays
C
      IP2 = 4*NCORSE + 20
      IPT = IP2 + 4*NFINE + 20
      CALL DCFFTI(NCORSE+1,WORK(1))
      CALL DCFFTI(NFINE+1,WORK(IP2))
      RETURN
      END
C
      SUBROUTINE FINTERC(FCORSE,FFINE,NCORSE,NFINE,WORK,LW,IP2,IPT)
C
C     Fourier interpolation subroutine for periodic function.
C
C     INPUT:
C
C     fcorse   = array of function values on coarse grid
C     ncorse+1 = number of coarse grid points
C                IMPORTANT: ncorse is assumed to be odd
C     work     = workspace of length at least
C                4*ncorse + 6*nfine + 50
C     nfine    = number of fine grid points, with
C                (nfine+1) an integral multiple of (ncorse+1)
C
C     OUTPUT:
C
C     ffine = array of function values on fine grid with
C             nfine points. The zero entry positions of
C             the two arrays are chosen to coincide
C
      IMPLICIT REAL *8 (A-H,O-Z)
      INTEGER *4 NCORSE,NFINE,LW,IPT,IP2
      COMPLEX *16 FCORSE(0:NCORSE),FFINE(0:NFINE),WORK(LW)
C
      NC = NCORSE+1
      NF = NFINE+1
c
c---- complexify FCORSE data
c
      DO 100 I = 0,NCORSE
	 WORK(IPT+I) = FCORSE(I)
100   CONTINUE
c
      CALL DCFFTF(NC,WORK(IPT),WORK(1))
ccc      call prin2(' work(ipt) = *',work(ipt),2*nc)
      DO 120 I = 0,NC/2
	 WORK(IPT+I) = WORK(IPT+I)/NC
120   CONTINUE
c
      DO 140 I = 0,(NC/2)-2
	 WORK(IPT+NFINE-I) = WORK(IPT+NCORSE-I)/NC
140   CONTINUE
      do 160 i = (NC/2)+1,NFINE-(NC/2)+1
	 WORK(IPT+I) = 0
160   CONTINUE
      CALL DCFFTB(NF,WORK(IPT),WORK(IP2))
      DO 200 I = 0,NFINE
	 FFINE(I) = WORK(IPT+I)
200   CONTINUE
      RETURN
      END
c
      SUBROUTINE FTRPINC(WORK,LW,IP2,IPT,NCORSE,NFINE)
      IMPLICIT REAL *8 (A-H,O-Z)
      INTEGER *4 NCORSE,NFINE,LW,IPT,IP2
      COMPLEX *16 WORK(LW)
C
C---- initialize the two fft work arrays
C
      IP2 = 2*NCORSE + 20
      IPT = IP2 + 2*NFINE + 20
ccc      call prinf(' calling DCFFTI *',ncorse,1)
      NC = NCORSE+1
      NF = NFINE+1
      CALL DCFFTI(NC,WORK(1))
      CALL DCFFTI(NF,WORK(IP2))
      RETURN
      END

