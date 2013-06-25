      PROGRAM VORTEX_MOTION
c     ------------------------------------------------------------------
c
c
c  Computes the motion of point vortices in the presence of multiple 
c  islands
c
c  the problem is set up so that there is circulation around the island 
c  at the north pole, only
c 
c  Solve:
c     laplace_beltrami(u) = 0 (except at point vortex locations)
c       u = 0 on \gamma_1
c       u = A_k on \gamma_k
c
c Vortex position is advanced by forward Euler method
c
c The program outputs various matlab files for plotting:
c  geo_3d.m		shows geometry on the sphere
c  geo_stereo.m		shows geometry in stereographic plane
c  targets_stereo.m	location of target points for checking accuracy
c  vort_path.m		shows vortex path (not sure it works for 
c			multiple vortices)
c  vortex_stereo.m	initial location of vortices in stereographic plane
c
c In addition, run surgrid.m to plot solution on grid
c This requires the following files:
c  igrid.dat - shows whether grid point in or out of domain
c  xgrid.dat, ygrid.dat, zgrid.dat - surface grid on sphere
c  xzeta_grid.dat, yzeta_grid.dat - grid in stereographic plane
c  ugrid.dat - solution on grid
c
c     ------------------------------------------------------------------
      implicit real*8 (a-h, o-z)
      parameter (kmax = 10, npmax = 2048, nmax = kmax*npmax)
c
c Geometry of holes
      dimension ak(kmax), bk(kmax), th_k(kmax), phi_k(kmax), cx(kmax),
     1          cy(kmax), cz(kmax)
      complex*16 zeta_k(kmax)
c
c Information about source points on the sphere
      dimension xs(nmax), ys(nmax), zs(nmax), dx(nmax), dy(nmax),
     1          dz(nmax), d2x(nmax), d2y(nmax), d2z(nmax)
c
c Information about source points in stereographic plane  
c (i.e. the discretization points)
      dimension diag(nmax), x_zeta(nmax), y_zeta(nmax)
      complex*16 dzeta(nmax), zeta(nmax)
c
c  Grid variables
      parameter (nth_max = 1000, nphi_max = 1000, 
     1          ng_max = nth_max*nphi_max)
      dimension igrid(ng_max), th_gr(ng_max), phi_gr(ng_max),
     1          u_gr(ng_max), x_gr(ng_max), y_gr(ng_max), z_gr(ng_max),
     2          xzeta_gr(ng_max), yzeta_gr(ng_max), 
     3          alph_gr(ng_max)
      complex*16 zeta_gr(ng_max)
c
c target points are used to check accuracy
      parameter (ntar_max = 1000)
      dimension xz_tar(ntar_max), yz_tar(ntar_max), u_tar(ntar_max)
      complex*16 zeta_tar(ntar_max) 
c
c boundary conditions
      dimension rhs(nmax), A_k(kmax)   
c
c Point Vortices
      parameter (nvortmax = 100)
      dimension vort_k(nvortmax), x1_vort(nvortmax), x2_vort(nvortmax), 
     1          x3_vort(nvortmax)
      complex*16 zk_vort(nvortmax), zvel(nvortmax)
c
c  Matrix equation variables for GMRES
c  MAXL is the maximum nubmer of GMRES iterations performed
c       before restarting.
c  LRWORK is the dimension of a real workspace needed by DGMRES.
c  LIWORK is the dimension of an integer workspace needed by DGMRES.
c  GMWORK and IGWORK are work arrays used by DGMRES
c
      parameter (maxl = 50,liwork=30,  
     1           lrwork=10+(nmax+kmax)*(maxl+6)+maxl*(maxl+3))
      dimension gmwork(lrwork), igwork(liwork),density(nmax)
c
c Fast Multipole Arrays
c Location of dipoles - dipvec(normal to dz) 
c pottarg is the potential at the target,
c gradtarg, its gradient and hesstarg, its gradient
c pot, grad, hess - at the source points
c Charges , qa are essentially 0 
c dipstr - dipole strength is essentially density

      dimension xat(nmax+ng_max), yat(nmax+ng_max)
      complex*16 qa(nmax+ng_max),
     $           pot(nmax+ng_max), grad(2*nmax+ng_max),
     $           hess(3*nmax+ng_max), pottarg(nth_max*nphi_max),
     $		 gradtarg(2*(nth_max*nphi_max)), hesstarg(3*(nth_max*nphi_max)),
     $		 dipstr(nmax+ng_max)
	
      dimension dipvec(2*(nmax+ng_max)),
     $          source(2*(nmax+ng_max)),targ(2*(nth_max*nphi_max))
     
c
c Logicals
      logical make_movie, debug
c
c Other arrays
      dimension alpha(nmax), w(nmax), u(nmax)
      REAL*4 TIMEP(2), ETIME
c
c common blocks
      common /geometry/ x_zeta, y_zeta, zeta, dzeta
      common /inteqn/ diag, zeta_k
      common /sys_size/ k, nd, nbk
      common /fasblk2/ schur,wb,ipvtbf
c
c Open output file for vortex path
         open (unit = 41, file = 'vort_path.m')
c
c if making a movie (this is not debugged right now!)
         make_movie = .false.
         if (make_movie) then
            open (unit = 37, file = 'movie/u_movie.m')
         end if
c
c set debug = .true. if running the spherical cap case
         debug = .true.
c
c Read hole geometry data
         call READ_DATA (k, nd, nbk, nth, nphi, ak, bk, cx, cy, cz, 
     1                   th_k, phi_k, nvort, x1_vort, x2_vort, x3_vort,
     2                   vort_k, gamma_tot, r0, zeta_k, dt, ntime)
c
c Construct boundary geometry on surface of sphere
         call MAKE_GEO (k, nd, nbk, ak, bk, th_k, phi_k, xs, ys, zs,
     1                  dx, dy, dz, d2x, d2y, d2z)
c
c Get stereo graphic projection
         call STEREO (k, nd, nbk, xs, ys, zs, dx, dy, dz, d2x, d2y, d2z,
     1                zeta, dzeta, x_zeta, y_zeta, diag, nvort, x1_vort, 
     2                x2_vort, x3_vort, zk_vort) 
         call RSCPLOT (zk_vort, nvort, 1, 41)
c
c Construct grid on surface of sphere
         call SURFACE_GRID (k, nd, nbk, nth, nphi, ak, bk, th_k, phi_k,  
     1                      th_gr, phi_gr, x_gr, y_gr, z_gr, zeta_gr, 
     2                      xzeta_gr, yzeta_gr, igrid, alph_gr)
         if (debug) then
            call TARGET_POINTS (k, nd, nbk, ak, bk, th_k, phi_k, ntar,
     2                          xz_tar, yz_tar, zeta_tar)
            call prinf (' ntar = *', ntar, 1)
            call prin2 (' xz_tar = *', xz_tar, ntar)
            call prin2 (' yz_tar = *', yz_tar, ntar)
            ntime = 1
         end if
c
c Time loop for vortex path
         tbeg = etime(timep)
         do it = 1, ntime
            time = it*dt  
            call PRIN2 (' TIME = *', time, 1)       
c
c Construct the RHS and solve
            call PRINI (6,13)
            call GETRHS (k, nd, nbk, zeta_k, zeta, rhs, nvort, vort_k, 
     1                   zk_vort, gamma_tot)
            call PRIN2 (' rhs = *', rhs, nbk)
            call SOLVE (nd, k, kmax, nbk, rhs, density, A_k,  
     1                  gmwork, lrwork, igwork, liwork, maxl)
            call prin2 (' density = *', density, nbk)
c
c Construct solution on surface grid
ccc         call SOL_GRID (nd, k, nbk, nth, nphi, density, A_k, zeta_k,   
ccc     1                  zeta, dzeta, igrid, zeta_gr, u_gr)
         nplot = mod(it,100)
         call PRINF (' nplot = *', nplot, 1)
ccc         if (mod(it,100).eq.0) then
         call SOL_GRID_FMM (nd, k, nbk, nth, nphi, density, zeta_k,   
     1                      zeta, dzeta, igrid, zeta_gr, u_gr,
     2                      qa,dipstr,grad,pot,gradtarg,pottarg,hess,
     3		            hesstarg,source,targ,dipvec, 
     4                      nvort, vort_k, zk_vort, gamma_tot)
         if (debug) then
            call SOL_TAR_FMM (nd, k, nbk, ntar, density, zeta_k,   
     1                       zeta, dzeta, xz_tar, yz_tar, u_tar, 
     2                       qa,dipstr, grad, pot,pottarg,gradtarg,hess,
     3			     hesstarg, source,targ,dipvec,nvort, 
     4                       vort_k, zk_vort, gamma_tot)
c
c          for a vortex in presence of cap with radius r0, check solution
            call CHECK_ERROR_TAR (nd, k, nbk, ntar, zeta_k, zeta_tar,  
     1                            u_tar, nvort, vort_k, zk_vort, r0)
         end if
         if (make_movie) then
            call DUMP_MOVIE_ALL (nth, nphi, time, u_gr, it, 37)
            call DUMP_MOVIE_VORT (nth, nphi, time, zk_vort(1), u_gr, 
     1                            it, 37)
         end if
c
c Calculate velocity at a point
         call CALC_VEL (k, nd, nbk, nvort, density, gamma_tot, zeta,  
     1                  dzeta, zeta_k, vort_k, zk_vort, zvel)
         do ivort = 1, nvort
            zk_vort(ivort) = zk_vort(ivort) + dt*zvel(ivort)
         end do
         call PRIn2 (' zk_vort = *', zk_vort, 2*nvort)
         if (mod(it,1).eq.0) then
            call RSCPLOT (zk_vort, nvort, 1, 41)
         end if
         end do
         tend = etime(timep)
         call PRIN2 (' TOTAL CPU TIME = *', tend-tbeg, 1)
c
      stop
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine READ_DATA (k, nd, nbk, nth, nphi, ak, bk, cx, cy, cz, 
     1                      th_k, phi_k, nvort, x1_vort, x2_vort, 
     2                      x3_vort, vort_k, gamma_tot, r0, zeta_k, dt,
     3                      ntime)
c---------------
c OUTPUT
c	k	= number of contours
c	nd	= number of points per contour
c	nbk	= total size of system	
c	nth,nphi
c		= dimensions of grid for plotting
c	ak,bk	= major/minor axis of ellipse
c	(th_k,	= location of hole centre in spherical coordinates
c       phi_k)
c	(cx,cy,cz)
c		= location of centre in physical coordinates
c       vort_k	= strength of kth vortex
c	gamma_tot
c		= sum of vortex strengths
c	r0	= radius of spherical cap (for debugging)
c	zeta_k	= hole centres in stereographic plane
c	dt	= time step size	
c	ntime	= number of time steps to take
c
      implicit real*8 (a-h,o-z)
      dimension ak(*), bk(*), cx(*), cy(*), cz(*), th_k(*), phi_k(*)
      dimension x1_vort(*), x2_vort(*), x3_vort(*), vort_k(*)
      complex*16 zeta_k(*), eye
c
         eye = dcmplx(0.d0,1.d0)
c
         open (unit = 12, file = 'input_vort.data')
         call PRINI (6,13)
c
         read (12,*) k, nd, nvort
         read (12,*) dt, ntime
         nbk = k*nd
         call PRINF (' nbk = *', nbk, 1)
         call PRINF (' Number of time steps = *', ntime, 1)
         call PRIN2 (' dt = *', dt, 1)
         read(12,*) nth, nphi
         do kbod = 1, k
            read(12,*) ak(kbod), bk(kbod), th_k(kbod), phi_k(kbod)
            call SPH2CART (th_k(kbod),phi_k(kbod), 1.d0, cx(kbod),
     1                     cy(kbod), cz(kbod))
         end do
         call PRINF (' nvort = *', nvort, 1)
         gamma_tot = 0.d0
         do ivort = 1, nvort
            read(12,*) theta, phi, vort_k(ivort)
            call SPH2CART (theta, phi, 1.d0, x1_vort(ivort),
     1                     x2_vort(ivort), x3_vort(ivort))
            gamma_tot = gamma_tot + vort_k(ivort)
         end do
         close(12)
c
c initialize plotting outputs (6 is to screen, 13 is to fort.13)
         call PRINI (6,13) 
         call PRIN2 (' ak = *', ak, k)
         call PRIN2 (' bk = *', bk, k)
         call PRIN2 (' cx = *', cx, k)
         call PRIN2 (' cy = *', cy, k)
         call PRIN2 (' cz = *', cz, k)
         call PRINF (' nth = *', nth, 1)
         call PRINF (' nphi = *', nphi, 1)
         call PRIN2 ('    vort_k = *', vort_k, nvort)
         call PRIN2 ('    x1_vort = *', x1_vort, nvort)
         call PRIN2 ('    x2_vort = *', x2_vort, nvort)
         call PRIN2 ('    x3_vort = *', x3_vort, nvort)
         r0 = ak(1)/(1.d0-dsqrt(1.d0-(ak(1))**2))
         call PRIN2 (' r0 = *', r0, 1)
c
c be careful if one of the hole centres is at the north pole
c if it is, nudge it a little
         eps = 1.d-6
         do kbod = 1, k
            check = dabs(cz(kbod)-1.d0)
            if (check.lt.eps) then
               cz_s = 0.999d0
               cx_s = dsqrt(0.5d0*(1-cz_s**2))
               cy_s = cx_s
               zeta_k(kbod) = (cx_s + eye*cy_s)/(1.d0-cz_s)
               write (6,*) 'Fudging hole centre a little'
               call prin2 (' old centre, cx = *', cx(kbod), 1)
               call prin2 (' old centre, cy = *', cy(kbod), 1)
               call prin2 (' old centre, cz = *', cz(kbod), 1)
               call prin2 (' new centre, cx = *', cx_s, 1)
               call prin2 (' new centre, cy = *', cy_s, 1)
               call prin2 (' new centre, cz = *', cz_s, 1)
            else
               zeta_k(kbod) = (cx(kbod) + eye*cy(kbod))/(1.d0-cz(kbod))
            end if
         end do
         call PRIN2 (' zeta_k = *', zeta_k, 2*k)
c
c
      return
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine SPH2CART (theta, phi, r, x, y, z)
c---------------
      implicit real*8 (a-h,o-z)
c
         x = r*dcos(phi)*dcos(theta)
         y = r*dcos(phi)*dsin(theta)
         z = r*dsin(phi)
c
      return
      end      
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine R_FUNC (alpha, A, B, th, phi, x, y, z)
c---------------
      implicit real*8 (a-h,o-z)
      dimension zaxis(3), xaxis(3), yaxis(3)
c
         pi = 4.d0*datan(1.d0)
c
         call SPH2CART (th, phi, 1.d0, zaxis(1), zaxis(2), zaxis(3))
         call SPH2CART (th, phi-0.5d0*pi, 1.d0, xaxis(1), xaxis(2), 
     1                  xaxis(3))
         call CROSS (zaxis, xaxis, yaxis)
         xp = A*dcos(alpha)
         yp = B*dsin(alpha)
         zp = dsqrt(1.d0 - xp**2 - yp**2)
         x = xp*xaxis(1) + yp*yaxis(1) + zp*zaxis(1)
         y = xp*xaxis(2) + yp*yaxis(2) + zp*zaxis(2)
         z = xp*xaxis(3) + yp*yaxis(3) + zp*zaxis(3)
c
      return
      end      
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine CROSS (u, v, w)
c---------------
      implicit real*8 (a-h,o-z)
      dimension u(3), v(3), w(3)
c
         w(1) = u(2)*v(3) - v(2)*u(3)
         w(2) = v(1)*u(3) - u(1)*v(3)
         w(3) = u(1)*v(2) - v(1)*u(2)
c
      return
      end      
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine DOT (u, v, u_dot_v)
c---------------
      implicit real*8 (a-h,o-z)
      dimension u(3), v(3)
c
         u_dot_v = 0.d0
         do i = 1, 3
            u_dot_v = u_dot_v + u(i)*v(i)
         end do
c
      return
      end      
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine DR_FUNC (alpha, A, B, th, phi, dx, dy, dz)
c---------------
c DR_FUNC and D2R_FUNC should really not be in two separate subroutines
c
      implicit real*8 (a-h,o-z)
      dimension zaxis(3), xaxis(3), yaxis(3)
c
         pi = 4.d0*datan(1.d0)
c
         call SPH2CART (th, phi, 1.d0, zaxis(1), zaxis(2), zaxis(3))
         call SPH2CART (th, phi-0.5d0*pi, 1.d0, xaxis(1), xaxis(2), 
     1                  xaxis(3))
         call CROSS (zaxis, xaxis, yaxis)
         xp = A*dcos(alpha)
          dxp = -A*dsin(alpha)
         yp = B*dsin(alpha)
          dyp = B*dcos(alpha)
         zp = dsqrt(1.d0 - xp**2 - yp**2)
          dzp = (-xp*dxp-yp*dyp)/zp
         dx = dxp*xaxis(1) + dyp*yaxis(1) + dzp*zaxis(1)
         dy = dxp*xaxis(2) + dyp*yaxis(2) + dzp*zaxis(2)
         dz = dxp*xaxis(3) + dyp*yaxis(3) + dzp*zaxis(3)
c
      return
      end      
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine D2R_FUNC (alpha, A, B, th, phi, d2x, d2y, d2z)
c---------------
      implicit real*8 (a-h,o-z)
      dimension zaxis(3), xaxis(3), yaxis(3)
c
         pi = 4.d0*datan(1.d0)
c
         call SPH2CART (th, phi, 1.d0, zaxis(1), zaxis(2), zaxis(3))
         call SPH2CART (th, phi-0.5d0*pi, 1.d0, xaxis(1), xaxis(2), 
     1                  xaxis(3))
         call CROSS (zaxis, xaxis, yaxis)
         xp = A*dcos(alpha)
          dxp = -A*dsin(alpha)
          d2xp = -A*dcos(alpha)
         yp = B*dsin(alpha)
          dyp = B*dcos(alpha)
          d2yp = -B*dsin(alpha)
         zp = dsqrt(1.d0 - xp**2 - yp**2)
          dzp = (-xp*dxp-yp*dyp)/zp
          d2zp = (-dxp**2 - xp*d2xp - dyp**2 - yp*d2yp - dzp**2)/zp
         d2x = d2xp*xaxis(1) + d2yp*yaxis(1) + d2zp*zaxis(1)
         d2y = d2xp*xaxis(2) + d2yp*yaxis(2) + d2zp*zaxis(2)
         d2z = d2xp*xaxis(3) + d2yp*yaxis(3) + d2zp*zaxis(3)
c
      return
      end      
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine MAKE_GEO (k, nd, nbk, ak, bk, th_k, phi_k, xs, ys, zs,
     1                     dx, dy, dz, d2x, d2y, d2z)
c---------------
c INPUT
c	k	= number of contours
c	nd	= number of points per contour
c	nbk	= total size of system	
c	ak,bk	= major/minor axis of ellipse
c	(th_k,	= location of hole centre in spherical coordinates
c       phi_k)
c OUTPUT
c	(xs,ys,zs)	
c		= coordinates of each point on boundary
c	(dx,dy,dz)			
c		= derivative wrt parametrization
c	(d2x,d2y,d2z)
c		= 2nd derivative wrt parametrization
c
      implicit real*8 (a-h,o-z)
      dimension ak(k), bk(k), th_k(k), phi_k(k), d2x(nbk), d2y(nbk), 
     1          d2z(nbk)
      dimension xs(nbk), ys(nbk), zs(nbk), dx(nbk), dy(nbk), dz(nbk)
c
         pi = 4.d0*datan(1.d0)
c
         dalph = 2.d0*pi/nd
         istart = 0
         do kbod = 1, k
            do i = 1, nd
               alpha = dalph*(i-1.d0)
               call R_FUNC (alpha, ak(kbod), bk(kbod), th_k(kbod), 
     1                      phi_k(kbod), xs(istart+i), ys(istart+i), 
     2                      zs(istart+i))
               call DR_FUNC (alpha, ak(kbod), bk(kbod), th_k(kbod), 
     1                      phi_k(kbod), dx(istart+i), dy(istart+i), 
     2                      dz(istart+i))
               call D2R_FUNC (alpha, ak(kbod), bk(kbod), th_k(kbod), 
     1                        phi_k(kbod), d2x(istart+i), d2y(istart+i),
     2                        d2z(istart+i))
            end do
            istart = istart + nd
         end do
c
c dump out for matlab plotting
         open (unit = 42, file = 'geo_3d.m')
         is = 1
         do kbod = 1, k
            call RS_3D_PLOT (xs(is),ys(is),zs(is), nd, 1, 42)
            is = is + nd
         end do
         close(42)
c
      return
      end      
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine IN_OR_OUT (x, y, z, eps, ak, bk, th_k, phi_k, itest)
c---------------
      implicit real*8 (a-h,o-z)
      dimension p(3), x_ax(3), y_ax(3), z_ax(3)
c
         pi = 4.d0*datan(1.d0)
c
         p(1) = x
         p(2) = y
         p(3) = z
         call SPH2CART (th_k, phi_k, 1.d0, z_ax(1), z_ax(2), z_ax(3))
         call SPH2CART (th_k, phi_k-0.5d0*pi, 1.d0, x_ax(1), x_ax(2), 
     1                  x_ax(3))
         call CROSS (z_ax, x_ax, y_ax)
         call DOT (p, x_ax, x1)
         call DOT (p, y_ax, y1)
         call DOT (p, z_ax, z1)
         rad = dsqrt((x1/ak)**2 + (y1/bk)**2)
         if ((rad<1.d0+eps).and.(z1>0.d0)) then 
            itest = 1
           else
            itest = 0
         end if
c
      return
      end      
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine SURFACE_GRID (k, nd, nbk, nth, nphi, ak, bk, th_k, 
     1                         phi_k, th_gr, phi_gr, x_gr, y_gr, z_gr, 
     2                         zeta_gr, xzeta_gr, yzeta_gr, igrid, 
     3                         alph_gr)
c---------------
c Constructs grid points on surface of sphere for plotting solution
c in domain
c INPUT
c	k	= number of contours
c	nd	= number of points per contour
c	nbk	= total size of system	
c	(nth,nphi)
c		= number of grid points in theta and phi directions
c	ak,bk	= major/minor axes of ellipses
c	(th_k,phi_k)	
c		= hole centres in spherical coordinates
c OUTPUT
c	(th_gr,phi_gr)	
c		= (theta,phi) values at grid points
c	(x_gr,y_gr,z_gr)			
c		= (x,y,z) values at grid points (on sphere)
c	zeta_gr	= grid point locations in stereographic plane	
c		= xzeta_gr + i yzeta_gr
c	igrid(i,j)	
c		= 1 if (i,j)th point is in domain, 0 otherwise
c	alph_gr	= ignore this - might be used for matlab plotting
c 
      implicit real*8 (a-h,o-z)
      dimension ak(k), bk(k), th_k(k), phi_k(k)
      dimension igrid(nth,nphi), th_gr(nth,nphi), phi_gr(nth,nphi),
     1          x_gr(nth,nphi), y_gr(nth,nphi), z_gr(nth,nphi),
     2          xzeta_gr(nth,nphi), yzeta_gr(nth,nphi), 
     3          alph_gr(nth,nphi)
      complex*16 zeta_gr(nth,nphi), eye
c
         pi = 4.d0*datan(1.d0)
         eye = dcmplx(0.d0,1.d0)
c
c Calculate epsilon, which determines a buffer zone between boundary
c and grid points considered in the domain (i.e. it ensures target 
c points don't get too close to boundary where quadrature breaks down)
         radmax = 0.d0
         do kbod = 1, k
            radmax = max(radmax, dabs(ak(kbod)))
            radmax = max(radmax, dabs(bk(kbod)))
         end do
         call PRIN2 (' maximum radius = *', radmax, 1)
c
c for accuracy, choose nfac to be 5 or higher. Smaller values are better
c for plotting purposes only
         nfac = 2
         eps = nfac*2.d0*pi*radmax/nd
         call PRIN2 (' Epsilon = *', eps, 1) 
c
         dth = 2.d0*pi/nth
         dphi = pi/nphi
         ntar = 0
         do i = 1, nth
            theta = (i-1)*dth
            do j = 1, nphi
               phi = (j-0.5d0)*dphi - 0.5d0*pi
               th_gr(i,j) = theta
               phi_gr(i,j) = phi
               call SPH2CART (th_gr(i,j), phi_gr(i,j), 1.d0, 
     1                        x_gr(i,j), y_gr(i,j), z_gr(i,j))
               zeta_gr(i,j) = (x_gr(i,j) + eye*y_gr(i,j))/
     1                             (1.d0 - z_gr(i,j))
               xzeta_gr(i,j) = dreal(zeta_gr(i,j))
               yzeta_gr(i,j) = dimag(zeta_gr(i,j))
               in_out = 0
               in_out2 = 0
               do kbod = 1, k
                  call IN_OR_OUT (x_gr(i,j), y_gr(i,j), z_gr(i,j), eps,  
     1                            ak(kbod), bk(kbod), th_k(kbod), 
     2                            phi_k(kbod), itest)
                  call IN_OR_OUT (x_gr(i,j), y_gr(i,j), z_gr(i,j), 0.d0,  
     1                            ak(kbod), bk(kbod), th_k(kbod), 
     2                            phi_k(kbod), itest2)
                  in_out = in_out + itest
                  in_out2 = in_out2 + itest2
               end do
               if (in_out>0) then
                  igrid(i,j) = 0
                 else
                  igrid(i,j) = 1
               end if
               if (in_out2>0) then
                  alph_gr(i,j) = 0.9
                 else
                  alph_gr(i,j) = 0.1
               end if
            end do
         end do 
c
c dump out for matlab plotting
         open (unit = 31, file = 'igrid.dat')
         open (unit = 32, file = 'xgrid.dat')
         open (unit = 33, file = 'ygrid.dat')
         open (unit = 34, file = 'zgrid.dat')
         open (unit = 35, file = 'xzeta_grid.dat')
         open (unit = 36, file = 'yzeta_grid.dat')
            call DUMP (nth, nphi, x_gr, igrid, 0, 31)
            call DUMP (nth, nphi, x_gr, igrid, 1, 32)
            call DUMP (nth, nphi, y_gr, igrid, 1, 33)
            call DUMP (nth, nphi, z_gr, igrid, 1, 34)
            call DUMP (nth, nphi, xzeta_gr, igrid, 1, 35)
            call DUMP (nth, nphi, yzeta_gr, igrid, 1, 36)
         close(31)
         close(32)
         close(33)
         close(34)
         close(35)
         close(36)
c
      return
      end      
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine TARGET_POINTS (k, nd, nbk, ak, bk, th_k, phi_k, ntar,
     1                          xz_tar, yz_tar, zeta_tar)
c---------------
c Same idea as SURFACE_GRID, except far fewer points and epsilon
c is fixed (and relatively large)
c INPUT
c	k	= number of contours
c	nd	= number of points per contour
c	nbk	= total size of system	
c	(nth,nphi)
c		= number of grid points in theta and phi directions
c	ak,bk	= major/minor axes of ellipses
c	(th_k,phi_k)	
c		= hole centres in spherical coordinates
c OUTPUT
c       ntar	= number of target points
c	zeta_tar	
c		= target point locations in stereographic plane	
c		= xzeta_gr + i yzeta_gr
c 
      implicit real*8 (a-h,o-z)
      dimension ak(k), bk(k), th_k(k), phi_k(k)
      dimension xz_tar(*), yz_tar(*)
      complex*16 zeta, eye, zeta_tar(*)
c
         pi = 4.d0*datan(1.d0)
         eye = dcmplx(0.d0,1.d0)
c
c Calculate epsilon
         radmax = 0.d0
         do kbod = 1, k
            radmax = max(radmax, dabs(ak(kbod)))
            radmax = max(radmax, dabs(bk(kbod)))
         end do
         call PRIN2 (' maximum radius = *', radmax, 1)
ccc         eps = 20*2.d0*pi*radmax/nd
         eps = 0.5d0
         call PRIN2 (' Epsilon = *', eps, 1) 
c
         nth = 10
         nphi = 15
         dth = 2.d0*pi/nth
         dphi = pi/nphi
         ntar = 0
         do i = 1, nth
            theta = (i-1)*dth
            do j = 1, nphi
               phi = (j-0.5d0)*dphi - 0.5d0*pi
               call SPH2CART (theta, phi, 1.d0, x, y, z)
               zeta = (x + eye*y)/(1.d0 - z)
               xzeta = dreal(zeta)
               yzeta = dimag(zeta)
               in_out = 0
               do kbod = 1, k
                  call IN_OR_OUT (x, y, z, eps,  
     1                            ak(kbod), bk(kbod), th_k(kbod), 
     2                            phi_k(kbod), itest)
                  in_out = in_out + itest
               end do
               if (in_out.eq.0) then
                  ntar = ntar+1
                  zeta_tar(ntar) = zeta
                  xz_tar(ntar) = xzeta
                  yz_tar(ntar) = yzeta
               end if
            end do
         end do 
         call PRINF (' ntar = *', ntar, 1)
cccc
cccc Can also pick a circular contour within the domain as the target 
cccc points
ccc         ntar = 20
ccc         dth = 2.d0*pi/ntar
ccc         do i = 1, ntar
ccc            theta = dth*(i-1)
ccc            zeta_tar(i) = 0.5d0*cdexp(eye*theta)
ccc            xz_tar(i) = dreal(zeta_tar(i))
ccc            yz_tar(i) = dimag(zeta_tar(i))
ccc         end do
c
         open (unit = 53, file = 'targets_stereo.m')
         call RSCPLOT (zeta_tar, ntar, 1, 53)
         close(53)
c
      return
      end      
c
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine STEREO (k, nd, nbk, xs, ys, zs, dx, dy, dz, d2x, d2y,  
     1                   d2z, zeta, dzeta, x_zeta, y_zeta, diag, 
     2                   nvort, x1_vort, x2_vort, x3_vort, zk_vort)
c---------------
c Constructs geometry of boundaries in stereographic plane
c INPUT
c	k	= number of contours
c	nd	= number of points per contour
c	nbk	= total size of system	
c	(xs,ys,zs)	
c		= coordinates of each point on boundary
c	(dx,dy,dz)			
c		= derivative wrt parametrization
c	(d2x,d2y,d2z)
c		= 2nd derivative wrt parametrization
c	nvort	= number of point vortices
c	(x1_vort,y1_vort,z1_vort)	
c		= coordinates of vortices on sphere
c OUTPUT
c	zeta = x_zeta + i y_zeta
c		= boundary points in stereographic plane
c	dzeta	= derivative of zeta wrt parametrization oriented
c		  carefully for evaluation of Cauchy integrals
c	diag	= self-interacting term in integral operator
c
      implicit real*8 (a-h,o-z)
      dimension xs(nbk), ys(nbk), zs(nbk), dx(nbk), dy(nbk), dz(nbk),
     1          x_zeta(nbk), y_zeta(nbk), d2x(nbk), d2y(nbk), d2z(nbk),
     2          diag(nbk), x1_vort(nvort), x2_vort(nvort), 
     3          x3_vort(nvort)
      complex*16 zeta(nbk), dzeta(nbk), eye, d2zeta, zextra, 
     1           zk_vort(nvort)
c
         pi = 4.d0*datan(1.d0)
         eye = dcmplx(0.d0,1.d0)
c
         do i = 1, nbk
            zeta(i) = (xs(i) + eye*ys(i))/(1.d0-zs(i))
            x_zeta(i) = dreal(zeta(i))
            y_zeta(i) = dimag(zeta(i))
            dzeta(i) = (dx(i) + eye*dy(i) + zeta(i)*dz(i))/(1.d0-zs(i))
            xdot = dreal(dzeta(i))
            ydot = dimag(dzeta(i))
            ds = cdabs(dzeta(i))
            d2zeta = (d2x(i) + eye*d2y(i) + 2.d0*dzeta(i)*dz(i) + 
     1                zeta(i)*d2z(i))/(1.d0-zs(i))
            xddot = dreal(d2zeta)
            yddot = dimag(d2zeta)
            rkappa = (xdot*yddot-ydot*xddot)/ds**3
            zextra = dzeta(i)*dconjg(zeta(i))/(1.d0+cdabs(zeta(i))**2)
            zextra = zextra/(2.d0*pi)
            diag(i) = 0.25d0*rkappa*ds/pi - dimag(zextra)
         end do
         do ivort = 1, nvort
            zk_vort(ivort) = (x1_vort(ivort) + eye*x2_vort(ivort))/
     1                       (1.d0 - x3_vort(ivort))
         end do
c
c dump out for matlab plotting
         open (unit = 11, file = 'geo_stereo.m')
         open (unit = 22, file = 'vortex_stereo.m')
         call RSCPLOT (zk_vort, nvort, 1, 22)
         do kbod = 1, k
            call RSCPLOT (zeta((kbod-1)*nd+1), nd, 1, 11)
         end do
         close(11)
         close(22)
c
      return
      end      
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine POINT_VORTEX (zeta, zeta_k, psi)
c---------------
c  evaluates G(zeta,zeta_k) from paper (i.e. Green's function in 
c  complex plane)
c
      implicit real*8 (a-h,o-z) 
      complex*16 zeta, zeta_k, zdis 
c   
         pi = 4.d0*datan(1.d0)
c
         zdis = zeta - zeta_k
         az1 = (cdabs(zeta))**2
         az2 = (cdabs(zeta_k))**2
         arg = 2.d0*dreal(zdis*conjg(zdis)/((1.d0+az1)*(1.d0+az2)))
         psi = -dlog(arg)/(4.d0*pi)
c
      return
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine GRAD_POINT_VORTEX (zeta, zeta_k, zgrad)
c---------------
c  evaluates d G(zeta,zeta_k)/dz for velocity calculation 
c
      implicit real*8 (a-h,o-z) 
      complex*16 zeta, zeta_k, zdis , zgrad
c   
         pi = 4.d0*datan(1.d0)
c
         zdis = zeta - zeta_k
         dis = cdabs(zdis)
         az1 = (cdabs(zeta))**2
         zgrad = -1.d0/zdis + dconjg(zeta)/(1.d0+az1)
         zgrad = zgrad/(4.d0*pi)
c
      return
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine GETRHS (k, nd, nbk, zeta_k, zeta, rhs, nvort, vort_k, 
     1                   zk_vort, gamma_tot)
c---------------
c represent psi = psi^* + sum vortices + A G(zeta,zeta_k(1))
c boundary conditions for psi^* are 
c    psi^* = -sum vortices - A G(zeta,zeta_k(1))
c INPUT
c	k	= number of contours
c	nd	= number of points per contour
c	nbk	= total size of system	
c       zeta_k	= hole centres in stereographic plane
c	nvort	= number of vortices
c	vort_k	= vortex strength
c	zk_vort	= vortex locations in stereographic plane
c	gamma_tot
c		= sum of vortex strengths
c OUTPUT
c	rhs	= boundary conditions according to expression above
c
      implicit real*8 (a-h,o-z)
      dimension rhs(nbk), vort_k(nvort)
      complex*16 zeta_k(k), eye, zeta(nbk), zk_vort(nvort)
c
         eye = dcmplx(0.d0,1.d0)
c
         A = -gamma_tot
ccc         call PRIN2 (' in GETRHS, zk_vort = *', zk_vort, 2*k)
ccc         call PRIN2 ('            vort_k = *', vort_k, k)
ccc         call PRIN2 ('            gamma_tot = *', gamma_tot, 1)
         istart = 0
         do kbod = 1, k
            do j = 1, nd
               psi_vort = 0.d0
               call POINT_VORTEX (zeta(istart+j), zeta_k(1), circ)
               do ivort = 1, nvort
                  call POINT_VORTEX (zeta(istart+j), zk_vort(ivort), 
     1                               psi)
                  psi_vort = psi_vort + vort_k(ivort)*psi
               end do
               rhs(istart+j) =  - psi_vort - A*circ
            end do
            istart = istart+nd
         end do
c
      return
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine  FASMVP (k, nd, nbk, x_zeta, y_zeta, zeta,    
     1                    dzeta, zeta_k, diag, A_k, u, w, charge,   
     2                    grad, pot,hess,source,dipstr,dipvec)
c---------------
c
      implicit real*8 (a-h, o-z)
      integer*4 inform(10), ier, iprec, ifcharge, ifpot,
     1	       ifgrad,ifhess,ifdipole,k,nd,nbk	
      complex*16 zeta(nbk), dzeta(nbk), charge(nbk), grad(2,nbk), 
     1           zQsum, zQ2sum, eye, zeta_k(k), zdis, 
     1           pot(nbk), hess(3,nbk),dipstr(nbk)
      dimension u(*),w(*), x_zeta(nbk), 
     1          y_zeta(nbk), diag(nbk), A_k(k),
     1          source(2,nbk),dipvec(2,nbk)
      REAL*4 TIMEP(2), ETIME
c
         pi = 4.D0*DATAN(1.D0)
         eye = DCMPLX(0.D0,1.D0)
         dalph = 2.d0*pi/nd
c
	do i = 1,nbk
		dipstr(i) = dcmplx(u(i),0.d0)
	end do
c Calculate A_k
         istart = 0
         do kbod = 1, k
            A_k(kbod) = 0.d0
            do i = 1, nd
ccc               A_k(kbod) = A_k(kbod) 
ccc     1               - u(istart+i)*dimag(dzeta(istart+i))
               A_k(kbod) = A_k(kbod) + u(istart+i)
            end do
ccc            A_k(kbod) = A_k(kbod)*dalph/(2.d0*pi)
            A_k(kbod) = A_k(kbod)*dalph
            istart = istart+nd
         end do
         A_k(1) = 0.d0
c
         zQsum = 0.d0
         do i = 1, nbk
            charge(i) = 0.d0
            zQ2sum = dalph*u(i)*dzeta(i)*dconjg(zeta(i))
     1                   /(1.d0+cdabs(zeta(i))**2)
            zQsum = zQsum - dreal(zQ2sum/(2.d0*pi*eye))
         end do
ccc         call PRIn2 (' zQsum = *', zQsum, 2)
c
c         tbeg = etime(timep)
c         iout(1) = 0
c         iout(2) = 13
c         iflag7 = 3
c         napb = 20
c         ninire = 2
c         mex = 300
c         eps7 = 1.d-14
c         tol = 1.d-14
c         nnn = nbk
c         call DAPIF2 (iout, iflag7, nnn, napb, ninire, mex, ierr, 
c     &                inform, tol, eps7, x_zeta, y_zeta, qa, poten,  
c     &                cfield, wksp, nsp, CLOSE)

	
c Set parameters for FMM call
	iprec    = 5
	ifcharge = 0
	ifdipole = 1
	ifpot    = 1
	ifgrad   = 0
	ifhess   = 0
	
	do i = 1,nbk
	   source(1,i) = x_zeta(i)
	   source(2,i) = y_zeta(i)
	   dipvec(1,i) = -dimag(dzeta(i))
	   dipvec(2,i) = dreal(dzeta(i))  		
	end do


	call lfmm2dpartself(ier,iprec,nbk,source,ifcharge, 
     &			charge,ifdipole,dipstr,dipvec,ifpot,
     &			pot,ifgrad,grad,ifhess,hess) 	
         call PRINI (6, 13)
ccc         call PRIN2 (' qa = *', qa, 2*nnn)
ccc         call PRIN2 (' cfielf = *', cfield, 2*nnn)
c         if (ier.ne.0) then
c            write (6,*) '  ERROR IN FMM '
c            stop
c         end if

	  if (ier.eq.4) then
            print *, 'ERROR IN FMM: Cannot allocate tree workspace'
            stop
	   else if(ier.eq.8) then
		print *, 'ERROR IN FMM: Cannot allocate bulk FMM workspace'
		stop
	   else if(ier.eq.16) then
		print *, 'ERROR IN FMM: Cannot allocate multipole expansion workspace 
     1			in FMM' 
		stop
         end if

ccc         call PRINF (' Number of Levels used = *', inform(3), 1)
ccc         call PRIN2 (' cfield = *', cfield, 2*nbk)
c
c Fix up field
         istart = 0
         do kbod = 1, k
            do i = 1, nd
               zQ2sum = dalph*u(istart+i)*dzeta(istart+i)
     1                   *dconjg(zeta(istart+i))
     2                   /(1.d0+cdabs(zeta(istart+i))**2)
               zQ2sum = - dreal(zQ2sum/(2.d0*pi*eye))
               pot(istart+i) = dreal(dalph*pot(istart+i)/(2.d0*pi)) 
     1				- zQsum + zQ2sum
               w(istart+i) = 0.5d0*u(istart+i) + pot(istart+i) 
     1                        - dalph*diag(istart+i)*u(istart+i)
     2                        - A_k(kbod)
            end do
            istart = istart+nd
         end do
c
      return
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine  BUILD_MAT_STEREO (k, nd, nbk, x_zeta, y_zeta,     
     1                              zeta, dzeta, zeta_k, diag, amat)
c---------------
c
      implicit real*8 (a-h, o-z)
      integer*4 iout(2), inform(10), ierr(10)
      complex*16 zeta(nbk), dzeta(nbk),  
     1           zQsum, zQ2sum, eye, zeta_k(k), zdis, zsrc
      dimension x_zeta(nbk), 
     1          y_zeta(nbk), diag(nbk), amat(nbk+k,nbk+k)
      REAL*4 TIMEP(2), ETIME
c
         pi = 4.D0*DATAN(1.D0)
         eye = DCMPLX(0.D0,1.D0)
         dalph = 2.d0*pi/nd
c
ccc         call PRIn2 (' zQsum = *', zQsum, 2)
c
c Fix up field
         do i = 1, nbk
            do j = 1, nbk
               if (i.ne.j) then
                  zdis = zeta(i)-zeta(j)
                  zsrc = dzeta(j)/zdis
                  zsrc = dalph*zsrc/(2.d0*pi*eye)
                  zQsum = dzeta(j)*dconjg(zeta(j))
     1                   /(1.d0+cdabs(zeta(j))**2)
                  zQsum = dalph*zQsum/(2.d0*pi*eye)
                  amat(i,j) = dreal(zsrc+zQsum)
                else
                  amat(i,i) = 0.5d0-dalph*diag(i) 
               end if
            end do
c
c Add on log singularities
            do kbod = 1, k
ccc            do kbod = 1, 1
               zdis = zeta(i) - zeta_k(kbod)
               rad = 2.d0*(cdabs(zdis))**2/((1+(cdabs(zeta(i)))**2)
     1                  *((1+(cdabs(zeta_k(kbod)))**2)))
               amat(i,nbk+kbod)= 0.5d0*dlog(rad)
            end do
         end do
c
c
c Constraints for multiple log but with same-valued streamlines
         do kbod = 1, k
            amat(nbk+1,nbk+kbod) = 1.d0 
         end do
         do kbod = 2, k
            do i = (kbod-1)*nd+1, kbod*nd
ccc               amat(nbk+kbod,i) = cdabs(dzeta(i))
               amat(nbk+kbod,i) = 1.d0
            end do
         end do 
c
c Dump it out
         open (unit = 24, file = 'amat_stereo.dat')
         do i = 1, nbk+k
            do j = 1, nbk+k
               write(24,'(e20.13,$)')(amat(i,j))
               write (24,'(a)')  ''
            end do
         end do
         close (24) 
c
      return
      end
c
c
c---------------
      subroutine SOLVE (nd, k, kmax, nbk, rhs, soln, A_k, rwork, 
     1                  lrwork, iwork, liwork, maxl)
c---------------
c
      implicit real*8 (a-h,o-z)
      external MATVEC_LAPL, MSOLVE
c
c  System
      dimension soln(*), rhs(*), A_k(k)
c
c  DGMRES work arrays
      dimension rwork(lrwork), iwork(liwork)
c
c  Timings
c
      real*4 timep(2), etime
c
         pi = 4.d0*datan(1.d0)
c
c  solve linear system using GMRES.
c
c
c     parameters for DGMRES
         itol = 0
         tol = 1.0d-11
         isym = 0
         iwork(1) = maxl
         do i=2,liwork
            iwork(i) = 0
         enddo
c
c  Preconditioner flag
c
         iwork(4) = 0
c
c  Restart flag
c  
         iwork(5) = 5      
c
c     provide initial guess soln
         norder = nbk
ccc         norder = nbk
         do i=1,norder
            soln(i) = rhs(i)
         enddo
c
         t0 = etime(timep)
         call DGMRES (norder, rhs, soln, nelt, ia, ja, a, isym,
     1               MATVEC_LAPL, MSOLVE, itol, tol, itmax, iter, err,  
     1               ierr, 0, sb, sx, rwork, lrwork, iwork, 
     1               liwork, rw, iw)
         call Prin2 (' after laplace solve, err = *', err, 1)
         call PrinF (' after laplace solve, ierr = *', ierr, 1)
         call PRINI (6,13)
         call PRINF ('  # GMRES ITERATIONS = *',iter,1)
         if (ierr.gt.2) then
            call PRINF ('  SOMETHING WRONG IN GMRES, IERR = *',ierr,1)
            call PRINF ('  iwork = *',iwork,10)
            stop
           elseif (ierr.ge.0) then
            t1 = etime(timep)
            tsec = t1 - t0
            call PRIn2 (' time in solve = *', tsec, 1)
c
c  calculate A_k 
            istart = nd
            dth = 2.d0*pi/nd
            A_k(1) = 0.d0
            do kbod = 2, k
               A_k(kbod) = 0.d0
               do i = 1, nd
                  A_k(kbod) = A_k(kbod) 
     1               + soln(istart+i)
               end do
               A_k(kbod) = A_k(kbod)*dth
               istart = istart+nd
            end do
            call PRIN2 (' A_k = *', A_k, k)
         end if
c
      return
      end
c
c*************************************************
c
c      subroutine RESAMPLE (ns,nsk,k0,kk,nsamp,h,z,dz,xat,yat,u,x,y,
c     *                     wksp,nsp)
c
c  Overresolve the data on the boundary
c
c      implicit real*8 (a-h,o-z)
c      dimension h(k0:kk),ns(k0:kk),xat(*),yat(*),x(*),y(*)
c      complex*16 z(*),dz(*),u(*),wksp(nsp)
c
c         pi = 4.d0*datan(1.d0)
c
c  do z first
c
c         istart = 0
c         istart2 = 0
c         do nbod = k0,kk
c           nd2 = nsamp*ns(nbod)
c            ndm1 = ns(nbod)-1
c            nd2m1 = nd2-1
c            call FTRPIN (wksp,nsp,ip2,ipt,ndm1,nd2m1)
c            call FINTER (xat(istart+1),x(istart2+1),ndm1,nd2m1,wksp,
c     *                   nsp,ip2,ipt)
c            call FINTER (yat(istart+1),y(istart2+1),ndm1,nd2m1,wksp,
c     *                   nsp,ip2,ipt)
c            do i = 1,nsamp*ns(nbod)
c               z(istart2+i) = dcmplx(x(istart2+i),y(istart2+i))
c            end do
c            istart = istart + ns(nbod)
c            istart2 = istart2 + nsamp*ns(nbod)
c         end do
c
c  now do dz
c
c         do i = 1,nsk
c            xat(i) = dreal(dz(i))
c            yat(i) = dimag(dz(i))
c         end do
c         istart = 0
c         istart2 = 0
c         do nbod = k0,kk
c            nd2 = nsamp*ns(nbod)
c            ndm1 = ns(nbod)-1
c           nd2m1 = nd2-1
c            call FTRPIN (wksp,nsp,ip2,ipt,ndm1,nd2m1)
c            call FINTER (xat(istart+1),x(istart2+1),ndm1,nd2m1,wksp,
c     *                   nsp,ip2,ipt)
c            call FINTER (yat(istart+1),y(istart2+1),ndm1,nd2m1,wksp,
c     *                   nsp,ip2,ipt)
c            do i = 1,nsamp*ns(nbod)
c               dz(istart2+i) = dcmplx(x(istart2+i),y(istart2+i))
c            end do
c            istart = istart + ns(nbod)
c            istart2 = istart2 + nsamp*ns(nbod)
c         end do
c
c  now do u
c
c         do i = 1,nsk
c            xat(i) = dreal(u(i))
c            yat(i) = dimag(u(i))
c         end do
c         istart = 0
c         istart2 = 0
c         do nbod = k0,kk
c            nd2 = nsamp*ns(nbod)
c            ndm1 = ns(nbod)-1
c            nd2m1 = nd2-1
c            call FTRPIN (wksp,nsp,ip2,ipt,ndm1,nd2m1)
c            call FINTER (xat(istart+1),x(istart2+1),ndm1,nd2m1,wksp,
c     *                   nsp,ip2,ipt)
c            call FINTER (yat(istart+1),y(istart2+1),ndm1,nd2m1,wksp,
c     *                   nsp,ip2,ipt)
c            do i = 1,nsamp*ns(nbod)
c               u(istart2+i) = dcmplx(x(istart2+i),y(istart2+i))
c           end do
c            istart = istart + ns(nbod)
c            istart2 = istart2 + nsamp*ns(nbod)
c         end do
c
c  Update points and stuff
c
c         nsk = nsamp*nsk
c         do nbod = k0,kk
c            ns(nbod) = nsamp*ns(nbod)
c            h(nbod) = 2.d0*pi/ns(nbod)
c         end do
c         do i = 1,nsk
c            xat(i) = dreal(z(i))
c            yat(i) = dimag(z(i))
c        end do
c
c      return
c      end
c
c
c---------------
      subroutine SOL_GRID (nd, k, nbk, nth, nphi, u, A_k, zeta_k, zeta,  
     1                     dzeta, igrid, zeta_gr, u_gr)
c---------------
c
      implicit real*8 (a-h,o-z)
      dimension u(nbk), igrid(nth,nphi), u_gr(nth,nphi), A_k(k)
      complex*16 zeta(nbk), dzeta(nbk), zeta_gr(nth,nphi), zkern, 
     1           zeta_k(k), zdis, eye
      REAL*4 TIMEP(2), ETIME
c
         pi = 4.d0*datan(1.d0)
         eye = dcmplx(0.d0, 1.d0)
         dalph = 2.d0*pi/nd
c
         tbeg = etime(timep)
         do i = 1, nth
            do j = 1, nphi
               u_gr(i,j) = 0.d0
               if (igrid(i,j).ne.0) then 
               do ip = 1, nbk
                  zdis = zeta(ip) - zeta_gr(i,j)
                  zkern = dzeta(ip)/zdis - dconjg(zeta(ip))*dzeta(ip)
     1                         /(1.d0+cdabs(zeta(ip))**2)
                  u_gr(i,j) = u_gr(i,j) 
     1                   - dalph*u(ip)*dimag(zkern)/(2.d0*pi)
               end do
               do kbod = 1, k
                  zdis = zeta_gr(i,j) - zeta_k(kbod)
                  rad = 
     1              2.d0*(cdabs(zdis))**2/((1+(cdabs(zeta_gr(i,j)))**2)
     2                  *((1+(cdabs(zeta_k(kbod)))**2)))
                  u_gr(i,j) = u_gr(i,j) + A_k(kbod)*0.5d0*dlog(rad)
               end do
               end if
            end do
         end do
         tend = etime(timep)
ccc         call PRIN2 (' poten = *', poten, n)
         call PRIN2 (' TIME FOR GRID = *',tend-tbeg,1)
c
      return
      end
c
c
c---------------
      subroutine SOL_GRID_FMM (nd, k, nbk, nth, nphi, u, zeta_k,   
     1                         zeta, dzeta, igrid, zeta_gr, u_gr,
     2                         qa,grad, pot,gradtarg,pottarg,
     3			       hess,hesstarg,source,targ,dipvec,
     4                         nvort, vort_k, zk_vort, gamma_tot)
c---------------
c
      implicit real*8 (a-h,o-z)
      dimension igrid(nth,nphi), u_gr(nth,nphi), vort_k(nvort)
      complex*16 zeta(nbk), dzeta(nbk), zeta_gr(nth,nphi), zkern, 
     1           zeta_k(k), zdis, eye, qa(*), grad(2,*), zQsum, 
     2           zQ2sum, zk_vort(nvort),hess(3,nbk),dipstr(nbk)
     3		,pot(nbk),pottarg(nth*nphi),gradtarg(2,nth*nphi),
     4		hesstarg(3,nth*nphi)
      dimension u(nbk),source(2,nbk), dipvec(2,nbk), targ(2,nth*nphi)
      integer*4 iout(2), inform(10), ier,iprec,ifcharge,
     1		ifpot,ifgrad,ifhess,ifdipole,ifpottarg,
     2		ifgradtarg,ifhesstarg,ntarg
      REAL*4 TIMEP(2), ETIME
c
         pi = 4.d0*datan(1.d0)
         eye = dcmplx(0.d0, 1.d0)
         dalph = 2.d0*pi/nd
         A = -gamma_tot
c
ccc         call prin2 (' in sol_grid_fmm, vort_k = *', vort_k, k)
ccc         call prin2 ('                  zk_vort = *',zk_vort, 2*k)
ccc         call prin2 ('                  gamma_tot = *', gamma_tot, 1)
ccc         call prin2 ('                 zeta_k = *', zeta_k, 2*k)
c
c pack zeta and zeta_gr into x_zeta and y_zeta
         zQsum = 0.d0
         do i = 1, nbk
            source(1,i) = dreal(zeta(i))
            source(2,i) = dimag(zeta(i))
            qa(i) = 0.d0
            zQ2sum = dalph*u(i)*dzeta(i)*dconjg(zeta(i))
     1                   /(1.d0+cdabs(zeta(i))**2)
            zQsum = zQsum - zQ2sum/(2.d0*pi*eye)
	   dipvec(1,i) = -dimag(dzeta(i))
	   dipvec(2,i) = dreal(dzeta(i))
	   dipstr(i) = dcmplx(u(i),0.d0)
         end do
	ntarg = nphi*nth
         ij = 0
         do i = 1, nth
            do j = 1, nphi
               if (igrid(i,j).ne.0) then
                  ij = ij + 1
                  targ(1,ij) = dreal(zeta_gr(i,j))
                  targ(2,ij) = dimag(zeta_gr(i,j))
               end if
            end do
         end do
         call PRINF (' ij = *', ij, 1)
ccc         call PRIn2 (' zQsum = *', zQsum, 2)
c
c         tbeg = etime(timep)
c         iout(1) = 0
c         iout(2) = 13
c         iflag7 = 3
c         napb = 50
c         ninire = 2
c         mex = 300
c         eps7 = 1.d-14
c         tol = 1.d-14
c         nnn = ij
cccc         nnn = nbk
ccc         nnn = nbk+100
c Set Parameters for FMM
			  	  
	iprec    = 5
	ifcharge = 0
	ifdipole = 1
	ifpot    = 1
	ifgrad   = 0
	ifhess   = 0
	ifpottarg = 1
	ifgradtarg = 0
	ifhesstarg = 0


	call lfmm2dparttarg(ier,iprec,nbk,source,ifcharge, 
     &			qa,ifdipole,dipstr,dipvec,ifpot,
     &			pot,ifgrad,grad,ifhess,hess,ntarg,targ,
     &			ifpottarg,pottarg,ifgradtarg,gradtarg,
     &			ifhesstarg,hesstarg) 	
         call PRINI (6, 13)	
c         call DAPIF2 (iout, iflag7, nnn, napb, ninire, mex, ierr, 
c     &                inform, tol, eps7, xat, yat, qa, poten,  
c     &                cfield, wksp, nsp, CLOSE)
c         call PRINI (6, 13)
ccc         call PRIN2 (' cfield = *', cfield, 2*nnn)
c         call PRINF (' Number of Levels = *', inform(3), 1)
c         if (ierr(1).ne.0) then
c            write (6,*) '  ERROR IN DAPIF2, IERR = ', (ierr(ii),ii=1,6)
c           write(6,*) '  INFORM = ', (inform(ii),ii=1,6)
c            stop
c         end if
ccc         call PRINF (' Number of Levels used = *', inform(3), 1)
ccc         call PRIN2 (' cfield = *', cfield, 2*nbk)
c

	if (ier.eq.4) then
            print *, 'ERROR IN FMM: Cannot allocate tree workspace'
            stop
	   else if(ier.eq.8) then
		print *, 'ERROR IN FMM: Cannot allocate bulk FMM workspace'
		stop
	   else if(ier.eq.16) then
		print *, 'ERROR IN FMM: Cannot allocate multipole expansion workspace 
     1			in FMM' 
		stop
         end if

c Fix up field
         ij = 0
         umax = -1.d10
         umin = 1.d10
         do i = 1, nth
            do j = 1, nphi
               u_gr(i,j) = -10.d0
               if (igrid(i,j).ne.0) then     
                  ij = ij + 1
		  pottarg(ij) = dreal(dalph*pottarg(ij)/(2*pi))           
                  u_gr(i,j) = pottarg(ij) - zQsum
                  psi_vort = 0.d0
                  call POINT_VORTEX (zeta_gr(i,j), zeta_k(1), circ)
                  do ivort = 1, nvort
                     call POINT_VORTEX (zeta_gr(i,j), zk_vort(ivort), 
     1                                  psi)
                     psi_vort = psi_vort + vort_k(ivort)*psi
                  end do
                  u_gr(i,j) = u_gr(i,j) + psi_vort + A*circ
ccc                  u_gr(i,j) = psi_vort + A*circ
                  umax = max(umax,u_gr(i,j))
                  umin = min(umin,u_gr(i,j))
                 else
                  u_gr(i,j) = -1000.d0
               end if
            end do
         end do
         call PRIN2 (' Max solution = *', umax, 1)
         call PRIN2 (' Min solution = *', umin, 1)
c
         tend = etime(timep)
ccc         call PRIN2 (' poten = *', poten, n)
         call PRIN2 (' TIME FOR FMM  ON GRID = *',tend-tbeg,1)
ccc         call PRIN2 (' cfield = *', cfield, 2*n)
         open (unit = 43, file = 'ugrid.dat')
            call DUMP (nth, nphi, u_gr, igrid, 1, 43)
            close(43)
c
      return
      end
c
c
c---------------
      subroutine SOL_TAR_FMM (nd, k, nbk, ntar, u, zeta_k, zeta,  
     1                        dzeta, xz_tar, yz_tar, u_tar, 
     2                        qa,dipstr,grad, pot,gradtarg,pottarg,hess,
     3			      hesstarg,source,targ,dipvec,nvort, 
     4                        vort_k, zk_vort, gamma_tot)
c---------------
c
      implicit real*8 (a-h,o-z)
      dimension u(nbk),u_tar(ntar), vort_k(nvort), xz_tar(ntar),
     1          yz_tar(ntar)
      complex*16 zeta(nbk), dzeta(nbk), zeta_k(k), zdis, eye, qa(nbk),
     1            grad(2,nbk), zQsum, zQ2sum, zk_vort(nvort), ztar,
     2		pot(nbk),hess(3,nbk),pottarg(ntar),dipstr(nbk),
     3		gradtarg(ntar),hesstarg(ntar)
      dimension source(2,nbk),targ(2,ntar),dipvec(2,nbk)
      integer*4 ier,iprec,ifcharge,ifdipole,ifpot,ifgrad,ifhess,
     1          ifpottarg,ifgradtarg,ifhesstarg
      REAL*4 TIMEP(2), ETIME
c
         pi = 4.d0*datan(1.d0)
         eye = dcmplx(0.d0, 1.d0)
         dalph = 2.d0*pi/nd
         A = -gamma_tot
c
c pack zeta and zeta_gr into x_zeta and y_zeta
         zQsum = 0.d0
         do i = 1, nbk
            source(1,i) = dreal(zeta(i))
            source(2,i) = dimag(zeta(i))
            zQ2sum = dalph*u(i)*dzeta(i)*dconjg(zeta(i))
     1                   /(1.d0+cdabs(zeta(i))**2)
            zQsum = zQsum - dreal(zQ2sum/(2.d0*pi*eye))
	   dipvec(1,i) = -dimag(dzeta(i))
	   dipvec(2,i) = dreal(dzeta(i))
	   qa(i) = 0.d0
	   dipstr(i) = dcmplx(u(i),0.d0)
         end do
         
         do i = 1, ntar
            targ(1,i) = xz_tar(i)
            targ(2,i) = yz_tar(i)
         end do
ccc         call PRIn2 (' zQsum = *', zQsum, 2)
c
c         tbeg = etime(timep)
c         iout(1) = 0
c         iout(2) = 13
c         iflag7 = 3
c         napb = 50
c         ninire = 2
c         mex = 300
c         eps7 = 1.d-14
c         tol = 1.d-14
c         nnn = ij
ccc         nnn = nbk
ccc         nnn = nbk+100
c Set parameters for FMM call
	iprec    = 5
	ifcharge = 0
	ifdipole = 1
	ifpot    = 1
	ifgrad   = 0
	ifhess   = 0
	ifpottarg = 1
	ifgradtarg = 0
	ifhesstarg = 0


	call lfmm2dparttarg(ier,iprec,nbk,source,ifcharge, 
     &			qa,ifdipole,dipstr,dipvec,ifpot,
     &			pot,ifgrad,grad,ifhess,hess,ntar,targ,
     &			ifpottarg,pottarg,ifgradtarg,gradtarg,
     &			ifhesstarg,hesstarg) 	


c         call DAPIF2 (iout, iflag7, nnn, napb, ninire, mex, ierr, 
c     &                inform, tol, eps7, xat, yat, qa, poten,  
c     &                cfield, wksp, nsp, CLOSE)
c         call PRINI (6, 13)
ccc         call PRIN2 (' cfield = *', cfield, 2*nnn)
c         call PRINF (' Number of Levels = *', inform(3), 1)
c         if (ierr(1).ne.0) then
c            write (6,*) '  ERROR IN DAPIF2, IERR = ', (ierr(ii),ii=1,6)
c            write(6,*) '  INFORM = ', (inform(ii),ii=1,6)
c            stop
c         end if
ccc         call PRINF (' Number of Levels used = *', inform(3), 1)
ccc         call PRIN2 (' cfield = *', cfield, 2*nbk)
c
ccc         call PRIN2 (' a_k in sol_GRID_FMM = *', A_k, k)

	if (ier.eq.4) then
            print *, 'ERROR IN FMM: Cannot allocate tree workspace'
            stop
	   else if(ier.eq.8) then
		print *, 'ERROR IN FMM: Cannot allocate bulk
     $			 FMM workspace'
		stop
	   else if(ier.eq.16) then
		print *, 'ERROR IN FMM: Cannot allocate multipole
     $  			 expansion workspace 
     1			in FMM' 
		stop
         end if

c Fix up field
         do i = 1, ntar
	   pottarg(i) = (dalph*pottarg(i)/(2*pi))          
            u_tar(i) = pottarg(i) - zQsum
            ztar = dcmplx(xz_tar(i),yz_tar(i))
            call POINT_VORTEX (ztar, zeta_k(1), circ)
            psi_vort = 0.d0
            do ivort = 1, nvort
               call POINT_VORTEX (ztar, zk_vort(ivort), psi)
               psi_vort = psi_vort + vort_k(ivort)*psi
            end do
            u_tar(i) = u_tar(i) + psi_vort + A*circ
         end do
ccc         call prin2 (' u_tar = *', u_tar, ntar)
c
         tend = etime(timep)
ccc         call PRIN2 (' poten = *', poten, n)
         call PRIN2 (' TIME FOR FMM  ON GRID = *',tend-tbeg,1)
ccc         call PRIN2 (' cfield = *', cfield, 2*n)
c
      return
      end
c
c
c---------------
      subroutine SOL_VORT_CHECK (nd, k, nbk, nth, nphi, nvort, zeta_gr,  
     1                           igrid, u_gr, uex_gr, zk_vort, r0)
c---------------
c
      implicit real*8 (a-h,o-z)
      dimension u_gr(nth,nphi), uex_gr(nth,nphi), igrid(nth,nphi)
      complex*16 zeta_gr(nth,nphi), zk_vort(nvort), zet, z, zvort, eye,
     1           zeta_vort, zdis
c
         pi = 4.d0*datan(1.d0)
         eye = dcmplx(0.d0, 1.d0)
c
         err = 0.d0
         call PRIN2 (' zk_vort = *', zk_vort, 2*nvort)
         do i = 1, 1
            do j = 1, 1
ccc               if (igrid(i,j).eq.1) then 
                  z = zeta_gr(i,j)
                  zeta_vort = eye*(r0-zk_vort(1))/(r0+zk_vort(1))
                  zet = eye*(r0-z)/(r0+z)
                  arg = cdabs((zet-zeta_vort)/(zet-dconjg(zeta_vort)))
                  uex_gr(i,j) = -dlog(arg)
                  err = max(err,dabs(uex_gr(i,j)-u_gr(i,j)))
ccc                  call PRIN2 (' uex = *', uex_gr(i,j), 1)
ccc                  call PRIn2 (' u = *', u_gr(i,j), 1)
ccc               end if
            end do
         end do
         call PRIN2 (' error in vorticity cap solution = *', err, 1)
c
      return
      end
c
c
c---------------
      subroutine CALC_VEL (k, nd, nbk, nvort, u, gamma_tot, zeta, dzeta, 
     1                     zeta_k, vort_k, zk_vort, zvel)
c---------------
c Calculate velocity at each of the point vortex locations according
c to formula (13) in Crowdy (2006)
c
      implicit real*8 (a-h,o-z)
      dimension u(nbk), vort_k(nvort)
      complex*16 zeta(nbk), dzeta(nbk), zeta_k(k), zgrad,  
     1           zk_vort(k), zvel(nvort), eye, zsum
c
         pi = 4.d0*datan(1.d0)
         eye = dcmplx(0.d0, 1.d0)
         dth = 2.d0*pi/nd
         A = -gamma_tot
c
         istart = 0
         do ivort = 1, nvort
            zvel(ivort) = 0.d0
            vfactor = 1.d0 + (cdabs(zk_vort(ivort)))**2
            zsum = 0.d0
c
c         add on contributions from integral operator
            do kbod = 1, k
               do i = 1, nd
                  zsum = zsum + u(istart+i)*dzeta(istart+i)/
     1                     (zk_vort(ivort)-zeta(istart+i))**2
               end do
               zsum = 0.25d0*eye*dth*zsum/pi
ccc               zvel(ivort) = zvel(ivort) + 0.5d0*eye*zf(1)/nd
               istart = istart + nd
            end do
            zvel(ivort) = zsum
ccc            call prin2 (' contribution from integral operator = *', 
ccc     1                      zsum, 2)
c
c         add on contributions from other point vortices
            do kvort = 1, nvort
               if (kvort.ne.ivort) then
                  call GRAD_POINT_VORTEX (zk_vort(ivort), 
     1                                    zk_vort(kvort), zgrad)
                  zvel(ivort) = zvel(ivort) + vort_k(kvort)*zgrad
               end if
            end do
c
c         add on contribution from log at north pole        
            call GRAD_POINT_VORTEX (zk_vort(ivort), zeta_k(1), zgrad)
            zvel(ivort) = zvel(ivort) + A*zgrad
ccc            call prin2 (' contribution from log = *', A*zgrad, 2)
c
c  Calculate final velocity according to Crowdy's formula
            zvel(ivort) = dconjg(-eye*0.5*zvel(ivort)*vfactor**2)
         end do
         call PRIN2 (' zvel = *', zvel, 2*nvort)
         call prin2 (' norm of velocity = *', cdabs(zvel(1)),1)
c
      return
      end
c
c
c---------------
      subroutine CHECK_ERROR (nd, k, nbk, nth, nphi, zeta_k, igrid, 
     1                        zeta_gr, u_gr)
c---------------
c
      implicit real*8 (a-h,o-z)
      dimension igrid(nth,nphi), u_gr(nth,nphi)
      complex*16 zeta_gr(nth,nphi), zeta_k(k), eye
c
         pi = 4.d0*datan(1.d0)
         eye = dcmplx(0.d0, 1.d0)
         dalph = 2.d0*pi/nd
c
         call PRIn2 (' zeta_k in check_ERROR = *', zeta_k, 2*k)
         call PRIN2 (' zeta_gr = *', zeta_gr(1,1), 2)
         err = 0.d0
         do i = 1, nth
            do j = 1, nphi
               if (igrid(i,j).eq.1) then
               u_ex = 0.d0
               do mbod = 1, k
                  u_ex = u_ex + 1.d0/(zeta_gr(i,j)-zeta_k(mbod)) 
     1                   + dconjg(1.d0/(zeta_gr(i,j)-zeta_k(mbod)))
               end do
ccc               u_ex = u_ex + 66.d0 
               err = max(err,dabs(u_ex-u_gr(i,j)))
               call PRIN2 ('### u_ex = *', u_ex, 1)
               call PRIN2 ('    u_gr  = *', u_gr(i,j), 1)
               end if
            end do
         end do
         call PRIN2 (' max error in solution on grid = *', err, 1)
c
c
      return
      end
c
c
c---------------
      subroutine CHECK_ERROR_TAR (nd, k, nbk, ntar, zeta_k, zeta_tar,  
     1                            u_tar, nvort, vort_k, zk_vort, r0)
c---------------
c
      implicit real*8 (a-h,o-z)
      dimension u_tar(ntar), vort_k(nvort)
      complex*16 zeta_tar(ntar), zeta_k(k), eye,
     1           zk_vort(nvort), z, zeta_vort, zet
c
         pi = 4.d0*datan(1.d0)
         eye = dcmplx(0.d0, 1.d0)
         dalph = 2.d0*pi/nd
c
         call PRIn2 (' zeta_k in check_ERROR = *', zeta_k, 2*k)
         call PRIN2 (' zeta_tar = *', zeta_tar, 2)
ccc         zeta_k(1) = dcmplx(-10.d0,-10.d0)
         call prin2 (' zeta_k in check error = *', zeta_k, 2)
         call prin2 (' zeta_tar = *', zeta_tar, 2*ntar)
         call PRINF (' In CHECK_ERROR_TAR, NTAR = *', ntar, 1)
         umax = 0.d0
         err = 0.d0
         do i = 1, ntar
            u_ex = 0.d0
            z = zeta_tar(i)
            zeta_vort = eye*(r0-zk_vort(1))/(r0+zk_vort(1))
            zet = eye*(r0-z)/(r0+z)
            arg = cdabs((zet-zeta_vort)/(zet-dconjg(zeta_vort)))
            u_ex = -vort_k(1)*dlog(arg)/(2.d0*pi)
            err = max(err,dabs(u_ex-u_tar(i)))
            umax = max(umax,u_ex)
            call PRIN2 ('### u_ex = *', u_ex, 1)
            call PRIN2 ('    u_tar  = *', u_tar(i), 1)
         end do
         call PRIN2 (' max abs error in solution = *', err, 1)
         call PRIN2 (' max rel error in solution = *', err/umax, 1)
c
c
      return
      end
c
C********************************************************************
      SUBROUTINE MSOLVE(N, R, Z, NELT, IA, JA, A, ISYM, RWORK, IWORK)
C********************************************************************
c
c     Another routine required by DGMRES. It allows the use of a
c     preconditioner.
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c
      DIMENSION R(N), Z(N)
c
      RETURN
      END
c
c
c---------------
      subroutine MATVEC_LAPL (N, XX, YY, NELT, IA, JA, A, ISYM)
c---------------
c
c  Required by DGMRES with this precise calling sequence.
c  We ignore most of the parameters except N, XX and YY
c  and call the fast multipole routine FASMVP to do the actual
c  work.
c
      implicit double precision (a-h,o-z)
      dimension xx(n), yy(n)
      parameter (kmax = 10, npmax = 2048, nmax = kmax*npmax)
      parameter (nth_max = 1000, nphi_max = 1000, 
     1          ng_max = nth_max*nphi_max)
      parameter (nsp = 20*nmax + 20*ng_max)
c      
      common /geometry/ x_zeta, y_zeta, zeta, dzeta
      common /inteqn/ diag, zeta_k
      common /sys_size/ k, nd, nbk
c
      dimension x_zeta(nmax), y_zeta(nmax),   
     1          diag(nmax)
      complex*16 zeta(nmax), dzeta(nmax), zeta_k(kmax)
c
c Fast Multipole Arrays
c
      complex*16 qa(nmax), cfield(2,nmax), poten(nmax),hess(3*nmax),
     1		 dipstr(nmax)
      dimension source(2*nmax),dipvec(2*nmax)
c
c  local work arrays
c
      real*4 timep(2), etime
      dimension A_k(kmax)
      complex*16 eye
c
         eye = dcmplx(0.d0,1.d0)
c
         t0 = etime(timep)
c
         call PRINI (6,13)
ccc         zeta_k(1) = dcmplx(-10.d0,-10.d0)
c
ccc         call MATVEC_SPHERE (k, nd, nbk, xs, ys, zs, xn, yn, zn,
ccc     1                       dsda, diag, cx, cy, cz, A_k, xx, yy)
         call FASMVP (k, nd, nbk, x_zeta, y_zeta, zeta, dzeta,   
     1                zeta_k, diag, A_k, xx, yy, qa,   
     2                cfield, poten,hess,source,dipstr,dipvec)
ccc         call FASMVP_TEST (k, nd, nbk, nsp, x_zeta, y_zeta, zeta,    
ccc     1                dzeta, zeta_k, diag, dsda, A_k, xx, yy, qa,   
ccc     2                cfield, poten, wksp)
ccc         call FASMVP (nd, nbk, k, x, y, dz, dsdth, rkappa, 
ccc     *                zk, xx, yy, h, za, pk, stress, dens)
ccc         call MSOLVE(N, yy, xx, NELT, IA, JA, A, ISYM, RWORK, IWORK)
ccc         call PRIn2 ('  xx = *', xx, n)
         t1 = etime(timep)
         tsec = t1 - t0
ccc         WRITE(13,*) 'TIME IN SECONDS FOR MATVEC = ',TSEC
ccc         WRITE(6,*) 'TIME IN SECONDS FOR MATVEC = ',TSEC
c
      RETURN
      END
c
c*********************
c
      subroutine DUMP (nx,ny,ugrid,igrid,ireal,if)
      implicit real*8 (a-h,o-z)
      dimension ugrid(nx,ny), igrid(nx,ny)
c
         DO i = 1,NX
            do j = 1, NY
               if (ireal.eq.1) then 
                  write(if,'(e20.13,$)')(ugrid(I,J))
                  write (if,'(a)')  ''
                 else
                  write(if,'(i4,$)') (igrid(i,j))
                  write (if,'(a)')  ''
               end if
            end do
         ENDDO
c
c periodic
            do j = 1, NY
               if (ireal.eq.1) then 
                  write(if,'(e20.13,$)')(ugrid(1,J))
                  write (if,'(a)')  ''
                 else
                  write(if,'(i4,$)') (igrid(1,j))
                  write (if,'(a)')  ''
               end if
            end do
c
      return
      end
c
c*********************
c
      subroutine DUMP_MOVIE (nx, ny, time, ugrid, iframe, if)
      implicit real*8 (a-h,o-z)
      dimension ugrid(nx,ny), igrid(nx,ny)
c
         write (if,*) 'NX = ',ny,';'
         write (if,*) 'NY = ',nx+1,';'
         write (if,*) 'a = zeros(NX,NY);'
         WRITE(if,*) 'sol = ['
         DO i = 1,NX
            do j = 1, NY
               write(if,'(e20.13,$)')(ugrid(I,J))
               write (if,'(a)')  ''
            end do
         ENDDO
c
c periodic
         do j = 1, NY
               write(if,'(e20.13,$)')(ugrid(1,J))
               write (if,'(a)')  ''
         end do
         write (if,*) '];'
         write (if,*) 'a(:) = sol(:);'
         write (if,1200) time
         write (if,*) 'figure(1); clf;'
         write (if,*) '   surf(xgrid,ygrid,zgrid,a)'
         write (if,*) '   view([64,-4])'
         write (if,*) '   shading interp'
         write (if,*) '   colormap(jet2)'
         write (if,*) '   lighting phong'
         write (if,*) '   material dull'
         write (if,*) '   camlight(''headlight'')'
         write (if,*) '   caxis([-3 5])'
         write (if,*) '   hold on'
c         write (if,*) '   geo_3d'
         write (if,*) '   axis([-1 1 -1 1 -1 1])'
         write (if,*) '   axis equal'
         write (if,*) '   axis off'
         write (if,*) '   title([ ''t = '', num2str(time)]);'
         write (if,*) 
     1          "fname = sprintf('pngfiles/frame%.4d',", iframe,");"
         write (if,*) "print('-dpng', '-r150', fname);"
c
1200  FORMAT ('time = ',F6.2)
c
c
      return
      end
c
c*********************
c
      subroutine DUMP_MOVIE_ALL (nx, ny, time, ugrid, iframe, if)
      implicit real*8 (a-h,o-z)
      dimension ugrid(nx,ny), igrid(nx,ny)
c
         write (if,*) 'NX = ',ny,';'
         write (if,*) 'NY = ',nx+1,';'
         write (if,*) 'a = zeros(NX,NY);'
         WRITE(if,*) 'sol = ['
         DO i = 1,NX
            do j = 1, NY
               write(if,'(e20.13,$)')(ugrid(I,J))
               write (if,'(a)')  ''
            end do
         ENDDO
c
c periodic
         do j = 1, NY
               write(if,'(e20.13,$)')(ugrid(1,J))
               write (if,'(a)')  ''
         end do
         write (if,*) '];'
         write (if,*) 'a(:) = sol(:);'
         write (if,1200) time
         write (if,*) 'figure(1); clf;'
         write (if,*) 'subplot(2,2,1)'
         write (if,*) '   surf(xgrid,ygrid,zgrid,a)'
         write (if,*) '   view([64,-4])'
         write (if,*) '   shading interp'
         write (if,*) '   colormap(jet2)'
         write (if,*) '   lighting phong'
         write (if,*) '   material dull'
         write (if,*) '   camlight(''headlight'')'
         write (if,*) '   caxis([-3 5])'
         write (if,*) '   hold on'
c         write (if,*) '   geo_3d'
         write (if,*) '   axis([-1 1 -1 1 -1 1])'
         write (if,*) '   axis equal'
         write (if,*) '   axis off'
         write (if,*) '   title([ ''t = '', num2str(time)]);'
         write (if,*) 'subplot(2,2,2)'
         write (if,*) '   surf(xgrid,ygrid,zgrid,a)'
         write (if,*) '   view([-116,-4])'
         write (if,*) '   shading interp'
         write (if,*) '   colormap(jet2)'
         write (if,*) '   lighting phong'
         write (if,*) '   material dull'
         write (if,*) '   camlight(''headlight'')'
         write (if,*) '   caxis([-3 5])'
         write (if,*) '   hold on'
c         write (if,*) '   geo_3d'
         write (if,*) '   axis([-1 1 -1 1 -1 1])'
         write (if,*) '   axis equal'
         write (if,*) '   axis off'
c         write (if,*) '   title([ ''t = '', num2str(time)]);'
         write (if,*) 'subplot(2,2,3)'
         write (if,*) '   v = [0:.1:6.2];'
         write (if,*) '   contour(xzeta_grid,yzeta_grid,a,v)'
         write (if,*) '   hold on'
         write (if,*) '   blob'
         write (if,*) '   axis([-2.5 2.5 -2.5 2.5])'
         write (if,*) '   axis equal'
         write (if,*) '   axis off'
c         write (if,*) '   title([ ''t = '', num2str(time)]);'
         write (if,*) 
     1          "fname = sprintf('pngfiles/frame%.4d',", iframe,");"
         write (if,*) "print('-dpng', '-r150', fname);"
c
1200  FORMAT ('time = ',F6.2)
c
c
      return
      end
c
c*********************
c
      subroutine DUMP_MOVIE_VORT (nx, ny, time, xi, ugrid, iframe, if)
      implicit real*8 (a-h,o-z)
      dimension ugrid(nx,ny), igrid(nx,ny)
      complex*16 xi, eye
c
         eye = dcmplx(0.d0,1.d0)
c
c compute view point
         factor = 1.d0+(cdabs(xi))**2
         x1 = (xi+dconjg(xi))/factor
         x2 = (xi-dconjg(xi))/(eye*factor)
         x3 = -(1.d0-(cdabs(xi))**2)/factor
c
         write (if,*) 'NX = ',ny,';'
         write (if,*) 'NY = ',nx+1,';'
         write (if,*) 'a = zeros(NX,NY);'
         WRITE(if,*) 'sol = ['
         DO i = 1,NX
            do j = 1, NY
               write(if,'(e20.13,$)')(ugrid(I,J))
               write (if,'(a)')  ''
            end do
         ENDDO
c
c periodic
         do j = 1, NY
               write(if,'(e20.13,$)')(ugrid(1,J))
               write (if,'(a)')  ''
         end do
         write (if,*) '];'
         write (if,*) 'a(:) = sol(:);'
         write (if,1200) time
         write (if,*) 'figure(1); clf;'
         write (if,*) 'subplot(1,2,1)'
         write (if,*) '   surf(xgrid,ygrid,zgrid,a)'
         write (if,*) '   shading interp'
         write (if,*) '   colormap(jet2)'
         write (if,*) '   caxis([-3 5])'
         write (if,*) '   hold on'
c         write (if,*) '   geo_3d'
         write (if,*) '   axis([-1 1 -1 1 -1 1])'
         write (if,*) '   axis equal'
         write (if,*) '   view([',x1,',',x2,',',x3,'])'
         write (if,*) '   lighting phong'
         write (if,*) '   material dull'
         write (if,*) '   camlight(''headlight'')'
         write (if,*) '   axis off'
         write (if,*) '   title([ ''t = '', num2str(time)]);'
         write (if,*) 'subplot(1,2,2)'
         write (if,*) '   v = [0:.1:6.2];'
         write (if,*) '   contour(xzeta_grid,yzeta_grid,a,v)'
         write (if,*) '   hold on'
         write (if,*) '   blob'
         write (if,*) '   axis([-2.5 2.5 -2.5 2.5])'
         write (if,*) '   axis equal'
         write (if,*) '   axis off'
c         write (if,*) '   title([ ''t = '', num2str(time)]);'
         write (if,*) 
     1          "fname = sprintf('pngfiles/frame%.4d',", iframe,");"
         write (if,*) "print('-dpng', '-r150', fname);"
c
1200  FORMAT ('time = ',F6.2)
c
c
      return
      end
c
C**********************************************************************
C
	SUBROUTINE CLOSE(EPS,I1,I2,Q1,Q2,X1,Y1,X2,Y2,
     *                   FI1,FI2,POT1,POT2)
C
C   *** DESCRIPTION :
C
C       Close approch subroutine provided by the user.
C   Called when two particles are closer to each other than EPS.
C
C   *** INPUT PARAMETERS :
C
C   EPS     =  close approach distance
C   I1      =  number of first  particle
C   I2      =  number of second particle
C   Q1      =  charge of first  particle
C   Q2      =  charge of second particle
C   X1      =  x-component of position of first  particle
C   Y1      =  y-component of position of first  particle
C   X2      =  x-component of position of second particle
C   Y2      =  y-component of position of second particle
C
C   *** OUTPUT PARAMETERS :
C
C   FI1     =  field contribution from first particle on second
C   FI2     =  field contribution from second particle on first
C   POT1    =  potential contribution from first particle on second
C   POT2    =  potential contribution from second particle on first
C
C
C**********************************************************************
C
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	INTEGER II(2)
	DOUBLE PRECISION POS1(2),POS2(2)
	DOUBLE COMPLEX FI1,FI2,OUT2,Q1,Q2
	DATA DONE/1.0/,DHALF/0.5/
	II(1) = I1
	II(2) = I2
	POS1(1) = X1
	POS1(2) = Y1
	POS2(1) = X2
	POS2(2) = Y2
	FI1 = 0
	FI2 = 0
CCC     POT1 = 0
CCC     POT2 = 0
	POT1=DLOG( (X2-X1)**2 + (Y2-Y1)**2 ) /2
	POT2=POT1
ccc	   CALL PRINF('CLOSE PARTICLES ARE:*',II,2)
ccc	   CALL PRIN2(' COORDINATES OF FIRST  PART=*',POS1,2)
ccc	   CALL PRIN2(' COORDINATES OF SECOND PART=*',POS2,2)
	RETURN
	END
