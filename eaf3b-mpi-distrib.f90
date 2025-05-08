! generated from eaf2b-mpi.f90 21th Feb 2023. (V2)
! generated from marmara2b-mpi.f90 on 19th Feb 2023.
! PROGRAM marmara1-mpi
! SIMULATION ON MAIN FAULT 
! created on 18 dec 2002 based upon izmit6-2-mpi.f90
! use friction2-2.f
!
IMPLICIT NONE
INCLUDE 'mpif.h'
INTEGER, PARAMETER :: itmx = 1200, master = 0
INTEGER :: ixmx, iymx, &
	ierr, nrank, nsize, mpmx, mpmy, ndata2, ndata3, icheck
REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: vel, pon
REAL(8), DIMENSION(:, :), ALLOCATABLE :: tau0, nor0, tp, a, dc, &
     		stress, normal,  sigma, w, tr, tau, dtau0, dw0, &
		c1, c3, csh, ctn, dnu1, dnu3, excess, &
        	vmpi, w_mpi, str_mpi, sig_mpi, tau_mpi, nor_mpi
INTEGER, DIMENSION(:, :), ALLOCATABLE :: iv, ivmpi
REAL(8), DIMENSION(:), ALLOCATABLE :: x, y, z, nu1, nu3,  &
		xdum, zdum, nu1dum, nu3dum
INTEGER, DIMENSION(:), ALLOCATABLE :: ix, nconnect
REAL(8) :: pi, mu, const, factor, ds, dt, vp, &
		p000, piece1, piece2, piece3, piece4, &
		ker33s, ker11s, ker31s, sh31, tn33, tn11, dtau, dnor, &
		chi1, chi3, c2, t, dep, dep2, depmax, cohesive, &
		fs, fd, arg1, th1, th, sigma1, sigma3, phi0, rad, dsunit, &
		dtratio, tp0, tr0, dc0, p2, radmin, travel, terminal, &
		dtp, dtr, da, dsigma, wtry, &
  		xwest, xeast, xhypo, yhypo, zhypo, xeuler, yeuler, &
		xlon, ylat, xmax, xmin, x0, y0, r0
INTEGER :: ndir, i, j, k, l, m, n, num, ndata, il, jm, jm2, kn, & 
		ixmain, ib, nmax, ihypo, jhypo, istatic, &
		imp, iproc, mpist(MPI_STATUS_SIZE), ncell, imp2i
EXTERNAL 	ker33s, ker11s, ker31s, imp2i
CHARACTER*50	infile1, infile2
CHARACTER*50 	dir, name0, name3, name6 
CHARACTER  	number*5
!
! INITIALIZATION FOR MPI
!
call MPI_Init( ierr )
call MPI_Comm_rank( MPI_COMM_WORLD, nrank, ierr )
call MPI_Comm_size( MPI_COMM_WORLD, nsize, ierr )
write(*, *) "I am alive! I'm ", nrank+1, " of ", nsize

!
! FIXED PARAMETER
!
pi = acos(-1.0d0)
mu = 32.40d0
vp = 6.00d0
factor = 2.0d0
const = sqrt(3.0d0)/(4.0d0*pi)*mu
terminal = 1.0d0/(2.0d0*sqrt(3.0d0))
p000 = ker31s(0.0d0, 0.0d0, 0.0d0, 0.0d0, factor)

!
! OUTPUT DIRECTORY
!
dir = 'zone1_X20_Opt0_T0.80'
ndir = index(dir, ' ')-1

!
! MODEL PARAMETERS
!
infile2 = dir(1:ndir)//'.prm'

! from input file
!fs = 0.60d0
!fd = 0.800d0*fs
!cohesive = 5.0d0
!dtratio = 0.750d0
dep2 = 12.d0

if (nrank.eq.master) write(*,*) infile2
!
! MPI: READING DATA
! FAULT GEOMETRY for EAF
!
IF( nrank.eq.master ) THEN
  open(21, file=infile2, status='old', err=23)
  read(21, '(a50)') infile1
  read(21,*) fs, fd, cohesive, dtratio
  write(6,*) infile1
  read(21, *) depmax, arg1
    arg1 = arg1*pi/180.
  read(21, *) xwest, xeast
  read(21, *) xhypo, yhypo, zhypo
  close(21)
  xmin = xwest
  xmax = xeast
  x0 = xhypo
  y0 = yhypo
  
  open(11, file=infile1, status='old')
  read(11,*) ixmx
  read(11,*) ixmain, ds
  xmin = xmin/ds
  xmax = xmax/ds
  write(6,*) xwest, xmin
  write(6,*) xeast, xmax

  ALLOCATE ( xdum(ixmx), zdum(ixmx), nu1dum(ixmx), nu3dum(ixmx) )
  l = 0
  do i = 1, ixmain
    read(11,*,err=23) num, xdum(i), zdum(i), nu1dum(i), nu3dum(i) 
    if(xdum(i).ge.xmin.and.xdum(i).le.xmax) then
      l = l + 1
    endif
    if(nrank.eq.master) write(*,*) i, num, xdum(i), zdum(i), nu1dum(i), nu3dum(i)
  enddo
  close(11)

  ixmx = l
  iymx = int(depmax/ds)
  write(6,*) ixmx, iymx
  dt = ds/factor
  x0 = x0/ds
  y0 = y0/ds
  ALLOCATE ( x(ixmx), y(iymx), z(ixmx), nu1(ixmx), nu3(ixmx), ix(ixmx) )
  l = 0
  radmin = 0.
  do i = 1, ixmain
    if(xdum(i).ge.xmin.and.xdum(i).le.xmax) then
      l = l + 1
      x(l) = xdum(i)
      z(l) = zdum(i)
      nu1(l) = nu1dum(i)
      nu3(l) = nu3dum(i)
      ix(l) = 0
      rad = sqrt((x(l)-x0)**2 + (z(l)-y0)**2)
      if(l.eq.1.or.rad.lt.radmin) then
	radmin = rad
	ihypo = l	 
      endif
    endif
  enddo
  DEALLOCATE ( xdum, zdum, nu1dum, nu3dum )

  jhypo = int(zhypo/ds)
  do j = 1, iymx
    y(j) = 0.50d0 + dble(j-1)
  enddo
  write(6,*) ihypo, jhypo

  call MPI_BCAST( ixmx,  1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr )
  call MPI_BCAST( iymx,  1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr )
  call MPI_BCAST( ds,    1, MPI_REAL8, master, MPI_COMM_WORLD, ierr )
  call MPI_BCAST( dt,    1, MPI_REAL8, master, MPI_COMM_WORLD, ierr )
  call MPI_BCAST( ihypo, 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr )
  call MPI_BCAST( jhypo, 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr )

  call MPI_BCAST( x,   ixmx, MPI_REAL8, master, MPI_COMM_WORLD, ierr )
  call MPI_BCAST( y,   iymx, MPI_REAL8, master, MPI_COMM_WORLD, ierr ) 
  call MPI_BCAST( z,   ixmx, MPI_REAL8, master, MPI_COMM_WORLD, ierr )
  call MPI_BCAST( nu1, ixmx, MPI_REAL8, master, MPI_COMM_WORLD, ierr )
  call MPI_BCAST( nu3, ixmx, MPI_REAL8, master, MPI_COMM_WORLD, ierr )
ELSE
  call MPI_BCAST( ixmx,  1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr )
  call MPI_BCAST( iymx,  1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr )
  call MPI_BCAST( ds,    1, MPI_REAL8, master, MPI_COMM_WORLD, ierr )
  call MPI_BCAST( dt,    1, MPI_REAL8, master, MPI_COMM_WORLD, ierr )
  call MPI_BCAST( ihypo, 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr )
  call MPI_BCAST( jhypo, 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr )

  ALLOCATE ( x(ixmx), y(iymx), z(ixmx), nu1(ixmx), nu3(ixmx) )
  call MPI_BCAST( x,   ixmx, MPI_REAL8, master, MPI_COMM_WORLD, ierr )
  call MPI_BCAST( y,   iymx, MPI_REAL8, master, MPI_COMM_WORLD, ierr ) 
  call MPI_BCAST( z,   ixmx, MPI_REAL8, master, MPI_COMM_WORLD, ierr )
  call MPI_BCAST( nu1, ixmx, MPI_REAL8, master, MPI_COMM_WORLD, ierr )
  call MPI_BCAST( nu3, ixmx, MPI_REAL8, master, MPI_COMM_WORLD, ierr )
ENDIF
!
! FAULT PROPERTY AND INITIAL CONDITION
!
ndata2 = ixmx*iymx
ndata3 = ndata2*(itmx + 1)

ALLOCATE (   tau0(ixmx, iymx), nor0(ixmx, iymx),    tp(ixmx, iymx), &
	        a(ixmx, iymx),   dc(ixmx, iymx), sigma(ixmx, iymx), &
		w(ixmx, iymx),   iv(ixmx, iymx),   dw0(ixmx, iymx), &
	       tr(ixmx, iymx),  tau(ixmx, iymx), dtau0(ixmx, iymx), &
           excess(ixmx, iymx) )

IF( nrank.eq.master ) THEN
  do i = 1, ixmx
! optimal plane azimuth from file
!    arg1 = 0.0d0
    th1 = atan2(-nu1(i), nu3(i))
     
    write(6,'(i4, 3f8.3)') i, xlon, arg1*180./pi, th1*180./pi
    do j = 1, iymx
      dep = ds*y(j)
      call coulomb(dep, dep2, fs, fd, dtratio, cohesive, &
    		 sigma1, sigma3, tp0, tr0, phi0, p2)
      call weakening(dep, dc0)
!      write(6,'(i4, 6f11.4)') j, dep, sigma1, sigma3, tp0, tr0, phi0

      th = phi0 - (th1 - arg1)
      w(i,j)  = 0.0d0
      dtau0(i,j) = 0.0d0
      dw0(i,j) = 0.0d0
      iv(i,j) = 0
      tp(i,j) = 5.0d0
      tr(i,j) = 0.0d0
      tau0(i,j) = tp(i,j)*dtratio
      nor0(i,j) = 0.0d0
      dc(i,j) = dc0
      tau0(i,j) = 0.5d0*(sigma1-sigma3)*sin(2.d0*th)
! normal stress satulates beneath dep = dep2
      if( dep.ge.dep2 ) then 
        nor0(i,j) = p2 - 0.5*(sigma1-sigma3)*cos(2.d0*th)
      else
        nor0(i,j) = 0.5d0*(sigma1+sigma3) &
                - 0.5*(sigma1-sigma3)*cos(2.d0*th)
      endif 
      tp(i,j) = cohesive + fs*nor0(i,j)
      tr(i,j) = fd*nor0(i,j)
      excess(i,j) = tp(i,j) - tau0(i,j)
      a(i,j) = ( tp(i,j) - tr(i,j) )/ dc0
      sigma(i,j) = tp(i,j)
	  if (j.eq.jhypo) write(6, '(2i5, 9f10.3)') &
	    i, j, x(i), dep, th1*180./pi, th*180./pi, phi0*180./pi, &
	    tp(i,j), tr(i,j), tau0(i,j), nor0(i,j)  
    enddo
  enddo
ENDIF

call MPI_BCAST( tau0, ndata2, MPI_REAL8, master, MPI_COMM_WORLD, ierr )
call MPI_BCAST( nor0, ndata2, MPI_REAL8, master, MPI_COMM_WORLD, ierr )
call MPI_BCAST(   tp, ndata2, MPI_REAL8, master, MPI_COMM_WORLD, ierr )
call MPI_BCAST(   tr, ndata2, MPI_REAL8, master, MPI_COMM_WORLD, ierr )
call MPI_BCAST(    a, ndata2, MPI_REAL8, master, MPI_COMM_WORLD, ierr )
call MPI_BCAST(   dc, ndata2, MPI_REAL8, master, MPI_COMM_WORLD, ierr )
call MPI_BCAST( sigma, ndata2, MPI_REAL8, master, MPI_COMM_WORLD, ierr )
call MPI_BCAST(    w, ndata2, MPI_REAL8, master, MPI_COMM_WORLD, ierr )
call MPI_BCAST(   iv, ndata2, MPI_INTEGER, master, MPI_COMM_WORLD, ierr )
call MPI_BCAST( dtau0, ndata2, MPI_REAL8, master, MPI_COMM_WORLD, ierr )
call MPI_BCAST(   dw0, ndata2, MPI_REAL8, master, MPI_COMM_WORLD, ierr )

!
! OUTPUT OF MODEL PARAMETERS
!
IF( nrank.eq.master ) THEN
  istatic = 1
  r0 = 3.00d0
!  name6 = dir(1:ndir)//'/hoge4.dat'
!  open(16, file=name6)
!  write(16, '(i5, f9.3)') istatic, r0
!  write(16, '(2i4, 3f11.5)') ihypo, jhypo, x(ihypo), y(jhypo), z(ihypo)

!
! INITIAL CRACK
!
  do i=1, ixmx
    do j=1, iymx
      rad = ds*sqrt((x(i)-x(ihypo))**2 + (y(j)-y(jhypo))**2 + &
	 	(z(i)-z(ihypo))**2)
      if ( rad.le.r0 ) then
	iv(i,j) = 2
        tp(i,j) =  tr(i,j)
	a(i,j) = (tp(i,j)-tr(i,j))/dc0
	sigma(i,j) = tp(i,j)
!	write(16,'(2i5, 2f8.3)') i, j, tp(i,j), tau0(i,j)
      endif
    enddo
  enddo
!  close(16)

  name0 = dir(1:ndir)//'/struct.dat'
  open(10, file=name0)
  ndata = 1
  do i=1, ixmx
    do j=1, iymx
      write(10, '(3i6, 6f10.3)') ndata, i, j, x(i), y(j), z(i), &
	tau0(i,j), sigma(i,j), nor0(i,j)
      ndata = ndata + 1
    enddo
  enddo
  close(10)

ENDIF

call MPI_BCAST( tau0, ndata2, MPI_REAL8, master, MPI_COMM_WORLD, ierr )
call MPI_BCAST( nor0, ndata2, MPI_REAL8, master, MPI_COMM_WORLD, ierr )
call MPI_BCAST(   tp, ndata2, MPI_REAL8, master, MPI_COMM_WORLD, ierr )
call MPI_BCAST(   tr, ndata2, MPI_REAL8, master, MPI_COMM_WORLD, ierr )
call MPI_BCAST(    a, ndata2, MPI_REAL8, master, MPI_COMM_WORLD, ierr )
call MPI_BCAST(   dc, ndata2, MPI_REAL8, master, MPI_COMM_WORLD, ierr )
call MPI_BCAST( sigma, ndata2, MPI_REAL8, master, MPI_COMM_WORLD, ierr )
call MPI_BCAST(    w, ndata2, MPI_REAL8, master, MPI_COMM_WORLD, ierr )
call MPI_BCAST(   iv, ndata2, MPI_INTEGER, master, MPI_COMM_WORLD, ierr )
call MPI_BCAST( dtau0, ndata2, MPI_REAL8, master, MPI_COMM_WORLD, ierr )
call MPI_BCAST(   dw0, ndata2, MPI_REAL8, master, MPI_COMM_WORLD, ierr )

!
! INITIALIZATION
!
ALLOCATE ( vel(ixmx, iymx, 0:itmx), &
		stress(ixmx, iymx), &
		normal(ixmx, iymx) ) 
call MPI_BCAST( tau0, ndata2, MPI_REAL8, master, MPI_COMM_WORLD, ierr )

mpmx = ixmx*iymx/nsize
if( mod((ixmx*iymx), nsize).ne.0 ) mpmx = mpmx + 1
mpmy = iymx

ALLOCATE ( pon(-ixmx:ixmx, -2*iymx:2*iymx, itmx) )
ALLOCATE ( c1(mpmx, ixmx),  c3(mpmx, ixmx), &
	 dnu1(mpmx, ixmx),dnu3(mpmx, ixmx), &
	  csh(mpmx, ixmx), ctn(mpmx, ixmx) )
ALLOCATE ( vmpi(mpmx, mpmy),   w_mpi(mpmx, mpmy), & 
	str_mpi(mpmx, mpmy), nor_mpi(mpmx, mpmy), sig_mpi(mpmx, mpmy), &
	tau_mpi(mpmx, mpmy),   ivmpi(mpmx, mpmy) )

do imp = 1, mpmx
!  i = (imp - 1)*nsize + nrank + 1
!  if( i.gt.ixmx ) exit
  ncell = imp2i(nsize, nrank, mpmx, imp)
  if( ncell.gt.ndata2 ) exit
  i = 1 + (ncell-1)/mpmy
  j = mod(ncell-1, mpmy) + 1
  if( j.eq.1) then
    do l=1, ixmx
      dnu1(imp, l) = -nu1(l)*nu3(i) + nu3(l)*nu1(i)
      dnu3(imp, l) =  nu3(l)*nu3(i) + nu1(l)*nu1(i)
      c1(imp,l) = nu3(l)*(x(i)-x(l)) - nu1(l)*(z(i)-z(l))
      c3(imp,l) = nu1(l)*(x(i)-x(l)) + nu3(l)*(z(i)-z(l))
      csh(imp,l) = dnu3(imp, l)**2 - dnu1(imp, l)**2
      ctn(imp,l) = dnu1(imp, l)*dnu3(imp, l)
    enddo
  endif
enddo
!	
! ITERATION
!
do k = 1, itmx
  call MPI_BCAST( vel,   ndata3, MPI_REAL8, master, MPI_COMM_WORLD, ierr )
  call MPI_BCAST( iv,    ndata2, MPI_INTEGER, master, MPI_COMM_WORLD, ierr )
  call MPI_BCAST( w,     ndata2, MPI_REAL8, master, MPI_COMM_WORLD, ierr )
  call MPI_BCAST( sigma, ndata2, MPI_REAL8, master, MPI_COMM_WORLD, ierr )
  call MPI_BCAST( dtau0, ndata2, MPI_REAL8, master, MPI_COMM_WORLD, ierr )
  call MPI_BCAST( dw0,   ndata2, MPI_REAL8, master, MPI_COMM_WORLD, ierr )

  icheck = 1
  if (k.ne.1) then
    icheck = 0
    do i=1, ixmx
      do j=1, iymx
        if (vel(i,j,k-1).ne.0.0) icheck = 1
      enddo
    enddo
  endif
  if (icheck.eq.0) exit

  do imp = 1, mpmx
!    i = (imp - 1)*nsize + nrank + 1
!    if( i.gt.ixmx ) exit
    ncell = imp2i(nsize, nrank, mpmx, imp)
    if( ncell.gt.ndata2 ) exit
    i = 1 + (ncell-1)/mpmy
    j = mod(ncell-1, mpmy) + 1

    dtau = 0.0d0
    dnor = 0.0d0
    if( k.ne.1 ) then
          do l = 1, ixmx
  	    il = i - l
	    chi1 = c1(imp, l)
	    chi3 = c3(imp, l)
	    do m = 1, iymx
	      jm = j - m
	      jm2 = j + m - 1
	      do n = 1, k-1
	        kn = k - n
	        if( vel(l,m,n).ne.0. ) then
		  if( chi3.eq.0. ) then
		    piece1 = pon(il, jm, kn)
                    piece2 = pon(il, jm2, kn)
		    piece3 = 0.d0
	     	    piece4 = 0.d0
		  else
                    t  = dble(kn)
                    c2 = dble(jm)
                    tn33 = ker33s(chi1, c2, chi3, t, factor)
                    tn11 = ker11s(chi1, c2, chi3, t, factor)
                    sh31 = ker31s(chi1, c2, chi3, t, factor)
                    piece1 = csh(imp, l)*sh31 - ctn(imp, l)*(tn33-tn11)
                    piece3 = tn33*dnu3(imp, l)**2 + &
                                tn11*dnu1(imp, l)**2 + 2.0d0*sh31*ctn(imp, l)

                    c2 = dble(jm2)
                    tn33 = ker33s(chi1, c2, chi3, t, factor)
                    tn11 = ker11s(chi1, c2, chi3, t, factor)
                    sh31 = ker31s(chi1, c2, chi3, t, factor)
                    piece2 = csh(imp, l)*sh31 - ctn(imp, l)*(tn33-tn11)
                    piece4 = tn33*dnu3(imp, l)**2 + &
                                tn11*dnu1(imp, l)**2 + 2.0d0*sh31*ctn(imp, l)
	  	  endif
                  if(piece1.eq.0.) exit
	          dtau = dtau + vel(l,m,n)*(piece1 + piece2)
              dnor = dnor + vel(l,m,n)*(piece3 + piece4)
	        endif
	      enddo
	    enddo
	  enddo
        endif
        dtau = tau0(i,j) + dtau0(i,j) + const*dtau
        dnor = nor0(i,j) + const*dnor

	wtry = w(i,j) + dw0(i,j)
!        dtp = cohesive + fs * dnor 
!        dtr = fd * dnor 
	dtp = tp(i,j)
	dtr = tr(i,j)
        da  = (dtp - dtr)/dc(i,j)
        if( wtry.ge.dc(i,j) ) then
          dsigma = dtr
        else
          dsigma = dtp - da * w(i,j)
        endif

        ivmpi(imp, j) = iv(i, j)
        if( ivmpi(imp, j).eq.0 ) then
          vmpi( imp, j ) = 0.0d0
          if( dtau.gt.dsigma ) ivmpi(imp, j) = 2
        else
          ivmpi(imp, j) = 1
          vmpi(imp, j) = (dsigma - dtau)/(const*p000 + da*dt)
          if( (wtry + vmpi(imp,j)*dt ).gt.dc(i,j) ) then
            vmpi(imp, j) = ( dtr - dtau)/(const*p000)
          endif
          if(vmpi(imp, j).lt.0.) then
            vmpi(imp, j)  = 0.0d0
            ivmpi(imp, j) = 0
            if( dtau.gt.dsigma ) ivmpi(imp, j) = 2
          endif
        endif

        str_mpi(imp, j) = dtau + const*p000*vmpi(imp, j)
        nor_mpi(imp, j) = dnor
        tau_mpi(imp, j) = dtau
        w_mpi(imp, j)   = w(i, j) + vmpi(imp, j)*dt

        if( (w_mpi(imp, j) + dw0(i,j)).ge.dc(i, j) ) then
          sig_mpi(imp, j) = dtr
        else
          sig_mpi(imp, j) = dtp - da*w_mpi(imp, j)
        endif
  enddo
!  write(6,*) "end iteration", nrank + 1
!
! ASSIMILATION BETWEEN MASTER CPU AND THE OTHERS
!
  IF( nrank.eq.master ) THEN
    do imp = 1, mpmx
      ncell = imp2i(nsize, 0, mpmx, imp)
      if( ncell.gt.ndata2 ) exit
      i = 1 + (ncell-1)/mpmy
      j = mod(ncell-1, mpmy) + 1

!      i = (imp - 1)*nsize + nrank + 1
!      if( i.gt.ixmx ) exit
      vel(i, j, k) = vmpi(imp, j)
      iv(i, j)     = ivmpi(imp, j)
      stress(i, j) = str_mpi(imp, j)
      normal(i, j) = nor_mpi(imp, j)
      tau(i, j)    = tau_mpi(imp, j)
      w(i, j)      = w_mpi(imp, j)
      sigma(i, j)  = sig_mpi(imp, j)
    enddo
    do iproc = 1, nsize - 1
      call MPI_RECV( vmpi,     mpmx*mpmy, MPI_REAL8, iproc, 1, &
                MPI_COMM_WORLD, mpist, ierr )
      call MPI_RECV( ivmpi,    mpmx*mpmy, MPI_INTEGER, iproc, 2, &
                MPI_COMM_WORLD, mpist, ierr )
      call MPI_RECV( str_mpi,  mpmx*mpmy, MPI_REAL8, iproc, 3, &
                MPI_COMM_WORLD, mpist, ierr )
      call MPI_RECV( nor_mpi,  mpmx*mpmy, MPI_REAL8, iproc, 4, &
                MPI_COMM_WORLD, mpist, ierr )
      call MPI_RECV( tau_mpi,  mpmx*mpmy, MPI_REAL8, iproc, 5, &
                MPI_COMM_WORLD, mpist, ierr )
      call MPI_RECV( w_mpi,    mpmx*mpmy, MPI_REAL8, iproc, 6, &
                MPI_COMM_WORLD, mpist, ierr )
      call MPI_RECV( sig_mpi,  mpmx*mpmy, MPI_REAL8, iproc, 7, &
                MPI_COMM_WORLD, mpist, ierr )

      do imp = 1, mpmx
!        i = (imp - 1)*nsize + iproc + 1
!        if(i.gt.ixmx) exit
	ncell = imp2i(nsize, iproc, mpmx, imp)
	if( ncell.gt.ndata2 ) exit
	i = 1 + (ncell-1)/mpmy
	j = mod(ncell-1, mpmy) + 1

        vel(i, j, k)  = vmpi(imp, j)
        iv(i, j)      = ivmpi(imp, j)
        stress(i, j)  = str_mpi(imp, j)
        normal(i, j)  = nor_mpi(imp, j)
        tau(i, j)     = tau_mpi(imp, j)
        w(i, j)       = w_mpi(imp, j)
        sigma(i, j)   = sig_mpi(imp, j)
      enddo
    enddo
!
! INITIAL RUPTURE PROPAGATION (ref DAY)
!
    if(istatic.eq.0) then
      travel = terminal*dble(k)/factor
      do i=1, ixmx
        rad = sqrt((x(i)-x(ihypo))**2 + (z(i)-z(ihypo))**2)
        if(rad.gt.travel.and.iv(i,jhypo).ne.0.and.travel.gt.1) then
          istatic = istatic + 1
        endif
      enddo
      if( k.lt.20.) then
	istatic = 0
      endif

      if(istatic.eq.0) then
        do i=1, ixmx
          do j=1, iymx
            rad = sqrt((x(i)-x(ihypo))**2  &
                   + (y(j)-y(jhypo))**2 + (z(i)-z(ihypo))**2)
            if(rad.le.travel.and.iv(i,j).eq.0)  then
              if(tp(i,j).ge.tr(i,j)) then
		dw0(i,j) = dc(i,j) 
		if( tau(i,j).le.tr(i,j).and.dtau0(i,j).eq.0. ) then
  		  dtau0(i,j) = tr(i,j) - tau(i,j) + 3.0d0
		endif
              else
                dtau0(i,j) = sigma(i,j) - tau(i,j)
              endif
              iv(i,j) = 2
            endif
          enddo
        enddo

      endif
    endif
  ELSE
    call MPI_SEND( vmpi,     mpmx*mpmy, MPI_REAL8, master, 1, &
                MPI_COMM_WORLD, ierr )
    call MPI_SEND( ivmpi,  mpmx*mpmy, MPI_INTEGER, master, 2, &
                MPI_COMM_WORLD, ierr )
    call MPI_SEND( str_mpi,  mpmx*mpmy, MPI_REAL8, master, 3, &
                MPI_COMM_WORLD, ierr )
    call MPI_SEND( nor_mpi,  mpmx*mpmy, MPI_REAL8, master, 4, &
                MPI_COMM_WORLD, ierr )
    call MPI_SEND( tau_mpi,  mpmx*mpmy, MPI_REAL8, master, 5, &
                MPI_COMM_WORLD, ierr )
    call MPI_SEND( w_mpi,    mpmx*mpmy, MPI_REAL8, master, 6, &
                MPI_COMM_WORLD, ierr )
    call MPI_SEND( sig_mpi,  mpmx*mpmy, MPI_REAL8, master, 7, &
                MPI_COMM_WORLD, ierr )
  ENDIF
  do i = -ixmx, ixmx
    t = dble(k)
    chi1 = dble(i)
    do j = -2*iymx, 2*iymx
      c2 = dble(j)
      pon(i,j,k) =  ker31s(chi1, c2, 0.0d0, t, factor)
    enddo
  enddo

!
! OUTPUT
!
  IF( nrank.eq.master ) THEN
    if (k.lt.1000) then
      write(number,'(i3.3)') k
      name3 = dir(1:ndir)//'/snap'//number(1:3)//'.dat'
    else
      write(number,'(i4.4)') k
      name3 = dir(1:ndir)//'/snap'//number(1:4)//'.dat'
    endif

    open(13, file=name3)
    do i=1, ixmx
      do j=1, iymx
        write(13, 105) real(i), real(j), vp*vel(i,j,k), w(i,j), 
     &          stress(i,j), sigma(i,j), normal(i,j)
 105    format(f8.1, 2x, f8.1, 2x, 5f11.6)
      enddo
    enddo
    close(13)
  ENDIF
enddo

23 continue


call MPI_FINALIZE( ierr )
END
!
!

INTEGER FUNCTION imp2i(nsize, nrank, mpmx, imp)
IMPLICIT NONE
INTEGER nsize, nrank, mpmx, imp

! order of 111222333444
! mpmx: number of data for each CPU
!imp2i = nrank*mpmx + imp

! order of 123412341234
! nsize: number of CPU
imp2i = (imp-1)*nsize + nrank + 1

return
END

