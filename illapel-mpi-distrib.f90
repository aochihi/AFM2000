! this is a sample program to understand the 3D boundary integral 
! equation method (BIEM) by Fukuyama and Madariaga (BSSA, 1995, 1998) and 
! Aochi, Fukuyama and Matsu'ura (Pageoph, 2000).
! This is used for preliminary simulations of 2011 Tohoku-oki  and 
! 2015 Illapel eathquakes. 
!
! For simplicity, a planar fault is supposed.
! Output is in physical unit.
!
! The code execute the spatio-temporal convolution directly. 
! Code is adjusted for distribution on the 5th December 2023.
!
! Please cite: 
! Aochi, H., E. Fukuyama and M. Matsu'ura (2000)
! Spontaneous rupture propagation on a non-planar fault in 3D 
! Elastic Medium, Pure appl. Geophy. 157, 2003-2027.
! https://doi.org/10.1007/PL00001072
!
IMPLICIT NONE
INCLUDE 'mpif.h'
INTEGER ixmn, ixmx, itmx, iymn, iymx
INTEGER master, nrank, nsize
REAL(8) :: xmax, ymax, ds, xmin, ymin
!!!
!!! MODIFIABLE PARAMETERS (*1)
!!!
! MODEL FAULT AREA (X,Y) = (xmin:xmax, ymin:ymax) in km
PARAMETER (xmin = -20.d0, xmax = 70.0d0, ymax = 20.0d0, ymin = -120.d0)
! BIEM element size in km
PARAMETER (ds = 2.0d0)
! BIEM time steps / master node = 0 for MPI
PARAMETER (itmx = 400, master = 0)
!!! END (*1)

REAL(8), DIMENSION(:, :, :), ALLOCATABLE :: vel, pon
REAL(8), DIMENSION(:, :), ALLOCATABLE :: tau0, tp, tr, a, dc, &
        stress, sigma, w, tau, &
        vmpi, w_mpi, str_mpi, sig_mpi, tau_mpi
REAL(8) :: pi, mu, const, factor, fac, &
        x, y, tp0, tr0, dc0, t0, prm1, &
        p000, ker31s, dtau, dsigma, alpha
INTEGER, DIMENSION(:, :), ALLOCATABLE :: iv, ivmpi
REAL(8) :: dt, rad, r0, x0, y0
REAL(8) :: c1, c2, t, sh31, piece1
REAL(8) :: a1, b1, xhypo, yhypo, r2, x2, y2, r3, x3, y3, r4, x4, y4, dist, &
        r5, x5, y5
INTEGER :: i, j, k, l, m, n, ndata, ndir, istatic, izone, nvel
INTEGER :: mpist(MPI_STATUS_SIZE), ndata3, ndata2, isend, mpmx, mpmy, &
        imp, jmp, iproc, ierr, ndatax, ndatay, idata, ldata
EXTERNAL ker31s
CHARACTER*40 name0, name2, name3, name4, name6, dir
CHARACTER*5  number

call MPI_Init( ierr )                          
call MPI_Comm_rank( MPI_COMM_WORLD, nrank, ierr )
call MPI_Comm_size( MPI_COMM_WORLD, nsize, ierr ) 
write(*, *) "I am alive! I'm ", nrank+1, " of ", nsize

ixmx = int(xmax/ds)
iymx = int(ymax/ds)
ixmn = int(xmin/ds)
iymn = int(ymin/ds)

!
! FIXED PARAMETER
!
        pi = acos(-1.0d0)
        factor = 2.0d0
!!!
!!! MODIFIABLE PARAMETERS (*2)
!!!
        dir = 'test1'
! ridigity in GPa
        mu = 32.4d0
! Vp in km/s
        alpha = 6.2d0
! peak strength of slip-weakening friction in MPa
        tp0 = 5.00d0
! residual stress of slip-weakening in MPa
        tr0 = 0.0d0
! initial shear stress in MPa
        t0 = 3.0d0
! Dc value in meter
        dc0 = 1.60d0
! initial crack size radius in km
        r0 = 8.0d0
!!! END (*2)

! unit (mu) = mu [GPa]/tb[MPa] dc0[m]/ds [km] = 1
! dc0 should be normalized by 0.001 ds.
!
!
! SETTING
!
! fix for calculation
dt = ds/factor
ndir = index(dir,' ')-1
p000 = ker31s(0.0d0, 0.0d0, 0.0d0, 0.0d0, factor)
const = sqrt(3.0d0)/(4.0d0*pi)*mu

if (nrank.eq.master) then
  write(*,*) dir
  write(*,*) ixmn, ixmx, iymn, iymx, itmx
endif

ndatax = ixmx - ixmn + 1
ndatay = iymx - iymn + 1
ndata2 = ndatax * ndatay
ndata3 = ndata2*(itmx + 1)

ALLOCATE( vel(ixmn:ixmx, iymn:iymx, 0:itmx) )
ALLOCATE( pon((ixmn-ixmx):(ixmx-ixmn), (iymn-iymx):(iymx-iymn), 0:itmx) )
ALLOCATE(tau0(ixmn:ixmx, iymn:iymx),     tp(ixmn:ixmx, iymn:iymx), &
       stress(ixmn:ixmx, iymn:iymx), &
           tr(ixmn:ixmx, iymn:iymx), &
            a(ixmn:ixmx, iymn:iymx),     dc(ixmn:ixmx, iymn:iymx), &
        sigma(ixmn:ixmx, iymn:iymx),      w(ixmn:ixmx, iymn:iymx), &
           iv(ixmn:ixmx, iymn:iymx),    tau(ixmn:ixmx, iymn:iymx) )
!!!
!!! MODIFIABLE PARAMETERS (*3)
!!!
! in the following, three paches of different size are superposed. 
! the rupture is expected to grow up from the hypocenter position, 
! according to the Dc-scaling with patch size (Ide & Aochi, JGR, 2005).
!
! hypocenter position (X, Y) on fault in km 
xhypo = 0.0d0
yhypo = 0.0d0

! To give a heterogeneity. Center of heterogeneity geometry (X, Y) in km
x0 = 30.0d0
y0 = -60.0d0

do idata = 1, ndata2
! mapping "idata => (i, j)"
  i = (idata-1)/ndatay + ixmn
  j = mod(idata-1, ndatay) + iymn

    x = ds*(i-0.5)
    y = ds*(j-0.5)
    w(i,j) = 0.0d0
! background condition
    tau0(i,j) = t0
    tp(i,j) = tp0
    tr(i,j) = tr0
    izone = -1
    dc(i,j) = 2.*dc0
! patch Rank 1 on (x0, y0)
    a1 = 50.0d0
    b1 = 50.0d0
    izone = 0
    dist = (x-x0)**2/a1**2 + (y-y0)**2/b1**2
    if (dist.le.1.) then
      izone = 1
      dc(i,j) = dc0
    endif

! patch rank 2 on (xhypo+10, yhypo+10) 
    r3 =  25.0d0
    x3 =  xhypo + 10.0d0
    y3 =  yhypo - 10.d0
    dist = sqrt((x-x3)**2 + (y-y3)**2)
      if (dist.le.r3) then
        izone = 2
        dc(i,j) = dc0/2.d0
        tp(i,j) = tp(i,j)*1.0
      endif

! patch rank 3 (new) added for rupture growth
    r3 =  12.5d0
    x3 =  xhypo 
    y3 =  yhypo 
      dist = sqrt((x-x3)**2 + (y-y3)**2)
      if (dist.le.r3) then
        izone = 3
        dc(i,j) = dc0/2.d0**2
        tp(i,j) = tp(i,j)*1.0
      endif

! background out of patches is barrior
    if( izone.lt.0 ) then
      tp(i,j) = 100.0d0
    endif

! for calculation
    a(i,j)  = (tp(i,j) - tr(i,j))/dc(i,j)
    sigma(i,j) = tp(i,j)
enddo
!!! END (*3)

do idata=1, ndata2
    i = (idata-1)/ndatay + ixmn
    j = mod(idata-1, ndatay) + iymn

    x = ds*(i-0.5)
    y = ds*(j-0.5)
    rad = sqrt((x - xhypo)**2 + (y - yhypo)**2)

    if ( rad.le.r0 ) then
      tp(i,j) = tr(i,j)
      dc(i,j) = 0.0d0
      a(i,j) = 0.0d0
      sigma(i,j) = tp(i,j)
      iv(i,j) = 2
    endif

enddo

IF( nrank.eq.master ) THEN
        name0 = dir(1:ndir)//'/struct.dat'
        name2 = dir(1:ndir)//'/hoge2.dat'

        open(12, file=name2)
        write(12, '(5i10)') ixmn, ixmx, iymn, iymx, itmx
        write(12, '(3f9.4)') tp0, tr0, t0
        write(12, '(4f12.5)') ds, dt, dc0, r0
        write(12, '(2f12.5)') mu, alpha
        write(12, '(2f12.2)') xhypo, yhypo
        write(12, '(2f12.2)') x0, y0
        close(12)

        open(10, file=name0)
        ndata = 1
        do i=ixmn, ixmx
          x = ds*(i-0.5)
          do j=iymn, iymx
            y = ds*(j-0.5)
            write(10,'(i5, 6f10.3)') ndata, x, y, &
                tau0(i,j), tp(i,j), tr(i,j), dc(i,j)
            ndata = ndata + 1
          enddo
        enddo
close(10)
ENDIF

mpmx = ndata2/nsize + 1
mpmy = 1

ALLOCATE (   vmpi(mpmx, mpmy),    ivmpi(mpmx, mpmy), &
          str_mpi(mpmx, mpmy), &
          tau_mpi(mpmx, mpmy),    w_mpi(mpmx, mpmy), &
          sig_mpi(mpmx, mpmy)  )

!	
! ITERATION
!
do k = 1, itmx

  if(nrank.eq.master) write(6,*) k
  do imp = 1, mpmx
    jmp = 1
    idata = nrank*mpmx + imp
    if ( idata.gt.ndata2) exit
    i = (idata-1)/ndatay + ixmn
    j = mod(idata-1, ndatay) + iymn
    fac = factor

      dtau = 0.0d0
      dsigma = 0.0d0
      if(k.ne.1) then

        do ldata = 1, ndata2
            l = (ldata-1)/ndatay + ixmn
            m = mod(ldata-1, ndatay) + iymn
            do n = 1, k-1
              if( vel(l,m,n).ne.0. ) then
                piece1 = pon(i-l, j-m, k-n)
                if(piece1.eq.0.) exit
                dtau = dtau +  vel(l,m,n)*piece1
              endif
            enddo
        enddo
      endif
      dtau   = tau0(i,j) + const*dtau
      dsigma = sigma(i,j)

      ivmpi(imp, jmp) = iv(i, j)
      if(ivmpi(imp,jmp).eq.0) then
        vmpi(imp,jmp) = 0.0d0
        if( dtau.gt.dsigma ) ivmpi(imp, jmp) = 2
      else
        ivmpi(imp, jmp) = 1
        vmpi(imp,jmp)  = (dsigma - dtau)/(const*p000 + a(i,j)*dt)
        if((w(i,j)+vmpi(imp,jmp)*dt).gt.dc(i,j)) then
          vmpi(imp,jmp) = ( tr(i,j) - dtau)/(const*p000)
        endif

!        write(*,*) nrank, imp, idata, i, j, vmpi(imp,jmp)
        if(vmpi(imp,jmp).lt.0.) then
          vmpi(imp, jmp)  = 0.0d0
          ivmpi(imp, jmp) = 0
          if( dtau.gt.dsigma ) ivmpi(imp, jmp) = 2
        endif
      endif

      str_mpi(imp, jmp)  = dtau + const*p000*vmpi(imp, jmp)
      tau_mpi(imp, jmp)  = dtau
      w_mpi(imp, jmp)    = w(i, j) + vmpi(imp, jmp)*dt
      if( w_mpi(imp, jmp).gt.dc(i, j) ) then
        sig_mpi(imp, jmp) = tr(i,j)
      else
        sig_mpi(imp, jmp) = tp(i,j) - a(i,j)*w_mpi(imp,jmp)
      endif
  enddo
  if(nrank.eq.master) write(6,*) "end convolution"
!
! SYNCHRONIZING
!
  IF( nrank.eq.master ) THEN
    do imp=1, mpmx
      jmp = 1
      idata = nrank*mpmx + imp
      if(idata.gt.ndata2) exit
      i = (idata-1)/ndatay + ixmn
      j = mod(idata-1, ndatay) + iymn
        vel(i, j, k)  = vmpi(imp, jmp)
        iv(i, j)      = ivmpi(imp, jmp)
        stress(i, j)  = str_mpi(imp, jmp)
        tau(i, j)     = tau_mpi(imp, jmp)
        w(i, j)       = w_mpi(imp, jmp)
        sigma(i, j)   = sig_mpi(imp, jmp)
    enddo
    do iproc = 1, nsize-1
      call MPI_RECV( vmpi,     mpmx*mpmy, MPI_REAL8, iproc, 1, &
                MPI_COMM_WORLD, mpist, ierr )
      call MPI_RECV( ivmpi,    mpmx*mpmy, MPI_INTEGER, iproc, 2, &
                MPI_COMM_WORLD, mpist, ierr )
      call MPI_RECV( str_mpi,  mpmx*mpmy, MPI_REAL8, iproc, 3, &
                MPI_COMM_WORLD, mpist, ierr )
      call MPI_RECV( tau_mpi,  mpmx*mpmy, MPI_REAL8, iproc, 5, &
                MPI_COMM_WORLD, mpist, ierr )
      call MPI_RECV( w_mpi,    mpmx*mpmy, MPI_REAL8, iproc, 6, &
                MPI_COMM_WORLD, mpist, ierr )
      call MPI_RECV( sig_mpi,  mpmx*mpmy, MPI_REAL8, iproc, 7, &
                MPI_COMM_WORLD, mpist, ierr )

      isend = iproc
      do imp = 1, mpmx
        jmp = 1
        idata = isend*mpmx + imp
        if(idata.gt.ndata2) exit
        i = (idata-1)/ndatay + ixmn 
        j = mod(idata-1, ndatay) + iymn
        vel(i, j, k)  = vmpi(imp, jmp)
        iv(i, j)      = ivmpi(imp, jmp)
        stress(i, j)  = str_mpi(imp, jmp)
        tau(i, j)     = tau_mpi(imp, jmp)
        w(i, j)       = w_mpi(imp, jmp)
        sigma(i, j)   = sig_mpi(imp, jmp)
      enddo
    enddo

  ELSE
    call MPI_SEND( vmpi,     mpmx*mpmy, MPI_REAL8, master, 1, &
                MPI_COMM_WORLD, ierr )
    call MPI_SEND( ivmpi,    mpmx*mpmy, MPI_INTEGER, master, 2, &
                MPI_COMM_WORLD, ierr )
    call MPI_SEND( str_mpi,  mpmx*mpmy, MPI_REAL8, master, 3, &
                MPI_COMM_WORLD, ierr )
    call MPI_SEND( tau_mpi,  mpmx*mpmy, MPI_REAL8, master, 5, &
                MPI_COMM_WORLD, ierr )
    call MPI_SEND( w_mpi,    mpmx*mpmy, MPI_REAL8, master, 6, &
                MPI_COMM_WORLD, ierr )
    call MPI_SEND( sig_mpi,  mpmx*mpmy, MPI_REAL8, master, 7, &
                MPI_COMM_WORLD, ierr )
  ENDIF

  t = dble(k)
  do i = ixmn-ixmx, ixmx-ixmn
    fac = factor
    c1 = dble(i)
    do j = iymn-iymx, iymx-iymn
      c2 = dble(j)
      pon(i,j,k) =  ker31s(c1, c2, 0.0d0, t, fac)
    enddo
  enddo
!
  call MPI_BCAST( vel,     ndata3, MPI_REAL8, master, MPI_COMM_WORLD, ierr ) 
  call MPI_BCAST( iv,     ndata2, MPI_INTEGER, master, MPI_COMM_WORLD, ierr ) 
  call MPI_BCAST( stress,  ndata2, MPI_REAL8, master, MPI_COMM_WORLD, ierr ) 
  call MPI_BCAST( tau,     ndata2, MPI_REAL8, master, MPI_COMM_WORLD, ierr ) 
  call MPI_BCAST( w,       ndata2, MPI_REAL8, master, MPI_COMM_WORLD, ierr ) 
  call MPI_BCAST( sigma,   ndata2, MPI_REAL8, master, MPI_COMM_WORLD, ierr )

  IF( nrank.eq.master ) THEN
          write(number,'(i3.3)') k
          name3 = dir(1:ndir)//'/output'//number(1:3)//'.dat'
          open(13, file=name3)
  ENDIF

  nvel = 0   
  do i=ixmn, ixmx
    do j=iymn, iymx 
 
      if (vel(i,j,k).ne.0.) nvel = nvel + 1 
      IF( nrank.eq.master ) THEN
        write(13, 105) ds*(i-0.5), ds*(j-0.5), vel(i,j,k)*alpha, w(i,j), 
                stress(i,j), sigma(i,j)
 105          format(f7.1, 1x, f7.1, 1x, 4f13.8)
      ENDIF
    enddo
  enddo

  IF( nrank.eq.master ) THEN
    close(13)
    close(14)
  ENDIF

  if (nvel.eq.0) exit
enddo

call MPI_FINALIZE( ierr )
END
