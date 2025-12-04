        DOUBLE PRECISION FUNCTION ker32s(c1, c2, c3, t, w)
        DOUBLE PRECISION c1, c2, c3, t, w, t1
        DOUBLE PRECISION kernel1, kernel2
        DOUBLE PRECISION comp1, comp2
        DOUBLE PRECISION p321, q321
        EXTERNAL p321, q321

        t1 = t + 1.0d0
              kernel1 = p321(c1+0.50d0, c2, c3, t1, w)
     &                - p321(c1-0.50d0, c2, c3, t1, w)
              kernel2 = p321(c1+0.50d0, c2, c3, t, w)
     &                - p321(c1-0.50d0, c2, c3, t, w)
              comp1 = kernel1 - kernel2

              kernel1 = q321(c1+0.50d0, c2, c3, t1, w)
     &                - q321(c1-0.50d0, c2, c3, t1, w)
              kernel2 = q321(c1+0.50d0, c2, c3, t, w)
     &                - q321(c1-0.50d0, c2, c3, t, w)
              comp2 = kernel1 - kernel2

              ker32s = comp1 + comp2 

        return
        END
c
        DOUBLE PRECISION FUNCTION p321(chi1, chi2, chi3, dt, w)
        DOUBLE PRECISION chi1, chi2, chi3, dt, w
	DOUBLE PRECISION c2a, c2b, r321b
        DOUBLE PRECISION c, rc2, rc, x21, x1, r2, r
	DOUBLE PRECISION chi1d2, m, m2

        c = 1.0d0/(w*sqrt(3.0d0))
        rc  = c*dt
        rc2 = rc**2.
	chi1d2 = chi1**2. + chi3**2.
        x21 = rc2 - chi1d2
	if(x21.le.0.) then
	  p321 = 0.0d0
	else
	  c2a = chi2 + 0.50d0
	  c2b = chi2 - 0.50d0
	  x1 = sqrt(x21)
	  if(abs(c2a).ge.x1.and.abs(c2b).ge.x1) then
	    p321 = 0.0d0
	    r321b = 0.0d0
	  else if(abs(c2a).le.x1.and.abs(c2b).le.x1) then
            r2 = chi1d2 + c2a**2.
	    m2 = rc2/r2
            r  = sqrt(r2)
	    m = rc/r
	    p321 = m*(1.0d0 - 2.0d0*m2/3.0d0)
            r321b = -(chi3/r)**2.*m*(1.0d0 - m2)

            r2 = chi1d2 + c2b**2.
	    m2 = rc2/r2
            r  = sqrt(r2)
            m = rc/r
            p321 = p321 - m*(1.0d0 - 2.0d0*m2/3.0d0)
            r321b = r321b + (chi3/r)**2.*m*(1.0d0 - m2)
            r321b = 2.0d0*r321b
	  else if(c2a.gt.x1.and.c2b.lt.x1) then
            r2 = chi1d2 + c2b**2.
            m2 = rc2/r2
            r  = sqrt(r2)
            m = rc/r
            p321 = 1.0d0/3.0d0 - m*(1.0d0 - 2.0d0*m2/3.0d0)
            r321b = 2.0d0*(chi3/r)**2.*m*(1.0d0 - m2)
	  else if(c2a.gt.(-x1).and.c2b.lt.(-x1)) then
            r2 = chi1d2 + c2a**2.
            m2 = rc2/r2
            r  = sqrt(r2)
            m = rc/r
            p321 = -1.0d0/3.0d0 + m*(1.0d0 - 2.0d0*m2/3.0d0)
            r321b = -2.0d0*(chi3/r)**2.*m*(1.0d0 - m2)
	  else 
            pause 'bad if'
          endif
	  if(chi3.ne.0.) p321 = p321 + r321b
        endif

	return
	END
c
        DOUBLE PRECISION FUNCTION q321(chi1, chi2, chi3, dt, w)
        DOUBLE PRECISION chi1, chi2, chi3, dt, w
	DOUBLE PRECISION c2a, c2b, r321a, const
        DOUBLE PRECISION c, p, p3, rc2, rc, x21, x1, r2, r
	DOUBLE PRECISION chi1d2, m, m2

        c = 1.0d0/w
        p = 1.0d0/sqrt(3.0d0)
	p3 = p**3.

        rc  = c*dt
        rc2 = rc**2.
	chi1d2 = chi1**2. + chi3**2.
        x21 = rc**2. - chi1d2
	if(x21.le.0.) then
	  q321 = 0.0d0
	else
	  c2a = chi2 + 0.50d0
	  c2b = chi2 - 0.50d0
	  x1 = sqrt(x21)
	  if(abs(c2a).ge.x1.and.abs(c2b).ge.x1) then
	    q321 = 0.0d0
	    r321a = 0.0d0
	  else if(abs(c2a).le.x1.and.abs(c2b).le.x1) then
            r2 = chi1d2 + c2a**2.
	    m2 = rc2/r2
            r  = sqrt(r2)
	    m = rc/r
	    q321 = -m*(1.0d0 - m2/3.0d0)
            r321a = -(chi3/r)**2.*m*(1.0d0 - m2)

            r2 = chi1d2 + c2b**2.
            m2 = rc2/r2
            r  = sqrt(r2)
            m = rc/r
	    const = 2.0d0*p3
            q321 = q321 + m*(1.0d0 - m2/3.0d0)
	    q321 = const*q321
            r321a = r321a + (chi3/r)**2.*m*(1.0d0 - m2)
            r321a = const*r321a
	  else if(c2a.gt.x1.and.c2b.lt.x1) then
	    r2 = chi1d2 + c2b**2.
            m2 = rc2/r2
            r  = sqrt(r2)
            m = rc/r
	    const = 2.0d0*p3
            q321 = -const*(2.0d0/3.0d0 - m*(1.0d0 - m2/3.0d0))
            r321a = const*(chi3/r)**2.*m*(1.0d0 - m2)
	  else if(c2a.gt.(-x1).and.c2b.lt.(-x1)) then
	    r2 = chi1d2 + c2a**2.
            m2 = rc2/r2
            r  = sqrt(r2)
            m = rc/r
	    const = 2.0d0*p3
            q321 = const*(2.0d0/3.0d0 - m*(1.0d0 - m2/3.0d0))
            r321a = -const*(chi3/r)**2.*m*(1.0d0 - m2)
	  else
	    pause 'bad if'
          endif
	  if(chi3.ne.0.) q321 = q321 - r321a
        endif

        return
        END
c
        DOUBLE PRECISION FUNCTION r321b(chi1, chi2, chi3, dt, w)
        DOUBLE PRECISION chi1, chi2, chi3, dt, w
	DOUBLE PRECISION c2a, c2b
        DOUBLE PRECISION c, rc2, rc, x21, x1, r2, r
        DOUBLE PRECISION chi1d2, m, m2

	if(chi3.eq.0.) then
	  r321b = 0.0d0
	  return
	endif
	 
        c = 1.0d0/(w*sqrt(3.0d0))
        rc  = c*dt
        rc2 = rc**2.
        chi1d2 = chi1**2. + chi3**2.
        x21 = rc2 - chi1d2
	if(x21.le.0.) then
	  r321b = 0.0d0
	else
	  c2a = chi2 + 0.50d0
	  c2b = chi2 - 0.50d0
	  x1 = sqrt(x21)
	  if(abs(c2a).ge.x1.and.abs(c2b).ge.x1) then
	    r321b = 0.0d0
	  else if(abs(c2a).le.x1.and.abs(c2b).le.x1) then
            r2 = chi1d2 + c2a**2.
	    m2 = rc2/r2
            r  = sqrt(r2)
	    m  = rc/r
	    r321b = -(chi3/r)**2.*m*(1.0d0 - m2)

            r2 = chi1d2 + c2b**2.
            m2 = rc2/r2
            r  = sqrt(r2)
            m  = rc/r
            r321b = r321b + (chi3/r)**2.*m*(1.0d0 - m2)
	    r321b = 2.0d0*r321b
	  else if(c2a.gt.x1.and.c2b.lt.x1) then
            r2 = chi1d2 + c2b**2.
            m2 = rc2/r2
            r  = sqrt(r2)
            m  = rc/r
            r321b = 2.0d0*(chi3/r)**2.*m*(1.0d0 - m2)
	  else if(c2a.gt.(-x1).and.c2b.lt.(-x1)) then
            r2 = chi1d2 + c2a**2.
            m2 = rc2/r2
            r  = sqrt(r2)
            m  = rc/r
            r321b = -2.0d0*(chi3/r)**2.*m*(1.0d0 - m2)
	  else
	    pause 'bad if'
          endif
        endif

        return
        END
c
        DOUBLE PRECISION FUNCTION r321a(chi1, chi2, chi3, dt, w)
        DOUBLE PRECISION chi1, chi2, chi3, dt, w
	DOUBLE PRECISION c2a, c2b
        DOUBLE PRECISION c, p, p3, rc2, rc, x21, x1, r2, r
        DOUBLE PRECISION chi1d2, m, m2

        if(chi3.eq.0.) then
          r321a = 0.0d0
          return
        endif

        c = 1.0d0/w
	p = 1.0d0/sqrt(3.0d0)
	p3 = p**3.

        rc  = c*dt
        rc2 = rc**2.
        chi1d2 = chi1**2. + chi3**2.
        x21 = rc2 - chi1d2
	if(x21.le.0.) then
	  r321a = 0.0d0
	else
	  c2a = chi2 + 0.50d0
	  c2b = chi2 - 0.50d0
	  x1 = sqrt(x21)
	  if(abs(c2a).ge.x1.and.abs(c2b).ge.x1) then
	    r321a = 0.0d0
	  else if(abs(c2a).le.x1.and.abs(c2b).le.x1) then
            r2 = chi1d2 + c2a**2.
	    m2 = rc2/r2
            r  = sqrt(r2)
	    m  = rc/r
	    r321a = -(chi3/r)**2.*m*(1.0d0 - m2)

            r2 = chi1d2 + c2b**2.
            m2 = rc2/r2
            r  = sqrt(r2)
            m  = rc/r
            r321a = r321a + (chi3/r)**2.*m*(1.0d0 - m2)
	    r321a = 2.0d0*p3*r321a
	  else if(c2a.gt.x1.and.c2b.lt.x1) then
            r2 = chi1d2 + c2b**2.
            m2 = rc2/r2
            r  = sqrt(r2)
            m  = rc/r
            r321a = 2.0d0*p3*(chi3/r)**2.*m*(1.0d0 - m2)
	  else if(c2a.gt.(-x1).and.c2b.lt.(-x1)) then
            r2 = chi1d2 + c2a**2.
            m2 = rc2/r2
            r  = sqrt(r2)
            m  = rc/r
            r321a = -2.0d0*p3*(chi3/r)**2.*m*(1.0d0 - m2)
	  else
	    pause 'bad if'
          endif
        endif

        return
        END
c
