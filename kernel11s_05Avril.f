        REAL(8) FUNCTION ker11s(c1, c2, c3, t, w)
	IMPLICIT NONE
        REAL(8) c1, c2, c3, t, w, t1
        REAL(8) kernel1, kernel2
        REAL(8) comp1, comp2
        REAL(8) p111, q111
        EXTERNAL p111, q111

        t1 = t + 1.0d0
              kernel1 = p111(c1+0.50d0, c2, c3, t1, w)
     &                - p111(c1-0.50d0, c2, c3, t1, w)
              kernel2 = p111(c1+0.50d0, c2, c3, t, w)
     &                - p111(c1-0.50d0, c2, c3, t, w)
              comp1 = kernel1 - kernel2

              kernel1 = q111(c1+0.50d0, c2, c3, t1, w)
     &                - q111(c1-0.50d0, c2, c3, t1, w)
              kernel2 = q111(c1+0.50d0, c2, c3, t, w)
     &                - q111(c1-0.50d0, c2, c3, t, w)
              comp2 = kernel1 - kernel2

              ker11s = comp1 + comp2 

        return
        END
c
        REAL(8) FUNCTION p111(chi1, chi2, chi3, dt, w)
	IMPLICIT NONE
        REAL(8) chi1, chi2, chi3, dt, w
        REAL(8) c2a, c2b, e111, r111b
        REAL(8) c, rc2, rc, x21, x1, r2, r
        REAL(8) chi1d2, fb1, fb2, fb3, fb4, fb5, fb6
        REAL(8) i31, j1, k1, l2, m, n11, const

        if(chi3.eq.0.) then
          p111  = 0.0d0
	  e111 =  0.0d0
	  r111b = 0.0d0
          return
        endif

        c = 1.0d0/(w*sqrt(3.0d0))
        rc  = c*dt
        rc2 = rc**2
        chi1d2 = chi1**2 + chi3**2
        x21 = rc2 - chi1d2
        if(x21.le.0.) then
          p111  = 0.0d0
	  e111  = 0.0d0
	  r111b = 0.0d0
        else
          c2a = chi2 + 0.50d0
          c2b = chi2 - 0.50d0
          x1 = sqrt(x21)
          if(abs(c2a).ge.x1.and.abs(c2b).ge.x1) then
            if((c2a*c2b).lt.0.) then
              i31 = chi3/chi1d2
              k1  = rc2/chi1d2
	      n11 = x1/chi1d2
              fb1 = 8.0d0*i31/3.0d0
              p111  = -x1*fb1*(1.0d0 - k1)
              r111b = chi1**2*fb1*n11*(1.0d0 - 4.0d0*k1)
	      e111  = -4.0d0*i31*x1
            else
              p111  = 0.0d0
	      r111b = 0.0d0
	      e111  = 0.0d0
            endif
          else if(abs(c2a).le.x1.and.abs(c2b).le.x1) then
            i31 = chi3/chi1d2
            k1  = rc2/chi1d2

            r2 = chi1**2 + c2a**2 + chi3**2
            r  = sqrt(r2)
            j1 = 1.0d0/r2 + 2.0d0/chi1d2
	    fb5 = c2a/r
            l2 = rc*fb5
	    m  = rc/r
            fb1 = l2*(1.0d0 - rc2*j1/3.0d0)
            fb3 = l2*(j1*(1.0d0 - 4.0d0/3.0d0*k1) - (m/r)**2)
            r2 = chi1**2 + c2b**2 + chi3**2
            r  = sqrt(r2)
            j1 = 1.0d0/r2 + 2.0d0/chi1d2
	    fb6 = c2b/r
            l2 = rc*fb6
	    m  = rc/r
            fb2 = l2*(1.0d0 - rc2*j1/3.0d0)
            fb4 = l2*(j1*(1.0d0 - 4.0d0/3.0d0*k1) - (m/r)**2)

	    const = 2.0d0*i31
            p111 = -const*(fb1 - fb2)
	    r111b = const*chi1**2*(fb3 - fb4)
	    e111 = -const*rc*(fb5 - fb6)
          else if(c2a.gt.x1.and.c2b.lt.x1) then
            i31 = chi3/chi1d2
            k1  = rc2/chi1d2
	    n11 = x1/chi1d2
            fb1 = 2.0d0*x1*(1.0d0 - k1)/3.0d0
            fb3 = 2.0d0/3.0d0*n11*(1.0d0 - 4.0d0*k1)

            r2 = chi1**2 + c2b**2 + chi3**2
            r  = sqrt(r2)
            j1 = 1.0d0/r2 + 2.0d0/chi1d2
            l2 = rc*c2b/r
            m  = rc/r
            fb2 = l2*(1.0d0 - rc2*j1/3.0d0)
            fb4 = l2*(j1*(1.0d0 - 4.0d0/3.0d0*k1) - (m/r)**2)

	    const = 2.0d0*i31
            p111 = const*(-fb1 + fb2)
	    r111b = const*chi1**2*(fb3 - fb4)
            e111 = const*(-x1 + l2)
          else if(c2a.gt.(-x1).and.c2b.lt.(-x1)) then
            i31 = chi3/chi1d2
            k1  = rc2/chi1d2
	    n11 = x1/chi1d2
            fb1 = 2.0d0*x1*(1.0d0 - k1)/3.0d0
            fb3 = 2.0d0/3.0d0*n11*(1.0d0 - 4.0d0*k1)

            r2 = chi1**2 + c2a**2 + chi3**2
            r  = sqrt(r2)
            j1 = 1.0d0/r2 + 2.0d0/chi1d2
            l2 = rc*c2a/r
	    m  = rc/r
            fb2 = l2*(1.0d0 - rc2*j1/3.0d0)
            fb4 = l2*(j1*(1.0d0 - 4.0d0/3.0d0*k1) - (m/r)**2)

	    const = 2.0d0*i31
            p111 = -const*(fb1 + fb2)
            r111b = const*chi1**2*(fb3 + fb4)
            e111 = -const*(x1 + l2)
          else
            pause 'bad condition'
          endif
	  p111 = p111 + r111b - e111
        endif

        return
        END
c
        REAL(8) FUNCTION q111(chi1, chi2, chi3, dt, w)
	IMPLICIT NONE
        REAL(8) chi1, chi2, chi3, dt, w
        REAL(8) c2a, c2b, r111a
        REAL(8) c, p, p2, p3, rc2, rc, x21, x1, r2, r
        REAL(8) chi1d2, fb1, fb2, fb3, fb4
        REAL(8) i31, j1, k1, l2, m, n11

        if(chi3.eq.0.) then
          q111 = 0.0d0
	  r111a = 0.0d0
          return
        endif

        c = 1.0d0/w
        p = 1.0d0/sqrt(3.0d0)
        p2 = 1.0d0/3.0d0
	p3 = p/3.0d0
        rc  = c*dt
        rc2 = rc**2
        chi1d2 = chi1**2 + chi3**2
        x21 = rc2 - chi1d2
        if(x21.le.0.) then
          q111  = 0.0d0
	  r111a = 0.0d0
        else
          c2a = chi2 + 0.50d0
          c2b = chi2 - 0.50d0
          x1  = sqrt(x21)
          if(abs(c2a).ge.x1.and.abs(c2b).ge.x1) then
            if((c2a*c2b).lt.0.) then
              i31 = chi3/chi1d2
              k1  = rc2/chi1d2
	      n11 = x1/chi1d2
              fb1 = (3.0d0 - 4.0d0*p2 - 2.0d0*p2*k1)/3.0d0
              fb3 = i31*n11*(1.0d0 - 4.0d0*k1)
              q111 = 4.0d0*p*i31*fb1
              r111a = p3*(8.0d0/3.0d0)*chi1**2*fb3
            else
              q111  = 0.0d0
	      r111a = 0.0d0
            endif
          else if(abs(c2a).le.x1.and.abs(c2b).le.x1) then
            i31 = chi3/chi1d2
            k1  = rc2/chi1d2

            r2 = chi1**2 + c2a**2 + chi3**2
            r  = sqrt(r2)
            j1 = 1.0d0/r2 + 2.0d0/chi1d2
            l2 = rc*c2a/r
            m  = rc/r
            fb1 = l2*(1.0d0 - p2 - p2*rc2*j1/3.0d0)
            fb3 = l2*(j1*(1.0d0 - 4.0d0/3.0d0*k1) - (m/r)**2)
            r2 = chi1**2 + c2b**2 + chi3**2
            r  = sqrt(r2)
            j1 = 1.0d0/r2 + 2.0d0/chi1d2
            l2 = rc*c2b/r
	    m  = rc/r
            fb2 = l2*(1.0d0 - p2 - p2*rc2*j1/3.0d0)
            fb4 = l2*(j1*(1.0d0 - 4.0d0/3.0d0*k1) - (m/r)**2)

            q111 = 2.0d0*p*i31*(fb1 - fb2)
            r111a = 2.0d0*p3*i31*chi1**2*(fb3 - fb4)

          else if(c2a.ge.x1.and.c2b.lt.x1) then
            i31 = chi3/chi1d2
            k1  = rc2/chi1d2
	    n11 = x1/chi1d2
            fb1 = x1*(3.0d0 - 4.0d0*p2 - 2.0d0*p2*k1)/3.0d0
            fb3 = 2.0d0/3.0d0*n11*(1.0d0 - 4.0d0*k1)

            r2 = chi1**2 + c2b**2 + chi3**2
            r  = sqrt(r2)
            j1 = 1.0d0/r2 + 2.0d0/chi1d2
            l2 = rc*c2b/r
	    m  = rc/r
            fb2 = l2*(1.0d0 - p2 - p2*rc2*j1/3.0d0)
            fb4 = l2*(j1*(1.0d0 - 4.0d0/3.0d0*k1) - (m/r)**2)

            q111  = 2.0d0*p*i31*(fb1 - fb2)
            r111a = 2.0d0*p3*i31*chi1**2*(fb3 - fb4)
          else if(c2a.ge.(-x1).and.c2b.lt.(-x1)) then
            i31 = chi3/chi1d2
            k1  = rc2/chi1d2
	    n11 = x1/chi1d2
            fb1 = x1*(3.0d0 - 4.0d0*p2 - 2.0d0*p2*k1)/3.0d0
            fb3 = 2.0d0/3.0d0*n11*(1.0d0 - 4.0d0*k1)

            r2 = chi1**2 + c2a**2 + chi3**2
            r  = sqrt(r2)
            j1 = 1.0d0/r2 + 2.0d0/chi1d2
            l2 = rc*c2a/r
	    m  = rc/r
            fb2 = l2*(1.0d0 - p2 - p2*rc2*j1/3.0d0)
            fb4 = l2*(j1*(1.0d0 - 4.0d0/3.0d0*k1) - (m/r)**2)

            q111  = 2.0d0*p*i31*(fb1 + fb2)
            r111a = 2.0d0*p3*i31*chi1**2*(fb3 + fb4)
          else
            pause 'bad if'
          endif
	  q111 = q111 - r111a
        endif

        return
        END
c
