        DOUBLE PRECISION FUNCTION ker22s(c1, c2, c3, t, w)
        DOUBLE PRECISION c1, c2, c3, t, w, t1
        DOUBLE PRECISION kernel1, kernel2
        DOUBLE PRECISION comp1, comp2
        DOUBLE PRECISION p221, q221
        EXTERNAL p221, q221

        t1 = t + 1.0d0
              kernel1 = p221(c1+0.50d0, c2, c3, t1, w)
     &                - p221(c1-0.50d0, c2, c3, t1, w)
              kernel2 = p221(c1+0.50d0, c2, c3, t, w)
     &                - p221(c1-0.50d0, c2, c3, t, w)
              comp1 = kernel1 - kernel2

              kernel1 = q221(c1+0.50d0, c2, c3, t1, w)
     &                - q221(c1-0.50d0, c2, c3, t1, w)
              kernel2 = q221(c1+0.50d0, c2, c3, t, w)
     &                - q221(c1-0.50d0, c2, c3, t, w)
              comp2 = kernel1 - kernel2

              ker22s = comp1 + comp2 

        return
        END
c
        DOUBLE PRECISION FUNCTION p221(chi1, chi2, chi3, dt, w)
        DOUBLE PRECISION chi1, chi2, chi3, dt, w
        DOUBLE PRECISION c2a, c2b, r221b
        DOUBLE PRECISION c, rc2, rc, x21, x1, r2, r
        DOUBLE PRECISION chi1d2, fb1, fb2, fb3, fb4
        DOUBLE PRECISION i31, j1, k1, l2, m2, const

        if(chi3.eq.0.) then
          p221  = 0.0d0
          return
        endif

        c = 1.0d0/(w*sqrt(3.0d0))
        rc  = c*dt
        rc2 = rc**2.
        chi1d2 = chi1**2. + chi3**2.
        x21 = rc2 - chi1d2
        if(x21.le.0.) then
          p221  = 0.0d0
        else
          c2a = chi2 + 0.50d0
          c2b = chi2 - 0.50d0
          x1 = sqrt(x21)
          if(abs(c2a).ge.x1.and.abs(c2b).ge.x1) then
	    p221  = 0.0d0
          else if(abs(c2a).le.x1.and.abs(c2b).le.x1) then
            i31 = chi3/chi1d2
            k1  = rc2/chi1d2

            r2 = chi1**2. + c2a**2. + chi3**2.
            r  = sqrt(r2)
            j1 = 1.0d0/r2 + 2.0d0/chi1d2
            l2 = rc*c2a/r
	    m2 = rc2/r2
            fb1 = l2*(1.0d0 - rc2*j1/3.0d0)
            fb3 = l2*c2a**2./r2*(1.0d0 - 2.0d0/3.0d0*k1 - m2)
            r2 = chi1**2. + c2b**2. + chi3**2.
            r  = sqrt(r2)
            j1 = 1.0d0/r2 + 2.0d0/chi1d2
            l2 = rc*c2b/r
	    m2 = rc2/r2
            fb2 = l2*(1.0d0 - rc2*j1/3.0d0)
            fb4 = l2*c2b**2./r2*(1.0d0 - 2.0d0/3.0d0*k1 - m2)

	    const = 2.0d0*i31
            p221 = -fb1 + fb2
            r221b = fb3 - fb4
	    p221 = const*(p221 + r221b)
          else if(c2a.gt.x1.and.c2b.lt.x1) then
            i31 = chi3/chi1d2
            k1  = rc2/chi1d2
            fb1 = 2.0d0*x1*(1.0d0 - k1)/3.0d0

            r2 = chi1**2. + c2b**2. + chi3**2.
            r  = sqrt(r2)
            j1 = 1.0d0/r2 + 2.0d0/chi1d2
            l2 = rc*c2b/r
            m2  = rc2/r2
            fb2 = l2*(1.0d0 - rc2*j1/3.0d0)
            fb4 = l2*c2b**2./r2*(1.0d0 - 2.0d0/3.0d0*k1 - m2)

	    const = 2.0d0*i31
            p221 = const*(fb2 - fb4)
          else if(c2a.gt.(-x1).and.c2b.lt.(-x1)) then
            i31 = chi3/chi1d2
            k1  = rc2/chi1d2
            fb1 = 2.0d0*x1*(1.0d0 - k1)/3.0d0

            r2 = chi1**2. + c2a**2. + chi3**2.
            r  = sqrt(r2)
            j1 = 1.0d0/r2 + 2.0d0/chi1d2
            l2 = rc*c2a/r
	    m2 = rc2/r2
            fb2 = l2*(1.0d0 - rc2*j1/3.0d0)
            fb4 = l2*c2a**2./r2*(1.0d0 - 2.0d0/3.0d0*k1 - m2)
	    const = 2.0d0*i31
            p221 = const*(-fb2 + fb4)
          else
            pause 'bad condition'
          endif
        endif

        return
        END
c
        DOUBLE PRECISION FUNCTION q221(chi1, chi2, chi3, dt, w)
        DOUBLE PRECISION chi1, chi2, chi3, dt, w
        DOUBLE PRECISION c2a, c2b, r221a
        DOUBLE PRECISION c, p, p2, rc2, rc, x21, x1, r2, r
        DOUBLE PRECISION chi1d2, fb1, fb2, fb3, fb4
        DOUBLE PRECISION i31, j1, k1, l2, m2, const

        if(chi3.eq.0.) then
          q221 = 0.0d0
          return
        endif

        c = 1.0d0/w
        p = 1.0d0/sqrt(3.0d0)
        p2 = 1.0d0/3.0d0
        rc  = c*dt
        rc2 = rc**2.
        chi1d2 = chi1**2. + chi3**2.
        x21 = rc**2. - chi1d2
        if(x21.le.0.) then
          q221  = 0.0d0
        else
          c2a = chi2 + 0.50d0
          c2b = chi2 - 0.50d0
          x1  = sqrt(x21)
          if(abs(c2a).ge.x1.and.abs(c2b).ge.x1) then
            if((c2a*c2b).lt.0.) then
              i31 = chi3/chi1d2
              k1  = rc2/chi1d2
              fb1 = 4.0d0*p*i31*x1/3.0d0
              q221 = fb1*(3.0d0 - 4.0d0*p2 - 2.0d0*p2*k1)
              r221a = 2.0d0*p2*fb1*(1.0d0 - k1)
            else
              q221  = 0.0d0
	      r221a = 0.0d0
            endif
          else if(abs(c2a).le.x1.and.abs(c2b).le.x1) then
            i31 = chi3/chi1d2
            k1  = rc2/chi1d2

            r2 = chi1**2. + c2a**2. + chi3**2.
            r  = sqrt(r2)
            j1 = 1.0d0/r2 + 2.0d0/chi1d2
            l2 = rc*c2a/r
            m2 = rc2/r2
            fb1 = l2*(1.0d0 - p2 - p2*rc2*j1/3.0d0)
            fb3 = l2*c2a**2./r2*(1.0d0 - 2.0d0/3.0d0*k1 - m2)
            r2 = chi1**2. + c2b**2. + chi3**2.
            r  = sqrt(r2)
            j1 = 1.0d0/r2 + 2.0d0/chi1d2
            l2 = rc*c2b/r
	    m2 = rc2/r2
            fb2 = l2*(1.0d0 - p2 - p2*rc2*j1/3.0d0)
            fb4 = l2*c2b**2./r2*(1.0d0 - 2.0d0/3.0d0*k1 - m2)

	    const = 2.0d0*p*i31
            q221 = const*(fb1 - fb2)
            r221a = const*p2*(fb3 - fb4)
          else if(c2a.gt.x1.and.c2b.lt.x1) then
            i31 = chi3/chi1d2
            k1  = rc2/chi1d2
            fb1 = x1*(3.0d0 - 4.0d0*p2 - 2.0d0*p2*k1)/3.0d0
            fb3 = 2.0d0*x1*(1.0d0 - k1)/3.0d0

            r2 = chi1**2. + c2b**2. + chi3**2.
            r  = sqrt(r2)
            j1 = 1.0d0/r2 + 2.0d0/chi1d2
            l2 = rc*c2b/r
	    m2 = rc2/r2
            fb2 = l2*(1.0d0 - p2 - p2*rc2*j1/3.0d0)
            fb4 = l2*c2b**2./r2*(1.0d0 - 2.0d0/3.0d0*k1 - m2)

	    const = 2.0d0*p*i31
            q221  = const*(fb1 - fb2)
            r221a = const*p2*(fb3 - fb4)
          else if(c2a.gt.(-x1).and.c2b.lt.(-x1)) then
            i31 = chi3/chi1d2
            k1  = rc2/chi1d2
            fb1 = x1*(3.0d0 - 4.0d0*p2 - 2.0d0*p2*k1)/3.0d0
            fb3 = 2.0d0*x1*(1.0d0 - k1)/3.0d0
            r2 = chi1**2. + c2a**2. + chi3**2.
            r  = sqrt(r2)
            j1 = 1.0d0/r2 + 2.0d0/chi1d2
            l2 = rc*c2a/r
	    m2 = rc2/r2
            fb2 = l2*(1.0d0 - p2 - p2*rc2*j1/3.0d0)
            fb4 = l2*c2a**2./r2*(1.0d0 - 2.0d0/3.0d0*k1 - m2)

	    const = 2.0d0*p*i31
            q221  = const*(fb1 + fb2)
            r221a = const*p2*(fb3 + fb4)
          else
            pause 'bad if'
          endif
	  q221 = q221 - r221a
        endif

        return
        END
c
        DOUBLE PRECISION FUNCTION r221b(chi1, chi2, chi3, dt, w)
        DOUBLE PRECISION chi1, chi2, chi3, dt, w
        DOUBLE PRECISION c2a, c2b
        DOUBLE PRECISION c, rc2, rc, x21, x1, r2, r
        DOUBLE PRECISION chi1d2, fb1, fb2
        DOUBLE PRECISION i31, j1, k1, l2, m2

        if(chi3.eq.0.) then
          r221b = 0.0d0
          return
        endif

        c = 1.0d0/(w*sqrt(3.0d0))
        rc  = c*dt
        rc2 = rc**2.
        chi1d2 = chi1**2. + chi3**2.
        x21 = rc2 - chi1d2
        if(x21.le.0.) then
          r221b = 0.0d0
        else
          c2a = chi2 + 0.50d0
          c2b = chi2 - 0.50d0
          x1  = sqrt(x21)
          if(abs(c2a).ge.x1.and.abs(c2b).ge.x1) then
            if((c2a*c2b).lt.0.) then
              i31 = chi3/chi1d2
              k1 = rc2/chi1d2
              r221b = (8.0d0/3.0d0)*i31*x1*(1.0d0 - k1)
            else
              r221b = 0.0d0
            endif
          else if(abs(c2a).le.x1.and.abs(c2b).le.x1) then
            i31 = chi3/chi1d2
             k1 = rc2/chi1d2

             r2 = chi1**2. + c2a**2. + chi3**2.
             r  = sqrt(r2)
             j1 = 1.0d0/r2 + 2.0d0/chi1d2
             l2 = rc*c2a/r
             m2 = rc2/r2
             fb1 = l2*c2a**2./r2*(1.0d0 - 2.0d0/3.0d0*k1 - m2)
             r2 = chi1**2. + c2b**2. + chi3**2.
             r  = sqrt(r2)
             j1 = 1.0d0/r2 + 2.0d0/chi1d2
             l2 = rc*c2b/r
             m2  = rc2/r2
             fb2 = l2*c2b**2./r2*(1.0d0 - 2.0d0/3.0d0*k1 - m2)
             r221b = 2.0d0*i31*(fb1 - fb2)
          else if(c2a.gt.x1.and.c2b.lt.x1) then
            i31 = chi3/chi1d2
            k1  = rc2/chi1d2
            fb1 = 2.0d0*x1*(1.0d0 - k1)/3.0d0

            r2 = chi1**2. + c2b**2. + chi3**2.
            r  = sqrt(r2)
            j1  = 1.0d0/r2 + 2.0d0/chi1d2
            l2  = rc*c2b/r
            m2   = rc2/r2
            fb2 = l2*c2b**2./r2*(1.0d0 - 2.0d0/3.0d0*k1 - m2)
            r221b = 2.0d0*i31*(fb1 - fb2)
          else if(c2a.gt.(-x1).and.c2b.lt.(-x1)) then
            i31 = chi3/chi1d2
            k1  = rc2/chi1d2
            fb1 = 2.0d0*x1*(1.0d0 - k1)/3.0d0

            r2 = chi1**2. + c2a**2. + chi3**2.
            r  = sqrt(r2)
            j1  = 1.0d0/r2 + 2.0d0/chi1d2
            l2  = rc*c2a/r
            m2   = rc2/r2
            fb2 = l2*c2a**2./r2*(1.0d0 - 2.0d0/3.0d0*k1 - m2)
            r221b = 2.0d0*i31*(fb1 + fb2)
          else
            pause 'bad condition'
          endif
        endif

        return
        END
c
