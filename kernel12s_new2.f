        DOUBLE PRECISION FUNCTION ker12s(c1, c2, c3, t, w)
        DOUBLE PRECISION c1, c2, c3, t, w, t1
        DOUBLE PRECISION kernel1, kernel2
        DOUBLE PRECISION comp1, comp4, comp5
        DOUBLE PRECISION p112, r112a, r112b
        EXTERNAL p112, r112a, r112b

        t1 = t + 1.0d0
              kernel1 = p112(c1, c2+0.50d0, c3, t1, w)
     &                - p112(c1, c2-0.50d0, c3, t1, w)
              kernel2 = p112(c1, c2+0.50d0, c3, t, w)
     &                - p112(c1, c2-0.50d0, c3, t, w)
              comp1 = kernel1 - kernel2

              kernel1 = r112a(c1+0.50d0, c2, c3, t1, w)
     &                - r112a(c1-0.50d0, c2, c3, t1, w)
              kernel2 = r112a(c1+0.50d0, c2, c3, t, w)
     &                - r112a(c1-0.50d0, c2, c3, t, w)
              comp4 = kernel1 - kernel2

              kernel1 = r112b(c1+0.50d0, c2, c3, t1, w)
     &                - r112b(c1-0.50d0, c2, c3, t1, w)
              kernel2 = r112b(c1+0.50d0, c2, c3, t, w)
     &                - r112b(c1-0.50d0, c2, c3, t, w)
              comp5 = kernel1 - kernel2

              ker12s = comp1 - comp4 + comp5

        return
        END
c
        DOUBLE PRECISION FUNCTION r112b(chi1, chi2, chi3, dt, w)
        DOUBLE PRECISION chi1, chi2, chi3, dt, w
	DOUBLE PRECISION c2a, c2b
        DOUBLE PRECISION c, rc2, rc, x21, r2, r
	DOUBLE PRECISION chi1d2, x1, l1, m2

	if(chi3.eq.0.) then
	  r112b = 0.0d0
	  return
	endif

        c = 1.0d0/(w*sqrt(3.0d0))
        rc  = c*dt
        rc2 = rc**2.
	chi1d2 = chi1**2. + chi3**2.
        x21 = rc2 - chi1d2
	if(x21.le.0.) then
	  r112b = 0.0d0
	else
	  c2a = chi2 + 0.50d0
	  c2b = chi2 - 0.50d0
	  x1 = sqrt(x21)
	  if(abs(c2a).ge.x1.and.abs(c2b).ge.x1) then
	    r112b = 0.0d0
	  else if(abs(c2a).le.x1.and.abs(c2b).le.x1) then
            r2 = chi1d2 + c2a**2.
	    m2 = rc2/r2
            r  = sqrt(r2)
            l1  = rc*chi1/r
            r112b = -2.0d0*chi3/r2*l1*(1.0d0 - m2)

            r2 = chi1d2 + c2b**2.
            m2 = rc2/r2
            r  = sqrt(r2)
            l1  = rc*chi1/r
            r112b = r112b + 2.0d0*chi3/r2*l1*(1.0d0 - m2)
	  else if(c2a.gt.x1.and.c2b.lt.x1) then
            r2 = chi1d2 + c2b**2.
            m2 = rc2/r2
            r  = sqrt(r2)
            l1  = rc*chi1/r
            r112b = 2.0d0*chi3/r2*l1*(1.0d0 - m2)
	  else if(c2a.gt.(-x1).and.c2b.lt.(-x1)) then
            r2 = chi1d2 + c2a**2.
            m2 = rc2/r2
            r  = sqrt(r2)
            l1  = rc*chi1/r
            r112b = -2.0d0*chi3/r2*l1*(1.0d0 - m2)
	  else
	    pause 'bad if'
          endif
        endif

	return
	END
c
        DOUBLE PRECISION FUNCTION r112a(chi1, chi2, chi3, dt, w)
        DOUBLE PRECISION chi1, chi2, chi3, dt, w
	DOUBLE PRECISION c2a, c2b, x1
        DOUBLE PRECISION c, p, p3, rc2, rc, x21, r2, r
	DOUBLE PRECISION chi1d2, l1, m2

	if(chi3.eq.0.) then
	  r112a = 0.0d0
	  return
	endif

        c = 1.0d0/w
        p = 1.0d0/sqrt(3.0d0)
	p3 = p/3.0d0

        rc  = c*dt
        rc2 = rc**2.
	chi1d2 = chi1**2. + chi3**2.
        x21 = rc2 - chi1d2
	if(x21.le.0.) then
	  r112a = 0.0d0
	else
	  c2a = chi2 + 0.50d0
	  c2b = chi2 - 0.50d0
	  x1 = sqrt(x21)
	  if(abs(c2a).ge.x1.and.abs(c2b).ge.x1) then
	    r112a = 0.0d0
	  else if(abs(c2a).le.x1.and.abs(c2b).le.x1) then
            r2 = chi1d2 + c2a**2.
            r  = sqrt(r2)
	    m2 = rc2/r2
	    l1 = rc*chi1/r
            r112a = -chi3/r2*l1*(1.0d0 - m2)

            r2 = chi1d2 + c2b**2. 
            r  = sqrt(r2)
            m2 = rc2/r2
            l1 = rc*chi1/r
            r112a = r112a + chi3/r2*l1*(1.0d0 - m2)
	    r112a = 2.0d0*p3*r112a
	  else if(c2a.gt.x1.and.c2b.lt.x1) then
	    r2 = chi1d2 + c2b**2. 
            r  = sqrt(r2)
            m2 = rc2/r2
            l1 = rc*chi1/r
            r112a = 2.0d0*p3*chi3/r2*l1*(1.0d0 - m2)
	  else if(c2a.gt.(-x1).and.c2b.lt.(-x1)) then
            r2 = chi1d2 + c2a**2. 
            r  = sqrt(r2)
            m2 = rc2/r2
            l1 = rc*chi1/r
            r112a = -2.0d0*p3*chi3/r2*l1*(1.0d0 - m2)
	  else
	    pause 'bad if'
          endif
        endif

        return
        END
c
c
        DOUBLE PRECISION FUNCTION p112(chi1, chi2, chi3, dt, w)
        DOUBLE PRECISION chi1, chi2, chi3, dt, w
	DOUBLE PRECISION c1a, c1b, fb1, fb2
        DOUBLE PRECISION c, rc2, rc, x22, x2, r2, r
        DOUBLE PRECISION chi2d2, i32, l1

	if(chi3.eq.0.) then
	  p112 = 0.0d0
	  return
	endif	
 
        c = 1.0d0/(w*sqrt(3.0d0))
        rc  = c*dt
        rc2 = rc**2.
        chi2d2 = chi2**2. + chi3**2.
        x22 = rc2 - chi2d2
	if(x22.le.0.) then
	  p112 = 0.0d0
	else
	  c1a = chi1 + 0.50d0
	  c1b = chi1 - 0.50d0
	  x2 = sqrt(x22) 
	  if(abs(c1a).ge.x2.and.abs(c1b).ge.x2) then
	    if((c1a*c1b).lt.0.) then
	      i32 = chi3/chi2d2
	      p112 = 2.0d0*i32*x2
	    else
	      p112 = 0.0d0
	    endif
	  else if(abs(c1a).le.x2.and.abs(c1b).le.x2) then
	    i32 = chi3/chi2d2
            r2 = c1a**2. + chi2d2
            r  = sqrt(r2)
	    fb1  = c1a/r

            r2 = c1b**2. + chi2d2
            r  = sqrt(r2)
            fb2  = c1b/r    
	    p112 = i32*rc*(fb1 - fb2)
	  else if(c1a.gt.x2.and.c1b.lt.x2) then
	    i32 = chi3/chi2d2
            r2 = c1b**2. + chi2d2
            r  = sqrt(r2)
            l1  = rc*c1b/r 
	    p112 = i32*(x2 - l1)
	  else if(c1a.gt.(-x2).and.c1b.lt.(-x2)) then
	    i32 = chi3/chi2d2
	    r2 = c1a**2. + chi2d2
            r  = sqrt(r2)
            l1  = rc*c1a/r 
            p112 = i32*(x2 + l1)
	  else
	    pause 'bad if'
          endif
        endif

        return
        END
