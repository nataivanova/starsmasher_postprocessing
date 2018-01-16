      subroutine getTemperature(q,r,x3)
c     Subroutine to solve 4th order equations to determine the temperature x3
c     for an equation of state with both ideal gas and radiation pressure.
c     Written by Scott Fleming 10/04/02 and James Lombardi 2002-2003
c     The fourth order equation comes from u_gas+ u_rad = u, with
c     u_gas proportional to T and u_rad proportional to T^4
c     In general, we can transform a 4th order equation to x^4+px^2+qx+r=0
c     (see pages 57-58 of Stillwell's "Mathematics and its history" text)
c     but we fortunately don't even have an x^2 term (that is, p=0).
c     Follow Stillwell, we can transform this into a cubic equation:
c     First solve for y by using the fact that B^2-4AC=0
c     equation is then:  y^3=ry+q^2/8
c     using the solution of cubic equations found in Stillwell page 55:
      implicit none
      real*8 q,r,k,B,piece1,piece2
      real*8 y1,y2,yy,AA,b2,c2,x3,kh
c      real*8 piece2old,x3old
      k = 0.125d0*q**2
      kh=0.5d0*k
      if(kh**2-(r/3.d0)**3.le.0.d0)then
         write(*,*) k,r,kh**2-(r/3.d0)**3
         stop 'bad input: imaginary results?'
      endif
      piece1 = kh+(kh**2-(r/3.d0)**3)**0.5d0
c      piece2old = kh-(kh**2-(r/3.d0)**3)**0.5d0
      piece2 = (r/3.d0)**3.d0/piece1
c      write(69,*)piece2old,piece2
      y1 = piece1**(1.d0/3.d0)
c     Fortran can't handle cube roots of neg. #'s
      y2 = -dabs(piece2)**(1.d0/3.d0)
      yy = y1+y2
c     Equation to solve: (x^2+p+y)^2=Ax^2+Bx+C
c     Now take square root of both sides with:
      AA = 2.d0*yy
      B = -q
c      C = -r+y**2
c     Re-writing Ax^2+Bx+C as a square then solving the equation we
c     obtain 2 results:
c     x^2 + (-(A^(1/2)))x + (-B/(2(A)^(1/2))+p+y) = 0 (1)
c     or
c     x^2 + (A^(1/2))x + (B/(2(A)^(1/2))+p+y) = 0     (2)
c     Our solution we're interested in:
      b2 = AA**0.5d0
      c2 = 0.5d0*B/b2 + yy
c     Therefore, we once again have x^2+bx+c=0, and our answer we want is:
c      x3old = 0.5d0*(-b2 + (b2**2-4.d0*c2)**0.5d0)
      x3 = -2.d0*c2/(b2 + (b2**2-4.d0*c2)**0.5d0)
c      write(69,*) 'Temperature Components', x3old, x3
      if(piece1.lt.0.d0) write(*,*)
     $     'piece 1 lt 0',k,r,piece1,piece2
c      if(piece1.eq.-piece2)then
      if(b2.eq.0.d0)then
         write(*,*)
     $        'piece 1 eq -piece 2 (rad pressure dominates)',
     $        k,r,piece1,piece2,b2,c2,x3,(-r)**0.25d0
         x3=(-r -q*( -r-q*(-r)**0.25d0 )**0.25d0)**0.25d0
         write(*,*) x3
         write(*,*) piece1,piece2
      endif
      if(piece2.ge.0.d0) x3=-(r+(r/q)**4)/q 
      end 
