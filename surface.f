      program surface
      integer i,j,k,n
      parameter(n=100)
      real*8 phi1
      real*8 x,y,z
      real*8 dx,dy,dz
      dx = 20.d0/n
      dy = 20.d0/n
      dz = 20.d0/n
      x = -10
      y = -10
      z = -10
 
      do i=1,100
        do j=1,100
           do k=1,100
              call phi(phi1,x,y,z,1.d0,1.d0,10.d0)
              write(10,*)x,y,z,phi1
              z = z + dz
           enddo
           y = y + dy
           z = -10
        enddo
        x = x +dx
        y = -10
      enddo

      end
