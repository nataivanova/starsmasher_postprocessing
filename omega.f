      subroutine omega_vec(Ixx,Iyy,Izz,Ixy,Ixz,Iyz,lx,ly,lz,omegax,omegay,omegaz)
      real*8 Ixx,Iyy,Izz,Ixy,Ixz,Iyz,lx,ly,lz,omegax,omegay,omegaz
      omegaz = ((lz-Ixz*lx/Ixx)*(Iyy-Ixy**2/Ixx)-(Iyz-Ixz*Ixy/Ixx)*(ly-Ixy*lx/Ixx))/
     &     ((Izz-Ixz**2/Ixx)*(Iyy-Ixy**2/Ixx)-(Iyz-Ixz*Ixy/Ixx)**2)
      omegay = (ly-Ixy*lx/Ixx)/(Iyy-Ixy**2/Ixx)-(Iyz-Ixz*Ixy/Ixx)/(Iyy-Ixy**2/Ixx)*omegaz
      omegax = lx/Ixx-Ixy/Ixx*omegay-Ixz/Ixx*omegaz
      end
