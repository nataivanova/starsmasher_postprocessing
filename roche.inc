      integer NPMAX
      parameter (NPMAX=200000) ! max number of SPH particiles to work with

c      constants
      real*8 x_PI, umr, c_const, g_const 
      real*8 M_SUN, R_SUN, L_SUN, LM_SUN, LR_SUN, LL_SUN

      parameter (x_PI    = 3.141592653589793238) ! Pi
      parameter (c_const = 2.99792458e10       ) ! speed of light
      parameter (g_const = 6.672607e-8         ) ! gravity const
      parameter (M_SUN   = 1.989e33            ) ! Mass of the Sun
      parameter (L_SUN   = 3.825e33            ) ! Luminosity of the Sun
      parameter (R_SUN   = 6.9599e10           ) ! Radius of the Sun


c   Roche lobe geometry

c   functions
      real*8 rl_eggleton, find_xl1, find_xl_out, find_xl_z
      real*8 phi_x_axis, phi_3d, phi_3d_sph, dphi_x_axis
      real*8 derivative_pot_3d_radial
      integer debug_rl


      real*8 mass1, mass2, a_sep
      real*8 q_rl1, rl_1, xl1_1, xlout_1, xlz_1, error_max_pot
      real*8 q_rl2, rl_2, xl1_2, xlout_2, xlz_2
      common /rlobes/ mass1, mass2, a_sep,
     &               q_rl1, rl_1, xl1_1, xlout_1, xlz_1,
     &               q_rl2, rl_2, xl1_2, xlout_2, xlz_2,
     &               error_max_pot, debug_rl

 




