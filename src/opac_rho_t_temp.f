      real*8 function opac_rho_t(rho,T)
c this function calculates the opacity for a given rho and T
c the rho and T are converted to T6 = T/10^6 and r = rho/T6^3
      implicit none
      real*8 rho,T,T6,R,logRi,logTi
      integer nummaxR,nummaxT,i
      integer id_T,id_R
      real*8 R1,R2,Q1,Q2,Q3,Q4,x1,x2,y1,y2
      real*8 x,y
      parameter(nummaxR=19,nummaxT=139)
      real*8 logR(nummaxR),logT(nummaxT)
      real*8 logopac(nummaxR,nummaxT)
      common/opac_low/logT,logR,logopac
      T6 = T*1.d-6
      r = rho/(T6**3.d0)
      logRi = log10(r)
      logTi = log10(T)
     
      x = logTi
      y = logRi      
C check that the value of logR is more than -8 otherwise, just give id_r = 1
      if (y .lt. -8.d0) then
          id_r = 1
      else2
C find the value for the given rho
          do i=1,nummaxR-1
             if ((logR(i).le.logRi).and.(logRi.le.logR(i+1))) then
                id_r = i
                exit
             else
                id_r = -1
             endif
          enddo
      endif
C do something similar for logT
      do i=1,nummaxT-1
         if ((logT(i).le.logTi).and.(logTi.le.logT(i+1))) then
            id_T = i
            exit
         else
            id_T = -1
         endif
      enddo
C if id_T and id_r were found, then get its opac value using bilinear interpolation
      if ((id_r.ne.-1).and.(id_t.ne.-1)) then 
          x1 = logT(id_T)
          x2 = logT(id_T+1)
          y1 = logR(id_r)
          y2 = logR(id_r+1)
          Q1 = logopac(id_r,id_T)
          Q2 = logopac(id_r,id_T+1)
          Q3 = logopac(id_r+1,id_T)
          Q4 = logopac(id_r+1,id_T+1)
          R1 = Q2*(x-x1)/(x2-x1)+Q1*(x2-x)/(x2-x1)
          R2 = Q4*(x-x1)/(x2-x1)+Q3*(x2-x)/(x2-x1)
          opac_rho_t = R2*(y-y1)/(y2-y1)+R1*(y2-y)/(y2-y1)

      else
          opac_rho_t=0.d0
      endif
 
      return
      end
