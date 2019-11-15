      subroutine get_vel_visc_terms(vel,beta,visc,dx,lo,hi)
      implicit none
      include 'spec.h'
      real*8 vel(-2:nfine+1)
      real*8 beta(-1:nfine)
      real*8 visc(-1:nfine)
      real*8 dx
      integer lo,hi
      
      integer i
      real*8 beta_lo,beta_hi
      real*8 flux_lo,flux_hi
      real*8 dxsqinv
      
c     Compute D(tau) = d/dx ( a . du/dx ), a=4.mu/3  

      dxsqinv = 1.d0/(dx*dx)
      do i=lo,hi
         if (coef_avg_harm.eq.1) then
            beta_lo = 2.d0 / (1.d0/beta(i)+1.d0/beta(i-1))
            beta_hi = 2.d0 / (1.d0/beta(i)+1.d0/beta(i+1))
         else
            beta_lo = 0.5*(beta(i) + beta(i-1))
            beta_hi = 0.5*(beta(i) + beta(i+1))
         endif

         flux_hi = beta_hi*(vel(i+1) - vel(i  ))
         flux_lo = beta_lo*(vel(i  ) - vel(i-1))
         visc(i) = (flux_hi - flux_lo) * dxsqinv
      enddo
      end
