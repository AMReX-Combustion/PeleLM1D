      subroutine calc_divu(scal,beta,I_R,divu,dx,lo,hi)
      implicit none
      include 'spec.h'

c     Quantities passed in
      real*8 scal(-2:nfine+1,nscal)
      real*8 beta(-1:nfine  ,nscal)
      real*8  I_R(-1:nfine  ,0:Nspec)
      real*8 divu(-1:nfine)
      real*8 dx
      integer lo,hi
      
      real*8 Y(Nspec)
      real*8 HK(Nspec)
      real*8 cpmix,mwmix
      
      real*8 diff(-1:nfine,nscal)
      real*8 diffdiff(-1:nfine)

      real*8 RWRK,rho,T
      integer IWRK,i,n

      real*8 gamma_lo(0:nfine-1,Nspec)
      real*8 gamma_hi(0:nfine-1,Nspec)

c     compute Gamma_m
      call get_spec_visc_terms(scal,beta,diff(:,FirstSpec:),
     $                         gamma_lo,gamma_hi,dx,lo,hi)

c     compute div lambda grad T
      diff(:,Temp) = 0.d0
      call addDivLambdaGradT(scal,beta,diff(:,Temp),dx,lo,hi)

c     compute div h_m Gamma_m
      call get_diffdiff_terms(scal,gamma_lo,gamma_hi,
     $                        diffdiff,dx,lo,hi)

c     combine div lambda grad T + div h_m Gamma_m
      do i=lo,hi
         diff(i,Temp) = diff(i,Temp) + diffdiff(i)
      end do

      do i=lo,hi
         rho = scal(i,Density)
         do n = 1,Nspec
            Y(n) = scal(i,FirstSpec + n - 1) / rho
         enddo
         T = scal(i,Temp)
         call CKMMWY(Y,IWRK,RWRK,mwmix)
         call CKCPBS(T,Y,IWRK,RWRK,cpmix)
         call CKHMS(T,IWRK,RWRK,HK)

         divu(i) = diff(i,Temp)/(rho*cpmix*T)

         do n=1,Nspec
            divu(i) = divu(i)
     &           + (diff(i,FirstSpec+n-1) + I_R(i,n))
     &           *(invmwt(n)*mwmix/rho - HK(n)/(rho*cpmix*T))
         enddo
       enddo

      end
