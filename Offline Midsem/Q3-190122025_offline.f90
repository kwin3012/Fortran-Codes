! rk 2.0

subroutine func(z,f,n,t)
    implicit none
    integer n
    real*16 z(n),f(n),t,k1,k3,k2,k4
    
    k1=3.d-12
    k2=1.2d-33
    k3=5.5d-4
    k4=6.9d-16

    f=0.d0
  
     f(1) = -k1*z(1) - k2*z(1)*z(2)*z(3) + k3*z(4) + 2*k4*z(2)*z(4)
     f(2) =  2*k1*z(1) - k2*z(3)*z(2)*z(1) + k3*z(4) - k4*z(2)*z(4)
     f(3) =  0
     f(4) =  k2*z(3)*z(2)*z(1) - k3*z(4) - k4*z(2)*z(4)
     
  
  end
  
  subroutine rk4(z,n,h,ndt,t)
      implicit none
      integer i,n,ndt
      real*16 h,k1(n),k2(n),k3(n),k4(n),k(n),f(n),z(n),t
      
      do i=1,ndt
  
              call func(z,f,n,t)
              k1=f*h
              call func(z+k1/2,f,n,t+h/2)
              k2=f*h
              call func(z+k2/2,f,n,t+h/2)
              k3=f*h
              call func(z+k3,f,n,t+h)
              k4=f*h
  
  
              k=(k1+k2*2+k3*2+k4)/6
              z=z+k
              t=t+h
  
             if(mod(i,10000)==0)then
                 write(6,*) "[O3]",":",z(4)
             endif

      enddo
  
  end
  
  
  program main
      implicit none
      integer,parameter::n=4 ! no of differential equations.
      real*16 h,z(n),t,t0,g
      integer ndt

      g=25.d0/100
      
      z(1)=1.d16*(1+g)   ![O2]
      z(2)=0             ![O]
      z(3)=z(1)/0.22d0   ![M]
      z(4)=0             ![O3]
  
      t0=0.d0
      t=10000000.d0
      h=2.d0

      ndt=(t-t0)/h 
      write(6,*)ndt
      call rk4(z,n,h,ndt,t0)
      
  end
