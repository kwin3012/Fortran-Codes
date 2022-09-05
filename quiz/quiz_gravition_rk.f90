! this is one of the variation that can be seen in the rk.
! rk 2.0

subroutine func(z,f,n,t)
    implicit none
    integer n
    real*16 z(n),f(n),t,G,Me,Ms

    G = 6.67408d-11
    Me = 5.972d24
    Ms = 1.d0

    f=0.d0
  
    f(1) = z(2)
    f(2) = -(z(1)*G*Me*Ms)/(z(1)**2 + z(3)**2)**(3.d0/2)
    f(3) = z(4)
    f(4) = -(z(3)*G*Me*Ms)/(z(1)**2 + z(3)**2)**(3.d0/2)
  
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
  
            
              write(8,*)i*h,z(1)/1000,z(2)/1000
      enddo

    !   write(6,*)i*h,z(1)/1000,z(3)/1000
  
  end
  
  
  program main
      implicit none
      integer,parameter::n=4 ! no of differential equations.
      real*16 h,z(n),t,t0,g,a,m
      integer ndt

      g=25.d0
      a=1+g/100
      
      z(1)=1.d6/a        ![x]
      z(2)=4.07d3        ![Px]
      z(3)=358.d5      ![y]
      z(4)=0             ![py]
  
      t0=0.d0
      t=7200.d0
      h=1.d0

      ndt=(t-t0)/h 
      write(6,*)ndt
      call rk4(z,n,h,ndt,t0)
      
  end
