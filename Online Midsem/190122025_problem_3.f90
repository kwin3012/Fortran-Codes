subroutine func(z,f,n,t)
    real*8 z(n),f(n),t,k1,ki1,k2

    k1=0.1
    ki1=0.051
    k2=0.21

    f(1)=-(k1*z(1)-ki1*z(2))
    f(2)=k1*z(1)-ki1*z(2)-k2*z(2)
    f(3)=k2*z(2)

end

subroutine rk4(z,n,h,ndt)
    implicit none
    integer i,n,ndt
    real*8 h,k1(n),k2(n),k3(n),k4(n),k(n),f(n),z(n),t
    t=0
    
    do i=1,ndt

            call func(z,f,n,t)
            k1=f*h
            call func(z,f+k1/2,n,t+h/2)
            k2=f*h
            call func(z,f+k2/2,n,t+h/2)
            k3=f*h
            call func(z,f+k3,n,t+h)
            k4=f*h

            k=(k1+k2*2+k3*2+k4)/6
            z=z+k
            t=t+h

            if((z(1)-0.5)<10.d-16)then
                goto 10
            endif

    enddo

    10 continue
    write(7,*) "[A]",":",z(1)
    write(7,*) "[I]",":",z(2)

end


program main
    implicit none
    integer,parameter::n=3
    real*8 h,z(n)
    integer t0,t,ndt
    
    z(1)=1![A]
    z(2)=0![B]
    z(3)=0![C]

    t0=0
    t=100
    h=0.1
    ndt=(t-t0)/h +1
    
    call rk4(z,n,h,ndt)
    

end