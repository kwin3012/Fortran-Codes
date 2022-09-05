subroutine fourier(f,g,x,k,n,dx)
    implicit none
    integer i,j,n
    real*8 x(n),k(n),dx,pi
    complex*16 f(n),g(n),zi

    zi=(0,1.d0)
    pi=dacos(-1.d0)


    do i=1,n    
        g(i)=0
        do j=1,n
            g(i)=g(i)+f(j)*exp(-zi*k(i)*x(j))
        enddo
        g(i)=g(i)*dx*(1.d0/sqrt(2*pi))
    enddo

end

subroutine inverse_fourier(f,g,x,k,n,dk)
    implicit none
    integer i,j,n
    real*8 x(n),k(n),pi,dk
    complex*16 f(n),g(n),zi

    zi=(0,1.d0)
    pi=dacos(-1.d0)

    do i=1,n
        f(i)=0
        do j=1,n
            f(i)=f(i)+g(j)*exp(zi*k(j)*x(i))
        enddo
        f(i)=f(i)*dk*(1.d0/sqrt(2*pi))
    enddo

end

program main
    implicit none
    integer,parameter::n=200
    real*8 ,parameter::pi=dacos(-1.d0)
    integer i,t
    real*8 x(n),k(n),xs,dx,ks,dk,dt
    complex*16 f(n),g(n),zi,h(n),s

    xs=-10.d0
    dx=0.1
    ks=-pi/dx
    dk=2*pi/(dx*n)
    dt=0.1
    zi=(0,1.d0)

    do i=1,n
        x(i)=xs+i*dx
        k(i)=ks+i*dk
        f(i)=exp(-x(i)*x(i))
        h(i)=exp(-x(i)*x(i))
    enddo

do t=1,30

    call fourier(f,g,x,k,n,dx)

    do i=1,n
        g(i)=g(i)*exp(-k(i)*k(i)*zi*dt*(1.d0/2))
    enddo

    call inverse_fourier(f,g,x,k,n,dk)
    ! this f will be initial f(dt)

    if(mod(t,10)==0)then
        ! find out integration of that part
        s=0.d0
        do i=1,n
            s=s+h(i)*f(i)*dx
        enddo
        write(6,*)s
    endif

enddo

end







