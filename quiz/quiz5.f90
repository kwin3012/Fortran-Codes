subroutine fourier_transform(f,g,x,k,n,dx,dk)
    implicit none
    complex*16 f(n),g(n),zi
    integer n,i,j
    real*8 dx,dk,k(n),x(n),pi
    
    zi=(0,1.d0)
    pi=dacos(-1.d0)
    g=0
    
    do i=1,n
    do j=1,n
    g(i)=g(i)+f(j)*exp(-zi*k(i)*x(j))
    enddo
    g(i)=g(i)*dx
    g(i)=g(i)/(sqrt(2*pi))
    enddo
    
end

subroutine inverse_fourier_tranform(h,g,x,k,n,dx,dk)
    implicit none
    complex*16 h(n),g(n),zi
    integer n,i,j
    real*8 dx,dk,k(n),x(n),pi
    
    zi=(0,1.d0)
    pi=dacos(-1.d0)
    h=0
    
    do i=1,n
    do j=1,n
    h(i)=h(i)+g(j)*exp(zi*k(j)*x(i))
    enddo
    h(i)=h(i)*dk
    h(i)=h(i)/(sqrt(2*pi))
    enddo

end

program main
    implicit none
    integer,parameter::n=200
    complex*16 f(n),g(n),h(n),zi,s
    integer i,j,l
    real*8 dx,dk,k(n),x(n),xs,ks,pi,ff(n),dt
    
    xs=-10.d0
    dx=0.1
    pi=dacos(-1.d0)
    ks=-pi/dx
    dk=(2*pi)/(dx*n)
    zi=(0,1.d0)
    dt=0.1

    ! function f(i) x(i) k(i) is been made
    do i=1,n
    x(i)=xs+i*dx
    f(i)=exp(-x(i)*x(i))
    ff(i)=exp(-x(i)*x(i))
    k(i)=ks+i*dk
    enddo

    do j=1,30

        call fourier_transform(f,g,x,k,n,dx,dk)
        
        do i=1,n
        g(i)=g(i)*exp(-0.5*zi*dt*k(i)*k(i))
        enddo

        call inverse_fourier_tranform(h,g,x,k,n,dx,dk)

        f=h
        if(mod(j,10)==0)then
        s=0
        do l=1,n
        s=s+ff(l)*f(l)
        enddo
        s=s*dx
        print*,j/10,s
        endif

    enddo
    
end
    
    
    
            
        
        