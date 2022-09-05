! program to verify f'' using fourier transform properties

subroutine d2p6(f,n,fd,h)
    implicit none
    integer n,i
    real*8 f(n),fd(n),h
    
    fd=0
    do i=4,n-3
    fd(i)=(f(i-3)*1/90-f(i-2)*3/20+f(i-1)*3/2-f(i)*49/18+f(i+1)*3/2-f(i+2)*3/20+f(i+3)*1/90)/h**2
    enddo
    
end

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
    complex*16 f(n),g(n),h(n)
    integer i
    real*8 dx,dk,k(n),x(n),xs,ks,pi,fd(n),ff(n)
    
    xs=-10.d0
    dx=0.1
    pi=dacos(-1.d0)
    ks=-pi/dx
    dk=(2*pi)/(dx*n)

    ! function f(i) x(i) k(i) is been made
    do i=1,n
    x(i)=xs+i*dx
    f(i)=exp(-x(i)*x(i))
    ff(i)=exp(-x(i)*x(i))
    k(i)=ks+i*dk
    enddo

    call d2p6(ff,n,fd,dx)
    
    call fourier_transform(f,g,x,k,n,dx,dk)

    do i=1,n
    g(i)=-g(i)*k(i)*k(i)
    enddo

    call inverse_fourier_tranform(h,g,x,k,n,dx,dk)
    
    do i=1,n
        write(7,*) x(i),real(h(i)),fd(i)
    enddo
    
end
    
    
    
            
        
        