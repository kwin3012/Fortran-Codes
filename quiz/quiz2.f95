subroutine d2p6(f,n,fd,h)
    implicit none
    integer n,i
    real*8 f(n),fd(n),h
    
    fd=0
    do i=4,n-3
    fd(i)=(f(i-3)*1/90-f(i-2)*3/20+f(i-1)*3/2-f(i)*49/18+f(i+1)*3/2-f(i+2)*3/20+f(i+3)*1/90)/h**2
    enddo
    
end

program quiz2
    implicit none
    integer,parameter::n=200
    real*8 h,x(n),f(n),fd(n),s,a
    integer i
    
    a=1.d0+25.d0/100
    h=0.1
    
    do i=1,n
    x(i)=-10+i*h
    f(i)=exp(-a*x(i)*x(i))
    enddo
    
    call d2p6(f,n,fd,h)
    
    print*,fd
    
    s=0
    do i=1,n
    s=s+f(i)*fd(i)
    enddo
    s=s*h
    
    print*,s
    
end
    
    
    
    
    
    
    

