subroutine d2p4(f,fd)
    implicit none
    integer,parameter::n=200
    real*8 f(n),fd(n),dx
    integer i
    dx=0.1
    
    do i=3,n-2
        fd(i)=(-f(i-2)*1/12+f(i-1)*4/3-f(i)*5/2+f(i+1)*4/3-f(i+2)*1/12)/dx**2
    enddo

end


program derivative_and_integration
implicit none
    integer,parameter::n=200
    real*8 x(n),f(n),fd(n),dx,s
    integer i
    dx=0.1
    do i=1,n
        x(i)=-10+i*dx
        f(i)=exp(-x(i)**2)
    enddo

    call d2p4(f,fd)

    do i=1,n
        write(6,*)x(i),f(i),fd(i)
    enddo

    s=0
    do i=1,n
        s=s+f(i)
    enddo
    s=s*dx

    write(6,*)s
    write(6,*)sqrt(dacos(-1.d0))

end




  


