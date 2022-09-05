subroutine mm(a,b,r,x,y,z)
    implicit none
    integer i,j,k,x,y,z
    real*8 a(x,y),b(y,z),r(x,z)
    
    do i=1,x
    do j=1,z
    r(i,j)=0
    do k=1,y
    r(i,j)=r(i,j)+a(i,k)*b(k,j)
    enddo
    enddo
    enddo
    
end
    
subroutine e_power_a(a,n,ea)
    implicit none
    integer i,n
    real*8 a(n,n),ea(n,n),b(n,n),p(n,n),c(n,n)
    
    do i=1,n
        p(i,i)=1
    enddo
        
    b=p     
    ea=p
    
    do i=1,100
    ! call mm(b,a,c,n,n,n)
    c = matmul(b,a)
    c=c/i
    ea=ea+c
    b=c
    enddo

end

program main
    implicit none
    integer,parameter::n=2
    real*8 a(n,n),ea(n,n)
    
    a(1,1)=1
    a(1,2)=0
    a(2,1)=0
    a(2,2)=1
    
    call e_power_a(a,n,ea)
    
    print*,ea
    
end
    
        
        
