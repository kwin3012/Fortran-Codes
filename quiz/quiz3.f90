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

subroutine cos_a(a,n,ca)
    implicit none
    integer i,j,n
    real*8 a(n,n),ca(n,n),b(n,n),p(n,n),c(n,n)

    !p is an indentity matrix
    do i=1,n
        p(i,i)=1
    enddo
        
    b=p     
    ca=p
    
    do i=1,100
    ! call mm(b,a,c,n,n,n)
    c=matmul(b,a)
    c=c/i
    b=c
    if(mod(i,2)==0)then
        if(mod(i/2,2)==1)then
            ca=ca-b
        else
            ca=ca+b
        endif
    endif
    enddo

end

program main
    implicit none
    integer,parameter::n=2
    integer i
    real*8 a(n,n),ca(n,n),pi

    pi=dacos(-1.d0)
    
    
    ! a(1,1)=0
    ! a(1,2)=1
    ! a(1,3)=1
    ! a(2,1)=0
    ! a(2,2)=0
    ! a(2,3)=1
    ! a(3,1)=0
    ! a(3,2)=0
    ! a(3,3)=0

    a(1,1)=1
    a(1,2)=2
    a(2,1)=3
    a(2,2)=4
    
    
    call cos_a(a,n,ca)
    

    do i=1,n
        write(6,*) ca(i,1),ca(i,2)
    enddo
        
    
end