subroutine matrix_inverse(a,u,n)
    implicit none
    integer i,j,n
    real*8 a(n,n),u(n,n),temp1,temp2

    u=0.d0
    do i=1,n
        u(i,i)=1.d0
    enddo

    do i=1,n
        temp1=a(i,i)
        a(i,:)=a(i,:)/temp1
        u(i,:)=u(i,:)/temp1

        do j=1,n
            if(j/=i)then
                temp2=a(j,i)
                a(j,:)=a(j,:)-a(i,:)*temp2
                u(j,:)=u(j,:)-u(i,:)*temp2
            endif
        enddo
    enddo

end


program main
    implicit none
    integer,parameter::n=3
    integer i,j
    real*8 a(n,n),ai(n,n),ans(n,1),col(n,1)

    a(1,:)=(/1,2,6/)
    a(2,:)=(/3,4,1/)
    a(3,:)=(/6,-1,-1/)

    col(1,1)=22
    col(2,1)=26
    col(3,1)=19

    call matrix_inverse(a,ai,n)

    ans = matmul(ai,col)
   
    do i=1,n
        write(6,*)ans(i,1)
    enddo

end
