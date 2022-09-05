subroutine bubble_sort(a,n)
    implicit none
    integer n,a(n),i,j,temp

    do i=0,n-1
        do j=2,n-i
            if(a(j-1)>a(j))then
                temp=a(j-1)
                a(j-1)=a(j)
                a(j)=temp
            endif
        enddo
    enddo

end



program main
    implicit none
    integer,parameter::n=10
    integer a(n),i

    do i=1,n
        read(5,*)a(i)
    enddo

    call bubble_sort(a,n)
    
    write(6,*)a

end 




