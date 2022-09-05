program problem1
    implicit none
    integer,parameter::n=50
    integer i,j,k
    real *8 x,a(3),s1(n),s2(n),ans

    x=0.8d0
    a(1)=20.2d0
    a(2)=30.d0
    a(3)=40.5d0
    
    do k=1,3
        s1(2)=x
        s2(2)=a(k)
        ans=1+x*a(k)
        do i=3,n
            s1(i)=s1(i-1)*(x/(i-1))
            s2(i)=s2(i-1)*(a(k)-(i-2))
            ans=ans+s2(i)*s1(i)
        enddo
        ! write(6,*)s1
        write(6,*)ans
        ! write(6,*)exp(x)
    enddo

end program