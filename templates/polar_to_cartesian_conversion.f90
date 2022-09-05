subroutine xy(r,th,x,y)
    implicit none
    real*8 x,y,r,th
    x=r*cos(th)
    y=r*sin(th)
end

subroutine rt(r,th,x,y)
    implicit none
    real*8 r,th,x,y,pi

    pi=dacos(-1.d0)
    r=sqrt(x**2 + y**2)
    th=datan(y/x)
    ! th=datan(abs(y/x))
    
    if(x<0 .and. y>0)then
        th=th+pi/2
    else if(x<0 .and. y<0)then
        th=th+pi
        
    else if(x>0 .and. y<0)then
        th=th+pi*3.d0/2
    endif

end

program main
    implicit none
    real*8 x,y,r,th,pi

    x=-2.d0
    y=1.d0
    pi = dacos(-1.d0)
    
    call rt(r,th,x,y)
    call xy(r,th,x,y)

    write(6,*)r,x,y
    

end




    




