real*8 function f(x,y)
real*8 x,y
f=x+2*y
end

program euler
    implicit none
    real*8,external::f
    real*8 h,x,y
    integer n,x0,y0,xt,i
    
    x0=0
    y0=0
    xt=1
    h=0.1
    n=(xt-x0)/h
    x=x0
    y=y0
    
    do i=1,n
    y=y+h*f(x,y)
    x=x+h
    enddo
    
    print*,x,y
    
end
    
    
    
    
    
    
    

