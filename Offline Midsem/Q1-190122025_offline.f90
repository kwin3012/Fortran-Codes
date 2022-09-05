subroutine e_power_a(a,n,ea)
    implicit none
    integer i,n
    complex*16 a(n,n),ea(n,n),b(n,n),p(n,n)
    
    p=0.d0
    do i=1,n
        p(i,i)=1.d0
    enddo
    
    b=p     
    ea=p
    
    do i=1,100
        ! call mm(b,a,c,n,n,n)
        b = matmul(b,a)
        b = b/i
        ea=ea+b
    enddo

end

subroutine func(phi,i,x)
    implicit none
    integer i
    real*8,parameter:: l=20.d0,pi = dacos(-1.d0),xst=-10.d0
    real*8 x,phi
    phi = sqrt(2.d0/l)*sin((i*pi*(x-xst))/l)
end

program main
    implicit none
    real*8,external::phi
    integer,parameter::ndt=200,nb=50,nx=400
    integer i,j,k
    real*8 t(nb,nb),dx,l,v(nb,nb),dt,xst,xed,pi,x,phii,phij,part(nb,nb),g,h(nb,nb)
    complex*16 c(nb),a(nb,nb),zi,ea(nb,nb),zai(nx),xx(0:ndt),nor,mm(1,1),xx1(0:ndt)

    xst=-10.d0
    xed=10.d0
    l=xed-xst
    dx=l/nx
    zi=(0.d0,1.d0)
    dt=0.1d0
    pi=dacos(-1.d0)
    g=25.d0/100

    ! finding T matrix
    t=0.d0
    do i=1,nb
        do j=1,nb
            do k=1,nx
                x=xst+k*dx
                call func(phii,i,x)
                call func(phij,j,x)
                t(i,j)=t(i,j)+phii*phij*j*j*pi*pi*(1.d0/2)*(1.d0/l)*((1.d0/l))*dx
            enddo
        enddo
    enddo

    ! t=0.d0
    ! do i=1,nb
    !     t(i,i)=i*i*pi*pi*(1.d0/2)*(1.d0/l)*(1.d0/l)
    ! enddo

    ! finding V matrix
    v=0
    do i=1,nb
        do j=1,nb
            v(i,j)=0
            do k=1,nx
                x=xst+k*dx
                call func(phii,i,x)
                call func(phij,j,x)
                v(i,j)=v(i,j)+phii*phij*x*x*(1.d0/2)*dx*(1.d0+g)
            enddo
        enddo
    enddo

    ! finding A matrix
    h=t+v
    a=-dt*(h)
    

    ! finding exponential of matrix A
    call e_power_a(a,nb,ea)

    ! zai function at t=0
    zai=0.d0
    do i=1,nx
        x=xst+i*dx
        zai(i)=exp(-0.5*x*x)*(x+1)
    enddo
 
    ! finding C column matrix at t=0
    c=0.d0
    do j=1,nb
        do k=1,nx
            x=xst+k*dx
            call func(phij,j,x)
            c(j)=c(j)+phij*zai(k)*dx
        enddo
    enddo

    ! helper matrix in the last integral to find x(t)
    part=0.d0
    do i=1,nb
        do j=1,nb
            do k=1,nx
                x=xst+k*dx
                call func(phii,i,x)
                call func(phij,j,x)
                part(i,j)=part(i,j)+phii*phij*x*dx
            enddo
        enddo
    enddo

    ! The time loop in which C Matrix updates itself        
    xx=0.d0 ! stores the final answer
    xx1=0.d0
    do i=0,ndt
        nor=0.d0
        do j=1,nb
            nor = nor + (conjg(c(j))*c(j))
        enddo
        c=c/sqrt(nor)
        xx(i)=sum(conjg(c)*matmul(h,c))
        xx1(i)=sum(conjg(c)*matmul(matmul(h,h),c))
        c = matmul(ea,c)
    enddo

    write(10,*)ndt
    do i=0,ndt
        if(mod(i,10)==0)then
            write(15,*)i*dt,real(xx(i)),real(xx1(i))
        endif
    enddo
end






    























