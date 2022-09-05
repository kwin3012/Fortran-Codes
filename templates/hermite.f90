
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

subroutine func(phi,i,j,her,nb,nx)!i for nb and j for nx
    implicit none
    integer i,j,nb,nx
    real*8 her(nx,nb),phi
    phi = her(j,i)
    
end

subroutine build_hermite(h,nx,nb,xst,xed,dx)
    implicit none
    integer i,j,nx,nb
    real*8 x(nx),h(nx,nb),xst,xed,dx,norm

    x = 0.d0
    do i=1,nx
        x(i) = xst + i*dx
    enddo

    h(:,1)=1.d0
    h(:,2)=2*x(:)

    do i=2,nb-1
        h(:,i+1)=2*x(:)*h(:,i)-2*(i-1)*h(:,i-1)
    enddo  
    
    do i=1,nx
        h(i,:)=h(i,:)*exp(-(x(i)**2)/2)
    enddo

    do i=1,nb
        norm = 0.d0
        do j=1,nx
            norm = norm + h(j,i)**2
        enddo
        norm = norm*dx
        h(:,i)=h(:,i)/sqrt(norm)
    enddo

    

end

program main
    implicit none
    real*8,external::phi
    integer,parameter::ndt=200,nb=30,nx=400
    integer i,j,k
    real*8 t(nb,nb),dx,l,v(nb,nb),dt,xst,xed,pi,x,phii,phij,part(nb,nb),g,h(nb,nb),her(nx,nb)
    complex*16 c(nb),a(nb,nb),zi,ea(nb,nb),zai(nx),xx(0:ndt),norm

    xst=-15.d0
    xed=15.d0
    l=xed-xst
    dx=l/nx
    zi=(0.d0,1.d0)
    dt=0.1d0
    pi=dacos(-1.d0)
    g=25.d0/100

    her=0.d0
    call build_hermite(her,nx,nb,xst,xed,dx)

    ! finding T matrix

   ! t=0.d0
    do i=1,nb
        do j=1,nb
            do k=1,nx
                x=xst+k*dx
                call func(phii,i,k,her,nb,nx)
                call func(phij,j,k,her,nb,nx)
                !t(i,j)=t(i,j)+phii*phij*j*j*pi*pi*(1.d0/2)*(1.d0/l)*((1.d0/l))*dx
                t(i,j)=t(i,j)+phii*phij*dx*(j-1.d0/2)
            enddo
        enddo
    enddo

    ! t=0.d0
    ! do i=1,nb
    !     t(i,i)=i*i*pi*pi*(1.d0/2)*(1.d0/l)*(1.d0/l)
    ! enddo

    ! finding V matrix
    v=0.d0
    do i=1,nb
        do j=1,nb
            do k=1,nx
                x=xst+k*dx
                call func(phii,i,k,her,nb,nx)
                call func(phij,j,k,her,nb,nx)
                v(i,j)=v(i,j)+phii*phij*x*x*(1.d0/2)*dx*g
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
            call func(phij,j,k,her,nb,nx)
            c(j)=c(j)+phij*zai(k)*dx
        enddo
    enddo

    
    ! helper matrix in the last integral to find x(t)
    part=0.d0
    do i=1,nb
        do j=1,nb
            do k=1,nx
                x=xst+k*dx
                call func(phii,i,k,her,nb,nx)
                call func(phij,j,k,her,nb,nx)
                part(i,j)=part(i,j)+phii*phij*x*dx
            enddo
        enddo
    enddo


    ! The time loop in which C Matrix updates itself
    xx=0.d0 ! stores the final answer
    zai = 0.d0


    do i=0,ndt

       zai = 0.d0
       do j=1,nx
            x = xst + j*dx
            do k=1,nb
                call func(phii,k,j,her,nb,nx)
                zai(j) = zai(j) + c(k)*phii
            enddo
        enddo


        norm=0.d0
        do j=1,nx
            x = xst + j*dx
            norm = norm + (conjg(zai(j))*zai(j))*dx
        enddo

        zai = zai/sqrt(norm)

        ! calculating and storing the answer
        ! do j=1,nx
        !     x = xst + j*dx
        !     xx(i) = xx(i) + conjg(zai(j))*x*zai(j)*dx
        ! enddo

        c=c/(sqrt(norm))

        xx(i)=sum(conjg(c)*matmul(h,c))

        c = matmul(ea,c)
        
    enddo

    write(6,*)ndt
    do i=0,ndt
        if(mod(i,10)==0)then
            ! write(14,*)i*dt,real(xx(i)),real(xx1(i))
            write(14,*)i*dt,real(xx(i))
        endif
    enddo
end





