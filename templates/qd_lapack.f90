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
    real*8,parameter:: l=20.d0,pi = dacos(-1.d0),xs=-10.d0
    real*8 x,phi
    phi = sqrt(2.d0/l)*sin((i*pi*(x-xs))/l)
end

program main
    implicit none
    real*8,external::phi
    integer,parameter::n=200,ndt=100,nb=25,nx=2000
    integer i,j,k
    real*8 t(nb,nb),dx,l,v(nb,nb),dt,xst,xed,pi,x,phii,phij,part(nb,nb)
    complex*16 c(nb),a(nb,nb),zi,ea(nb,nb),zai(nx),xx(ndt),norm,xxx(ndt)

    xst=-10.d0
    xed=10.d0
    l=xed-xst
    dx=(xed-xst)/nx
    zi=(0.d0,1.d0)
    dt=0.1d0
    pi=dacos(-1.d0)

    ! finding T matrix
    t=0
    do i=1,nb
        do j=1,nb
            t(i,j)=0
            do k=1,nx
                x=xst+k*dx
                call func(phii,i,x)
                call func(phij,j,x)
                t(i,j)=t(i,j)+phii*phij*j*j*pi*pi*(1.d0/2)*(1.d0/l)*((1.d0/l))*dx
            enddo
        enddo
    enddo

    ! Another way for T matrix
    t=0.d0
    do i=1,nb
        t(i,i)=i*i*pi*pi*(1.d0/2)*(1.d0/l)*(1.d0/l)
    enddo

    ! finding V matrix
    v=0
    do i=1,nb
        do j=1,nb
            v(i,j)=0
            do k=1,nx
                x=xst+k*dx
                call func(phii,i,x)
                call func(phij,j,x)
                v(i,j)=v(i,j)+phii*phij*x*x*(1.d0/2)*dx
            enddo
        enddo
    enddo

    ! finding A matrix
    a=-zi*dt*(t+v)

    ! finding exponential of matrix A
    call e_power_a(a,nb,ea)

    ! zai function at t=0
    zai=0
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

    ! norm=0.d0
    ! do j=1,nx
    !     x = xst + j*dx
    !     norm = norm + (conjg(zai(j))*zai(j))*dx
    ! enddo
               
    xx=0.d0 ! stores the final answer with zai column matrix.
    xxx=0.d0 ! stores final answer with C column matrix.
    zai = 0.d0

    do i=1,ndt
       ! calculating zai at t==0 
       zai=0.d0
        do j=1,nx
            x = xst + j*dx
            do k=1,nb
                call func(phii,k,x)
                zai(j) = zai(j) +  c(k)*phii
            enddo
        enddo

        !calculating normalizing constant
        norm=0.d0
        do j=1,nx
            x = xst + j*dx
            norm = norm + (conjg(zai(j))*zai(j))*dx
        enddo

        zai = zai/sqrt(norm)

        ! calculating and storing the answer
        do j=1,nx
            x = xst + j*dx
            xx(i) = xx(i) + conjg(zai(j))*x*zai(j)*dx
        enddo

        ! normalizing C column matrix
        c=c/sqrt(norm)

        ! Also checking the answer via C matrix.
        xxx(i)=sum(conjg(c)*matmul(part,c))

        ! making c to c(t+dt) from c(t)
        c = matmul(ea,c)
        
    enddo

    do i=1,ndt
        write(12,*)(i-1)*dt,real(xx(i)),real(xxx(i))
    enddo
    
end






    























