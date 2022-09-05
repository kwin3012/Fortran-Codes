subroutine fourier_transform(f,g,x,k,n,dx,dk)
    implicit none
    complex*16 f(n),g(n),zi
    integer n,i,j
    real*8 dx,dk,k(n),x(n),pi
    
    zi=(0,1.d0)
    pi=dacos(-1.d0)
    g=0
    do i=1,n
    do j=1,n
    g(i)=g(i)+f(j)*exp(-zi*k(i)*x(j))
    enddo
    g(i)=g(i)*dx
    g(i)=g(i)/(sqrt(2*pi))
    enddo
    
end

subroutine inverse_fourier_tranform(h,g,x,k,n,dx,dk)
    implicit none
    complex*16 h(n),g(n),zi
    integer n,i,j
    real*8 dx,dk,k(n),x(n),pi
    
    zi=(0,1.d0)
    pi=dacos(-1.d0)
    h=0
    
    do i=1,n
    do j=1,n
    h(i)=h(i)+g(j)*exp(zi*k(j)*x(i))
    enddo
    h(i)=h(i)*dk
    h(i)=h(i)/(sqrt(2*pi))
    enddo

end

subroutine split_operator_3(z,zai,dx,dt,nx,xst,s,lam)
    implicit none
    integer nx,i
    real*8 x(nx),k(nx),dx,xst,pi,ks,dk,dt,s,lam
    complex*16 zai(nx),z(nx),f(nx),g(nx),zi

    zi=(0.d0,1.d0)
    pi=dacos(-1.d0)
    ks=-pi/dx
    dk=(2*pi)/(dx*nx)

    do i=1,nx
        x(i) = xst+i*dx
        k(i) = ks+i*dk
        f(i) = zai(i)*exp((s*lam/2)*(1.d0/2)*x(i)*x(i))
    enddo

    call fourier_transform(f,g,x,k,nx,dx,dk)

    do i=1,nx
        g(i)=g(i)*exp(0.5*s*lam*k(i)*k(i))
    enddo

    call inverse_fourier_tranform(z,g,x,k,nx,dx,dk)

    do i=1,nx
        z(i) = z(i)*exp((s*lam/2)*(1.d0/2)*x(i)*x(i))
    enddo
end

subroutine split_operator_5(z,zai,dx,dt,nx,xst)
    implicit none
    integer nx,i
    real*8 x(nx),k(nx),dx,xst,pi,ks,dk,dt,s
    complex*16 zai(nx),z(nx),f(nx),g(nx),zi,lam

    zi=(0.d0,1.d0)
    pi=dacos(-1.d0)
    ks=-pi/dx
    dk=(2*pi)/(dx*nx)
    ! lam = -zi*dt
    lam = -dt
    s = (1.d0/3)*(2 + 1.d0/(2**(1.d0/3)) + (2**(1.d0/3)))
    ! s = 1.d0


    call split_operator_3(z,zai,dx,dt,nx,xst,s,lam)
    zai = z
    z = 0.d0
    call split_operator_3(z,zai,dx,dt,nx,xst,1-2*s,lam)
    zai = z
    z = 0.d0
    call split_operator_3(z,zai,dx,dt,nx,xst,s,lam)

end

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
    integer,parameter::ndt=200,nb=50,nx=400
    integer i,j,k
    real*8 t(nb,nb),dx,l,v(nb,nb),dt,xst,xed,pi,x,phii,phij,part(nb,nb),h(nb,nb),g
    complex*16 c(nb),a(nb,nb),zi,ea(nb,nb),zai(nx),x_zai(0:ndt),norm,norm2,x_c(0:ndt),z(nx)

    xst=-10.d0
    xed=10.d0
    l=xed-xst
    dx=(xed-xst)/nx
    zi=(0.d0,1.d0)
    dt=0.1d0
    pi=dacos(-1.d0)
    g = 25.d0/100

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
    h = t+v
    a= -dt*(h)

    ! finding exponential of matrix A
    ! call e_power_a(a,nb,ea)

    ! zai function at t=0
    zai=0
    do i=1,nx
        x=xst+i*dx
        zai(i)=exp(-0.5*x*x)*(x+1)
    enddo

    ! normalising zai
    norm=0.d0
    do j=1,nx
        x = xst + j*dx
        norm = norm + (conjg(zai(j))*zai(j))*dx
    enddo
    zai = zai/sqrt(norm)
 
    ! finding C column matrix at t=0
    c=0.d0
    do j=1,nb
        do k=1,nx
            x=xst+k*dx
            call func(phij,j,x)
            c(j)=c(j)+phij*zai(k)*dx
        enddo
    enddo

    ! normalizing c
    norm2=0.d0
    do j=1,nb
        norm2 = norm2 + (conjg(c(j))*c(j))
    enddo
    c = c/sqrt(norm2)

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
       
    x_zai = 0.d0 ! stores the final answer with zai column matrix.
    x_c = 0.d0 ! stores final answer with C column matrix.


    !zai = 0.d0

    do i=0,ndt

        norm2=0.d0
        do j=1,nb
            norm2 = norm2 + (conjg(c(j))*c(j))
        enddo

        ! normalizing C column matrix
        c=c/sqrt(norm2)

        ! Also checking the answer via C matrix.
        x_c(i)=sum(conjg(c)*matmul(h,c))

       ! calculating zai at t==0 through already calculated C. 
        ! do j=1,nx
        !     x = xst + j*dx
        !     do k=1,nb
        !         call func(phii,k,x)
        !         zai(j) = zai(j) + c(k)*phii
        !     enddo
        ! enddo

        !calculating normalizing constant
        norm=0.d0
        do j=1,nx
            x = xst + j*dx
            norm = norm + (conjg(zai(j))*zai(j))*dx
        enddo

        zai = zai/sqrt(norm)

        ! calculating and storing the answer
        ! for non matrix use either zai or make it matrix and use c
        ! do j=1,nx
        !     x = xst + j*dx
        !     x_zai(i) = x_zai(i) + conjg(zai(j))*x*zai(j)*dx
        ! enddo

        call split_operator_3(z,zai,dx,dt,nx,xst,1.d0,-dt)
        zai = z

        do j=1,nb
            do k=1,nx
                x=xst+k*dx
                call func(phij,j,x)
                c(j)=c(j)+phij*zai(k)*dx
            enddo
        enddo

        

        

        !!!!!!!!!!!!!!!!!!


        ! norm2=0.d0
        ! do j=1,nb
        !     norm2 = norm2 + (conjg(c(j))*c(j))
        ! enddo

        ! ! normalizing C column matrix
        ! c=c/sqrt(norm2)

        ! ! Also checking the answer via C matrix.
        ! x_c(i)=sum(conjg(c)*matmul(h,c))

        ! making c to c(t+dt) from c(t)
        ! c = matmul(ea,c)
        
    enddo

    do i=0,ndt
        if(mod(i,10)==0)then
            write(12,*)(i)*dt,real(x_c(i))
        endif
    enddo
    
end






    























