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
    real*8,parameter:: l=19.d0,pi = dacos(-1.d0),xst=1.d0
    real*8 x,phi
    phi = sqrt(2.d0/l)*sin((i*pi*(x-xst))/l)
end

program main
    implicit none
    real*8,external::phi
    integer,parameter::ndt=200,nb=25,nx=400
    integer i,j,k
    real*8 t(nb,nb),dx,l,v(nb,nb),dt,xst,xed,pi,x,phii,phij,part(nb,nb),xe,xo,vi(nb,nb)
    complex*16 c(nb),a(nb,nb),zi,ea(nb,nb),zai(nx,1),xx(ndt)

    xst=1.d0
    xed=20.d0
    l=xed-xst
    dx=l/nx
    zi=(0.d0,1.d0)
    dt=0.1d0
    pi=dacos(-1.d0)
    xe=3.d0
    xo=12.d0

    ! finding T matrix
    t=0.d0
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

    do i=1,nb
        write(9,*)t(i,:)
    enddo

    ! t=0.d0
    ! do i=1,nb
    !     t(i,i)=i*i*pi*pi*(1.d0/2)*(1.d0/l)*(1.d0/l)
    ! enddo

    ! finding V matrix
    v=0
    do i=1,nb
        do j=1,nb
            v(i,j)=0.d0
            vi(i,j)=0.d0
            do k=1,nx
                x=xst+k*dx
                call func(phii,i,x)
                call func(phij,j,x)
                v(i,j)=v(i,j)+phii*phij*dx*(exp(-2*(x-xe))-2*exp(-(x-xe)))
                ! write(8,*)x,xo
                if(x>xo)then
                    vi(i,j)=vi(i,j)+phii*phij*((x-xo)/(xed-xo))*dx
                endif
            enddo
        enddo
    enddo


    
    

    ! finding A matrix
    a=t+v-zi*vi

    
    

    ! finding exponential of matrix A
    call e_power_a(a,nb,ea)

    do i=1,nb
        write(10,*)ea(i,:)
    enddo

    

    ! zai function at t=0
    zai=0
    do i=1,nx
        x=xst+i*dx
        zai(i,1)=exp(-2*(x-xe)*(x-xe))*(x+1)
    enddo
 
    ! finding C column matrix at t=0
    c=0.d0
    do j=1,nb
        do k=1,nx
            x=xst+k*dx
            call func(phij,j,x)
            c(j)=c(j)+phij*zai(k,1)*dx
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
                part(i,j)=part(i,j)+phii*phij*v(i,j)*dx
            enddo
        enddo
    enddo

    do i=1,nb
        write(7,*)v(i,:)
    enddo

    ! The time loop in which C Matrix updates itself        
    xx=0.d0 ! stores the final answer
    do i=1,ndt
        c = matmul(ea,c)
        xx(i)=sum(conjg(c)*matmul(part,c))
    enddo

    do i=1,ndt
        ! if(i==20 .or. i==100 .or. i==200)then ! at t==2
        !     write(6,*)i*dt,real(xx(i))
        ! endif

        write(6,*)i*dt,(xx(i))

    enddo
end






    























