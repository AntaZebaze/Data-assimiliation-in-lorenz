!
module ad_lorenz98nl_mod
    !
    use constants
    use modern_tools
    use freerun_mod
    !
implicit none
    !
    private
    public :: ad_lorenz98nl
    !
contains
    !
    subroutine ad_lorenz98nl(bg_fname, ad_x, ad_xt, imp)
        !c
        !c declare static variables
        !c
        !c
        !c declare passed variables
        !c
        !integer(int_kind), intent(in) :: icycle, iter
        character(*), intent(in) :: bg_fname
        real(real_kind),dimension(:), intent(in out) :: ad_x
        real(real_kind),dimension(:,0:), intent(out) :: ad_xt
        real(real_kind),dimension(:,0:), intent(in) :: imp
        !c
        !c declare local variables
        !c
        real(real_kind),dimension(size(ad_x))      :: x, g1, g2, g3!, g4
        integer(int_kind) :: i
        real(real_kind),dimension(size(ad_x)) :: ad_k1,ad_k2,ad_k3,ad_k4,ad_g
        real(real_kind),dimension(size(ad_x)) :: ad_g1,ad_g2,ad_g3,tend
        !c
        !c background states
        !c
        !c      real(real_kind),dimension(NX,0:4*NT) :: xt
        real(real_kind),dimension(size(ad_xt,1),0:ubound(ad_xt,2)) :: xt  
        !c
        integer :: nx, nt
        !
        nx = size(ad_xt,1)
        nt = ubound(ad_xt,2)
        !c get background states
        !c
        !
        !c
        !c begin computations
        !c
        !c  if(icycle .eq. 1 .and. iter .eq. 1) then
        !bg_fname = make_cycle_fname('DATA/ubg',icycle)
        !write(*,*)"in, ad_lorenz98nl calling read_binary_data, with shape ", shape(xt)
        call read_binary_data(xt, bg_fname)

        !c
        !c begin computations
        !c
        ad_k4 = 0.0
        ad_k3 = 0.0
        ad_k2 = 0.0
        ad_k1 = 0.0
        ad_g3 = 0.0
        ad_g2 = 0.0
        ad_g1 = 0.0
        ad_g  = 0.0
        ad_xt = 0.0
        ad_x = imp(:,NT)
        do i = NT, 1, -1
            !! adjoint calculations

            x = xt(:,i-1)
            call nl_lorenz98_sub(x,tend)
            g1 = x + p5*dt*tend
            call nl_lorenz98_sub(g1,tend)
            g2 = x + p5*dt*tend
            call nl_lorenz98_sub(g2,tend)
            g3 = x + dt*tend


            ad_k4 = ad_k4 + ad_x/6.0
            ad_k3 = ad_k3 + 2.0*ad_x/6.0
            ad_k2 = ad_k2 + 2.0*ad_x/6.0
            ad_k1 = ad_k1 + ad_x/6.0

            ad_g  = ad_g + dt*ad_k4
            ad_k4 = 0.0
            call ad_lorenz98_sub(g3,ad_g3,ad_g)

            ad_k3 = ad_k3 + ad_g3
            ad_x = ad_x  + ad_g3
            ad_g3 = 0.0
            ad_g  = ad_g + dt*ad_k3
            ad_k3 = 0.0
            call ad_lorenz98_sub(g2,ad_g2,ad_g)

            ad_k2 = ad_k2 + p5*ad_g2
            ad_x  = ad_x  + ad_g2
            ad_g2 = 0.0
            ad_g  = ad_g + dt*ad_k2
            ad_k2 = 0.0
            call ad_lorenz98_sub(g1,ad_g1,ad_g)

            ad_k1 = ad_k1 + p5*ad_g1
            ad_x  = ad_x + ad_g1
            ad_g1 = 0.0
            ad_g  = ad_g + dt*ad_k1
            ad_k1 = 0.0
            call ad_lorenz98_sub(x,ad_x,ad_g)
            ad_x = ad_x + imp(:,i-1)
            ad_xt(:,i-1) = ad_x
            !c   write(*,*) 'adjoint max/min',i,maxval(ad_x),minval(ad_x)
        enddo
        ! just added
        ad_xt(:,0) = ad_xt(:,0) + ad_x
        !c
    end subroutine ad_lorenz98nl
    !c
    !c lorenz98_sub is next
    !c
    subroutine ad_lorenz98_sub( y, ad_y, ad_tend )
        !c
        !c declare passed variables
        !
        real(real_kind),dimension(:), intent(in) :: y
        real(real_kind),dimension(:), intent(in out) :: ad_y,ad_tend
        !c
        !c declare local variables
        !c
        integer(int_kind) :: i,j,l,l1,l2,l3,l4
        real(real_kind)   :: ad_ts
        !c
        integer :: nx
        !
        nx = size(y)
        !c
        !c calculate terms
        !c
        ad_ts = 0.0
        do l = 1 , NX
            ad_y(l) = ad_y(l) - ad_tend(l)
            ad_ts   = ad_ts + ad_tend(l)
            ad_tend(l) = 0.0
            do j = -NK/2, NK/2
                do i = -NK/2,NK/2
                    l1=l-2*NK-i
                    if (l1 .lt. 1) l1 = NX +l1-1
                    l2=l-NK-j
                    if (l2 .lt. 1) l2 = NX+l2-1
                    l3 = l-NK+j-i
                    if (l3 .lt. 1) l3 = NX +l3-1
                    l4 = l+NK+j
                    if(l4 .gt. NX) l4 = l4-NX
                    ad_y(l2) = ad_y(l2) - (y(l1)/(NK**2.))*ad_ts
                    ad_y(l1) = ad_y(l1) - (y(l2)/(NK**2.))*ad_ts
                    ad_y(l4) = ad_y(l4) + (y(l3)/(NK**2.))*ad_ts
                    ad_y(l3) = ad_y(l3) + (y(l4)/(NK**2.))*ad_ts
                    ad_ts    = ad_ts
                enddo
            enddo            
            ad_ts = 0.0
        enddo
    end subroutine ad_lorenz98_sub
    !
end module ad_lorenz98nl_mod