!
module rp_lorenz98nl_mod
    !
    use constants
    use modern_tools
    use freerun_mod
    !
implicit none
    !
    private
    public :: rp_lorenz98nl
    !
contains
    !
    subroutine rp_lorenz98nl(bg_fname, rp_xt, ad_xt)
        !c
        !c declare static variables
        !c
        !c
        !c declare passed variables
        !c
        character(*), intent(in) :: bg_fname
        real(kind=real_kind),dimension(:,0:), intent(in out) :: rp_xt
        real(kind=real_kind),dimension(:,0:), intent(in) :: ad_xt
        !c
        !c declare local variables
        !c
        real(kind=real_kind),dimension(size(rp_xt,1),0:ubound(rp_xt,2)) :: u_bg  
        real(kind=real_kind),dimension(size(rp_xt,1))      :: x,rp_x,tend
        integer(kind=int_kind) :: i
        real(kind=real_kind),dimension(size(rp_xt,1)) :: rp_k1,rp_k2
        real(kind=real_kind),dimension(size(rp_xt,1)) :: rp_k3,rp_k4,rp_g,g1,rp_g1,g2
        real(kind=real_kind),dimension(size(rp_xt,1)) :: rp_g2,g3,rp_g3
        !
        !c
        integer :: nx, nt
        !
        nx = size(u_bg,1)
        nt = ubound(u_bg,2)
        !c
        !c begin computations
        !c
        !
        !c  if(icycle .eq. 1 .and. iter .eq. 1) then
        ! bg_fname = make_cycle_fname('DATA/ubg',icycle)
        call read_binary_data(u_bg, bg_fname)
        !c   endif
        !c
        x = u_bg(:,0)
        rp_xt(:,0) = ad_xt(:,0)
        rp_x = rp_xt(:,0)
        !
        do i = 1, NT
            !
            ! C**     Computing background terms
            !
            x = u_bg(:,i-1) 
            call nl_lorenz98_sub(x,tend)
            g1 = x + p5*dt*tend
            call nl_lorenz98_sub(g1,tend)
            g2 = x + p5*dt*tend
            call nl_lorenz98_sub(g2,tend)
            g3 = x + dt*tend
            !
            !c         xx(:,i) = x !! save for adjoint
            call rp_lorenz98_sub(x,rp_x,rp_g)
            rp_k1 = dt*rp_g
            rp_g1 = rp_x + p5*rp_k1             
            call rp_lorenz98_sub(g1,rp_g1,rp_g)         
            rp_k2 = dt*rp_g
            rp_g2 = rp_x + p5*rp_k2
            call rp_lorenz98_sub(g2,rp_g2,rp_g)
            rp_k3 = dt*rp_g
            rp_g3 = rp_x + rp_k3
            call rp_lorenz98_sub(g3,rp_g3,rp_g)
            rp_k4 = dt*rp_g
            !
            rp_x = rp_x + (rp_k1+2.0*rp_k2+2.0*rp_k3+rp_k4)/6.0+ad_xt(:,i)
            !
            rp_xt(:,i) = rp_x
            !
        enddo
        !c
    end subroutine rp_lorenz98nl
    !
    !c
    !c lorenz98_sub is next
    !c
    subroutine rp_lorenz98_sub(y,rp_y,rp_tend)
        !c
        !c declare passed variables
        !c
        !c      integer ndim,kdim
        real(kind=real_kind),dimension(:), intent(in) :: y,rp_y
        real(kind=real_kind),dimension(:), intent(out) :: rp_tend
        !c
        !c declare local variables
        !c
        real(kind=real_kind),dimension(size(y)) :: tend
        integer(kind=int_kind) :: i,j,l,l1,l2,l3,l4
        real(kind=real_kind) :: ts,rp_ts
        !c
        integer :: nx
        !
        nx = size(y)
        !c
        !c calculate terms
        !c
        tend=0
        rp_tend=0
        do l = 1 , NX
            ts    = 0.0
            rp_ts = 0.0
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
                    rp_ts = rp_ts + (rp_y(l3)*y(l4) + y(l3)*rp_y(l4) - &
                            & rp_y(l1)*y(l2) - y(l1)*rp_y(l2))/(NK**2.0)
                enddo
            enddo
            rp_tend(l) = rp_ts    - rp_y(l)
        enddo
        !c
    end subroutine rp_lorenz98_sub
    !
end module rp_lorenz98nl_mod
