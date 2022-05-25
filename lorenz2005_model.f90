!> @file model.f90
!!
!<

module lorenz2005_model
    use constants, only : real_kind
implicit none
    private
    public :: lorenz2005_traj
    !
contains
    !
    !>@brief RHS of the Lorenz 2005 model
    !!
    !<
    subroutine lorenz2005_rhs(u, du, nk, F)
        real(real_kind),dimension(:), intent(in) :: u
        real(real_kind),dimension(:), intent(out) :: du
        integer, intent(in) :: nk
        real(real_kind), intent(in) :: F
        ! Local variables
        integer :: nx, i, j, l, l1, l2, l3, l4
        real(real_kind):: ts
        !
        nx = size(u)
        du = 0
        do l = 1 , NX
            ts    = 0.0
            do j = -nk/2, nk/2
            do i = -nk/2,nk/2
                l1=l-2*nk-i
                if (l1 .lt. 1) l1 = NX +l1-1
                l2=l-nk-j
                if (l2 .lt. 1) l2 = NX+l2-1
                l3 = l-nk+j-i
                if (l3 .lt. 1) l3 = NX +l3-1
                l4 = l+nk+j
                if(l4 .gt. NX) l4 = l4-NX
                ts = ts + (u(l3)*u(l4) - u(l1)*u(l2))/(nk**2.0)
            enddo
            enddo
            du(l) = ts - u(l) + F
        enddo
        !
    end subroutine lorenz2005_rhs
    !
    !>@brief Integration of the Lorenz 2005 model
    !!
    !<
    subroutine lorenz2005_traj(x0, trj, dt, nk, F)
        real(real_kind),dimension(:), intent(in) :: x0
        real(real_kind),dimension(:,0:), intent(out) :: trj
        integer, intent(in) :: nk
        real(real_kind), intent(in) :: dt, F
        !
        integer :: nx, nt, i
        real(real_kind),dimension(size(x0)) :: u, k1,k2,k3,k4,g,g1,g2,g3
        real(real_kind), parameter :: p5 = 0.5
        !
        nx = size(x0)
        nt = ubound(trj, 2)
        if( size(trj,1)/= nx )then
            write(*,"(A)") "In lorenz2005_traj, inconsistent array sizes"
            write(*,"(A)") "  Expecting the first dimension of trj to have"
            write(*,"(A)") "  the same size as x0. Instead, got:"
            write(*,"(A,I5.1)") "  size(trj,1) = ", size(trj,1)
            write(*,"(A,I5.1)") "  size(x0) = ", size(x0)
            stop
        end if
        !
        u = x0
        trj(:,0) = u
        do i = 1, nt
            ! c         xx(:,i) = u !! save for adjoint
            call lorenz2005_rhs( u, g, nk, F )
            k1 = dt*g
            g1 = u + p5*k1
            call lorenz2005_rhs( g1, g, nk, F )
            k2 = dt*g
            g2 = u + p5*k2
            call lorenz2005_rhs( g2, g, nk, F )
            k3 = dt*g
            g3 = u + k3
            call lorenz2005_rhs( g3, g, nk, F )
            k4 = dt*g
            !
            u = u + ( k1+2.0*k2+2.0*k3+k4 )/6.0
            trj(:,i) = u
            !
            if(modulo(i,50)==0)then
                write(*,"('lorenz2005_traj ',I5.1,5F13.6)"), i, u(10),u(15),u(20),u(25),u(30)
            end if
        enddo
        !
    end subroutine lorenz2005_traj
    !
end module lorenz2005_model