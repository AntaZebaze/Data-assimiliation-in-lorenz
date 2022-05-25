!> @file observations.f90
!!
!<
!> @brief observations module
!!
!<
module observations
    !
    use constants, only : real_kind, int_kind
    !
implicit none
    !
    ! Simple type for observations, to make life easier
    type obs_op
        integer (kind=int_kind),dimension(:), allocatable :: dloc,tloc
        integer :: nobs
        integer :: nothing
    contains
        procedure :: init => init_obs_op
        procedure :: allocate => allocate_obs_op
        procedure :: deallocate => deallocate_obs_op
        procedure :: get_nobs
        procedure :: apply
        procedure :: apply_t ! transpose
        !procedure :: write => write_obs_op
        !procedure :: read => read_obs_op
    end type obs_op
    !
    private
    public :: obs_op
    !
contains
    !
    subroutine do_nothing(nothing)
        integer, intent(in) :: nothing
        !
        integer, save :: internal_nothing
        !
        internal_nothing = nothing
        !
    end subroutine do_nothing
    !
    function get_nobs(this, nx, nt, xspacing, tspacing, x_starting, t_starting)result(nobs)
        class(obs_op), intent(in) :: this
        integer, intent(in) :: nx, nt, xspacing, tspacing
        integer, optional, intent(out) :: x_starting, t_starting
        ! local variables
        integer :: x_nobs, t_nobs, nobs, x_mod, t_mod, il_tstart, il_xstart, nx2
        !
        call do_nothing( this%nothing ) ! just to prevent warning at compiliation
        !
        ! Because of periodicity in x, we do not want to sample the first and the last
        ! we use a modified nx2
        nx2 = nx -1
        !
        x_mod = modulo(nx2, xspacing)
        t_mod = modulo(nt, tspacing)
        x_nobs = nx2/xspacing
        t_nobs = nt/tspacing
        !
        if(t_mod==0)then
            il_tstart = tspacing
        else
            il_tstart = t_mod
            t_nobs = t_nobs+1
        end if
        ! Because of periodicity in x, we do not want to sample the first and the last
        if(x_mod==0)then
            il_xstart = max(xspacing/2, 1)
        else
            il_xstart = x_mod
            x_nobs = x_nobs+1
        end if
        !
        nobs = t_nobs*x_nobs
        !
        if(present(x_starting)) x_starting = il_xstart
        if(present(t_starting)) t_starting = il_tstart
        !
    end function get_nobs
    !
    !> @brief Initialize the observations operator data structure
    !! @param[in, out] this observations operator data structure
    !! @param[in] var
    !! @param[in] nx, nt size of the vector that the observations operator appies to.
    !! @param[in] xspacing number of time steps between two observations
    !! @param[in] tspacing number of space steps between two observations
    !! @note, when calling the apply subroutine, we expect the state state
    !! to be one time step larger than nt. We do not want to sample t0,
    !! because in this twin experiment, we can end up using them twice
    !<
    subroutine init_obs_op( this, nx, nt, xspacing, tspacing )
        class(obs_op), intent(in out) :: this
        integer, intent(in) :: nx, nt, xspacing, tspacing
        !
        integer :: idata, total_nobs, t_idx, x_idx, x_start, t_start
        !
        total_nobs = this%get_nobs( nx, nt, xspacing, tspacing, x_start, t_start )
        call this%allocate( total_nobs )
        !
        !print*, "nx, nt, xspacing, tspacing, total_nobs = ", nt, nx, xspacing, tspacing, total_nobs
        !
        idata = 0
        do t_idx = t_start, nt, tspacing
            do x_idx = x_start, nx-1, xspacing ! because of periodicity we do not reach nx
                idata = idata+1
                !print*, "t_idx, x_idx, idata= ", t_idx, x_idx, idata
                this%dloc( idata ) = x_idx
                this%tloc( idata ) = t_idx
            end do
        end do
        !
        !print *, this%dloc
        !print *, this%tloc
        !pause
        !
    end subroutine init_obs_op
    !
    !> @brief Initialize the observations operator data structure
    !! @param[in, out] this observations operator data structure
    !! @param[in] n size of the vector that the observations vector.
    !<
    subroutine allocate_obs_op( this, n )
        class(obs_op), intent(in out) :: this
        integer, intent(in) :: n
        !
        if(allocated(this%dloc)) call this%deallocate()
        !
        allocate( this%dloc(n), this%tloc(n) )
        !
    end subroutine allocate_obs_op
    !
    !> @brief Initialize the observations operator data structure
    !! @param[in, out] this observations operator data structure
    !! @param[in] n size of the vector that the observations operator appies to.
    !<
    subroutine deallocate_obs_op( this )
        class(obs_op), intent(in out) :: this
        !
        if( allocated(this%dloc) )then
            deallocate( this%dloc, this%tloc )
        end if
        !
    end subroutine deallocate_obs_op
    !
    !> @brief
    !!
    !<
    subroutine apply( this, u_in, y_out )
        !
        class(obs_op), intent(in) :: this
        real (real_kind),dimension(:,0:),intent(in) :: u_in
        real (real_kind),dimension(:),intent(out)   :: y_out
        ! Local variables
        integer (int_kind) :: k, i1, j1
        !
        y_out = 0
        do k = 1, size(y_out)
            i1 = this%dloc(k)
            j1 = this%tloc(k)
            y_out(k) = u_in( i1, j1 )
        enddo
        !
    end subroutine apply
    !
    subroutine apply_t( this, y, implu) !, dtbig )
        !
        class(obs_op), intent(in) :: this
        real(real_kind), dimension(:)   , intent(in) :: y
        real(real_kind), dimension(:,0:), intent(in out) :: implu
        !real(real_kind), intent(in) :: dtbig
        !local variables
        integer(int_kind)::idata,i1,t1
        !
        do idata = 1, size(y)
            i1 = this%dloc(idata)
            t1 = this%tloc(idata)
            implu(i1,t1) = implu(i1,t1) + y(idata)
        enddo
        !c          implu  = implu/(M*dt)
        !implu  = implu/(dtbig) ! this must be taken out of here
        !
    end subroutine apply_t
    !
end module observations