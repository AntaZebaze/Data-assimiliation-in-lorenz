!> @file model_params.f90
!!
!<

!> @brief model parameters module
!!
!<
module model_params
    !
    use constants, only : real_kind
    !
implicit none
    !
    type model_param
        logical :: is_allocated = .false. ! 
        integer :: nx = -1 ! number of grid points in x
        integer :: nt = -1 ! number of time steps
        real(real_kind) :: dx
        real(real_kind) :: dt
        real(real_kind) :: dtbig
        character(256) :: bg_fname ! background file name for adjoint and TLM
        character(256) :: beta_fpref ! beta file prefix for the CG solver
        !
        ! Extra forcing for the adjoint, used in time parallel
        real(real_kind), dimension(:), allocatable :: uT_ad
        real(real_kind), dimension(:), allocatable :: u0_ad
        !
    contains
        procedure :: init
        procedure :: allocate => allocate_mp
        procedure :: deallocate => deallocate_mp
    end type model_param
    !
    private
    public :: model_param
    !
contains
    !
    subroutine init( this, nx, nt, dt, dtbig, bg_fname, beta_fpref )
        class(model_param), intent(in out) :: this
        integer, intent(in) :: nx, nt
        real(real_kind), intent(in) :: dt, dtbig
        character(*), intent(in) :: bg_fname, beta_fpref
        !
        this%nx = nx
        this%nt = nt
        this%dt = dt
        this%dtbig = dtbig
        this%bg_fname = trim( bg_fname )
        this%beta_fpref = trim( beta_fpref )
        call this%allocate(nx)
        !
    end subroutine init
    !
    subroutine allocate_mp( this, nx )
        class(model_param), intent(in out) :: this
        integer, intent(in) :: nx
        !
        if(this%is_allocated)then
            call this%deallocate()
        end if
        allocate( this%uT_ad(nx) )
        allocate( this%u0_ad(nx) )
        this%is_allocated = .true.
        this%uT_ad = 0
        this%u0_ad = 0
        !
    end subroutine allocate_mp
    !
    subroutine deallocate_mp( this )
        class(model_param), intent(in out) :: this
        !
        if(this%is_allocated)then
            deallocate( this%uT_ad )
            deallocate( this%u0_ad )
            this%is_allocated = .false.
        end if
        !
    end subroutine deallocate_mp
    !
end module model_params
