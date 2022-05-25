!> @file covariance.f90
!!
!<
!> @brief covariance module
!!
!<
module covariance
    use constants, only : real_kind
    use modern_tools
    use cov_tools
implicit none
    !
    ! simple derrived type for the covariance
    type cov_type
        integer :: n = -1 ! size of the covariance matrix
        real(real_kind), dimension(:,:), allocatable :: corr ! correlation matrix
        real(real_kind)   :: thetat, dtbig ! model error parameters
        real(real_kind), dimension(:), allocatable :: var  ! variance?
        real(real_kind), dimension(:), allocatable :: std  ! standard deviation
        ! This is not used as the variance should be used.
        ! What is its relationship to the normalization factors?
        real(real_kind), dimension(:), allocatable :: Normfac  ! What is the actual purpose of this?
        ! 1) I did not understand the logic of the computations.
        ! 2) It is not used, it is computed and then set to 1.
        ! 
        ! In the end, we need to remove var and Normfac or make sure
        ! they are consistent with the logic.
    contains
        procedure :: write => write_cov
        procedure :: read => read_cov
        procedure :: init => init_cov
        procedure :: set_param
        procedure :: compute_normalization
        procedure :: allocate => allocate_cov
        procedure :: deallocate => deallocate_cov
        !
        procedure :: apply_traj
        procedure :: apply_state
        procedure :: apply_state_inplace
        !
        generic :: apply => apply_traj, apply_state, apply_state_inplace
    end type cov_type
    !
    private
    public :: cov_type
    !
contains
    !
    !> @brief Initialize the covariance matrix data structure
    !! @param[in, out] this covariance data structure
    !! @param[in] var
    !! @param[in] n size of the vector that the covariance appies to.
    !<
    subroutine init_cov( this, var, dlenght, n, thetat, dtbig )
        class(cov_type), intent(in out) :: this
        real(real_kind), intent(in)    :: var
        real(real_kind), intent(in) :: dlenght
        integer, intent(in) :: n
        real(real_kind), optional, intent(in)    :: thetat, dtbig
        !
        call this%allocate( n )
        call sqrt_exp_cov( this%corr, dlenght, this%Normfac )
        !
        this%var = var
        !
        this%std = sqrt( this%var )
        !
        call this%set_param( thetat, dtbig )
        !
    end subroutine init_cov
    !
    subroutine set_param( this, thetat, dtbig )
        class(cov_type), intent(in out) :: this
        real(real_kind), optional, intent(in)    :: thetat, dtbig
        !
        if(present(thetat))then
            this%thetat = thetat
        else
            this%thetat = 1.0
        end if
        !
        if(present(dtbig))then
            this%dtbig = dtbig
        else
            this%dtbig = 1.0
        end if
        !
    end subroutine set_param
    !
    !> @brief compute normalization factors
    !! @
    !! @
    !<
    subroutine compute_normalization(this)
        class(cov_type), intent(in out) :: this
        ! Local variables
        integer :: k
        !
        do k = 1, this%n
            this%Normfac(k) = 1.0/sum( this%corr(:,k) )
        end do
        !
    end subroutine compute_normalization
    !
    !> @brief Initialize the covariance matrix data structure
    !! @param[in, out] this covariance data structure
    !! @param[in] n size of the vector that the covariance appies to.
    !<
    subroutine allocate_cov( this, n )
        class(cov_type), intent(in out) :: this
        integer, intent(in) :: n
        !
        this%n = n
        !
        if(allocated(this%corr)) call this%deallocate()
        !
        allocate( this%corr(n,n), this%Normfac(n), this%var(n), this%std(n) )
        !
    end subroutine allocate_cov
    !
    !> @brief Initialize the covariance matrix data structure
    !! @param[in, out] this covariance data structure
    !! @param[in] n size of the vector that the covariance appies to.
    !<
    subroutine deallocate_cov( this )
        class(cov_type), intent(in out) :: this
        !
        if( allocated(this%corr) )then
            deallocate( this%corr, this%Normfac, this%var, this%std )
        end if
        !
    end subroutine deallocate_cov
    !
    !> @brief read the covariance matrix from file
    !! @param[in, out] this covarance data structure
    !! @param[in] var
    !! @param[in] n size of the vector that the covariance appies to.
    !! @param[in] fName path to the covariance file.
    !<
    subroutine read_cov( this, n, fName, thetat, dtbig )
        class(cov_type), intent(in out) :: this
        integer, intent(in) :: n
        character(*), intent(in) :: fName
        real(real_kind), optional, intent(in)    :: thetat, dtbig
        !
        integer :: k, i, j
        !
        call this%allocate( n )
        !
        write(*,"('reading cov <= ', A)") trim(fName)
        call read_binary_data( this%corr, fName )
        ! Extract the variance
        do k = 1, this%n
            this%var(k) = abs(this%corr(k,k))
        end do
        this%std = sqrt(this%var)
        ! Extract the correlation matrix, that is what we want for operations
        do j = 1, this%n
            do i = 1, this%n
                this%corr(i,j) = this%corr(i,j)/(this%std(i)*this%std(j))
            end do
        end do
        ! Compute the normalization factors
        call this%compute_normalization()
        !
        call this%set_param( thetat, dtbig )
        !
    end subroutine read_cov
    !
    !> @brief write the covariance matrix to a file
    !! @param[in] this covarance data structure
    !! @param[in] fName path to the covariance file
    !<
    subroutine write_cov( this, fName )
        class(cov_type), intent(in) :: this
        character(*), intent(in) :: fName
        ! local variables
        integer :: i, j
        real(real_kind), dimension(:,:), allocatable :: cov
        !
        allocate( cov(this%n, this%n) )
        ! Apply the variance to save the covariance
        write(*,"('writing cov => ', A)") trim(fName)
        do j = 1, this%n
            do i = 1, this%n
                cov(i,j) = this%corr(i,j)*(this%std(i)*this%std(j))
            end do
        end do
        !
        call write_binary_data( cov, fName )
        deallocate ( cov )
        !
    end subroutine write_cov
    !
    !> @brief apply the covariance to a state vector
    !! @param[in] this covarance data structure
    subroutine apply_state(this, x, y)
        class(cov_type), intent(in) :: this
        real(real_kind), dimension(:), intent(in) :: x
        real(real_kind), dimension(:), intent(out) :: y
        ! local variables
        real(real_kind), dimension(size(x,1)) :: temp
        !!
        !temp = x*this%Normfac
        !y    = matmul( this%corr, temp )
        !y    = this%var*this%Normfac*y
        !!
        !
        temp = this%std*x
        y    = matmul( this%corr, temp )
        y    = this%std*this%Normfac*y
    end subroutine apply_state
    !
    !> @brief apply the covariance to a state vector
    !! @param[in] this covarance data structure
    subroutine apply_state_inplace(this, x)
        class(cov_type), intent(in) :: this
        real(real_kind), dimension(:), intent(in out) :: x
        ! local variables
        real(real_kind), dimension(size(x,1)) :: temp
        !!
        !temp = x*this%Normfac
        !x    = matmul( this%corr, temp )
        !x    = this%var*this%Normfac*x
        !!
        temp = this%std*x
        x    = matmul( this%corr, temp )
        x    = this%std*this%Normfac*x
        !
    end subroutine apply_state_inplace
    !
    !> @brief apply the covariance to a trajectory, used fo model error
    !! @param[in] this covarance data structure
    !! @param[in, out] x trajectory
    !<
    subroutine apply_traj(this, x)
        class(cov_type), intent(in) :: this
        real(real_kind), dimension(:,0:), intent(in out) :: x
        ! local variables
        !real(real_kind),dimension(size(x,1)) :: temp
        real(real_kind),allocatable,dimension(:,:) :: r1
        real(real_kind) :: thetat, dtbig
        integer :: nt, i, k
        ! shorten some variable names
        thetat = this%thetat
        dtbig  = this%dtbig
        !
        nt = ubound(x, 2)
        !
        allocate( r1(size(x,1),0:nt+1) )
        r1 = 0.0
        !
        ! Innocent change the beginning index from 2 to 1
        do i = 1 , NT+1
            r1(:,i) = r1(:,i-1) &
                &  - dtbig*( thetat*r1(:,i-1)+2.0*thetat*x(:,i-1) )
        enddo
        x(:,NT) = -r1(:,NT+1)/(2.0*thetat)&
            &       -dtbig*( thetat*(-r1(:,NT+1)/(2.0*thetat))+r1(:,NT+1) )
        ! Innocent change the ending index from 1 to 0
        do i = NT-1, 0, -1
            x(:,i) = x(:,i+1) &
                & -dtbig*( thetat*x(:,i+1)+r1(:,i+1) )
        enddo
        !
        r1 = 0
        !matrix
        do k = 1 , NT
            !x(:,k) = x(:,k)*this%Normfac(:)
            !temp = matmul( this%corr, x(:,k) )
            !r1(:,k) = temp*this%Normfac
            call this%apply(x(:,k))
        enddo
        !
        !x(:,1:NT) = this%var*r1( :,1:NT )
        !
        deallocate(r1)
        !
    end subroutine apply_traj
end module covariance