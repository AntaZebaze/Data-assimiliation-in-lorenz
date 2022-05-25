!> @file cov_tools.f90
!!
!<

!> @brief covariance tools
!!
!<
module cov_tools
    use constants, only : real_kind
implicit none
    !
    private
    public :: sqrt_exp_cov
    !
contains
    !
    subroutine sqrt_exp_cov(cov, dlenght, Normfac)
        real(real_kind), dimension(:,:), intent(out) :: cov
        real(real_kind), dimension(:), intent(out) :: Normfac
        real(real_kind), intent(in) :: dlenght
        ! Local variables
        !real (real_kind),dimension(size(Normfac)) :: delt
        integer :: i, j
        !
        !! square exponential covariance for the model error
        cov = 0
        do i = 1, size(cov,1)
            do j = 1 , size(cov,2)
                if( abs(i-j)<min(60, size(cov,1) ) )then  
                    cov(i,j) = 1.0*exp( -((i-j)/dlenght)**2.0 )
                end if
            end do
        end do
        !
        ! c delt = 0.0
        ! c delt(size(Normfac)/2) = 1.0
        do i = 1 , size(Normfac)
            ! delt = 0.0
            ! delt(i) = 1.0
            ! Normfac(i) = 0.0
            ! do j = 1 , size(Normfac)
            !     !Normfac(i) = Normfac(i) + cov(i,j)*delt(j)
            !     Normfac(i) = sum( cov(:,i) )
            ! end do
            ! ! c  write(*,*) i,delt(i),Normfac(i)
            ! ! 
            ! ! So the calculations are useless
            ! ! 
            ! Normfac(i) = 1.0 !/sqrt( Normfac(i) )
            !
            Normfac(i) = 1.0/sum( cov(:,i) )
        end do
        
        !
    end subroutine sqrt_exp_cov
end module cov_tools