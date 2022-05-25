module convol_mod
    use constants
    use modern_tools
    use cov_tools
implicit none
    !
    real(real_kind), dimension(:,:), allocatable :: covI, covM
    real(real_kind), dimension(:), allocatable   :: NormfacI, NormfacM
    logical :: initialized = .false.
    private
    public convol, make_cov, write_model_cov, write_ic_cov
    !
 contains
    !
    !
    subroutine convol( sigmau )

        !c/////////////////////////////////////////////////////////////////////////////
        !c//
        !c//     Time and spatial convolution 
        !c//
        !c/////////////////////////////////////////////////////////////////////////////
        !
        real (kind=real_kind),dimension(:,0:), intent(in out) :: sigmau
        ! Local variables
        real (kind=real_kind),dimension(NX) :: temp
        integer (kind=int_kind)::i,j,k
        ! real (kind=real_kind)::theta1

        real(kind=real_kind),allocatable,dimension(:,:)::r1
        !
        
        if( .not.initialized )then
            write(*,*), "In convol, covariance not initialize, call make_cov"
            stop
        end if
        !
        !
        if(weak_constraint)then
            allocate(r1(NX,NT+1))
            r1(:,:) = 0.0

            ! C**    time convolution of interior residuals

            ! r1(:,1) = 0.0  ! useless
            do i = 2 , NT+1
                r1(:,i) = r1(:,i-1) &
                    &  - dtbig*(thetat*r1(:,i-1)+2.0*thetat*sigmau(:,i-1))
            enddo
            sigmau(:,NT) = -r1(:,NT+1)/(2.0*thetat)&
                &       -dtbig*(thetat*(-r1(:,NT+1)/(2.0*thetat))+r1(:,NT+1))
            do i = NT-1, 1, -1
                sigmau(:,i) = sigmau(:,i+1) &
                    & -dtbig*( thetat*sigmau(:,i+1)+r1(:,i+1) )
            enddo
            !
            deallocate(r1)
            !
            allocate( r1(NX,0:NT) )
            r1 = 0
            !
            do k = 1 , NT
                sigmau(:,k) = sigmau(:,k)*NormfacM(:)
                temp = matmul( covM, sigmau(:,k) )
                !temp = 0
                !do i = 1 , NX
                !    temp(i) = 0.0
                !    do j = 1 , NX
                !        temp(i) = temp(i) + covM(i,j)*sigmau(j,k)
                !    enddo
                !enddo
                r1(:,k) = temp*NormfacM
            enddo
            !
            sigmau(:,1:NT) = varru*r1( :,1:NT )
            !
            deallocate(r1)
        end if
        !
        k = 0
        sigmau(:,k) = sigmau(:,k)*NormfacI(:)
        temp        = matmul( covI, sigmau(:,k) )
        sigmau(:,k) = varruI*NormfacI*temp
        !
    end subroutine convol
    !
    !> @brief compute the covariance and the normalization factors
    !! @param[in] dlenght decorrelation length scale
    !!
    !<
    subroutine make_cov( dlenght, cycle_num )
        real(real_kind), intent(in) :: dlenght
        integer, intent(in) :: cycle_num
        ! Local variables
        integer :: k!i, j
        character(256) :: fName
        !
        if( initialized )then
            deallocate( covI, covM, NormfacI, NormfacM )
        end if
        !
        allocate( covI(NX,NX), NormfacI(NX), covM(NX,NX),NormfacM(NX) )
        !
        ! Allocate 
        if(weak_constraint)then
            write(*,*) "============================================"
            write(*,*) "Running weak constraints"
            write(*,*) "============================================"
            !
            call sqrt_exp_cov(covM, dlenght, NormfacM)
        else
            write(*,*) "============================================"
            write(*,*) "Running Strong constraints"
            write(*,*) "============================================"
            covM = 0
            do k = 1, NX
                covM(k,k) = 1.0
            end do
            NormfacM = 1.0
        end if
        ! Use static covariance for cycle 1
        ! and when it is explicitly prescribed
        if(use_static_cov.or.(cycle_num==1))then
            write(*,*) "============================================"
            write(*,*) "Using the square exponential for IC error covariance"
            write(*,*) "============================================"
            call sqrt_exp_cov(covI, dlenght, NormfacI)
        else
            fName = make_cycle_fname( pert_dir, "cov", cycle_num )
            if( fileExist(fName) )then
                write(*,*) "============================================"
                write(*,*) "Using dynamic covariance from file <"//trim(fName)//">"
                write(*,*) "============================================"
                call read_binary_data( covI, fName )
                NormfacI = 1.0 ! /sqrt( NormfacI( k ) )
                varruI = 1.0
                write(*,"('min/max covI', 2f9.6)"), minval(covI), maxval(covI)
            else
                write(*,*) "Covariance file not found <"//trim(fName)//">"
                write(*,*) "Consider setting use_static_cov to True"
                write(*,*) "to use the square exponential covariance"
                stop
            end if
        end if
        !
        write(*,"('min/max covM = ',2f13.3)") minval(covM),maxval(covM)
        write(*,"('min/max covI = ',2f13.3)") minval(covI),maxval(covI)
        write(*,"('min/max NormfacM = ',2f13.3)") minval(NormfacM),maxval(NormfacM)
        write(*,"('min/max NormfacI = ',2f13.3)") minval(NormfacI),maxval(NormfacI)
        !
        initialized = .true.
        !
    end subroutine make_cov
    !
    subroutine write_model_cov( fName )
        !
        character(*), intent(in) :: fName
        !
        call write_binary_data( covM, fName )
        !
    end subroutine write_model_cov
    !
    subroutine write_ic_cov( fName )
        !
        character(*), intent(in) :: fName
        !
        call write_binary_data( covI, fName )
        !
    end subroutine write_ic_cov
    !
end module convol_mod
