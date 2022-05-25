!
!
module psas_solver
    use constants, only: real_kind
    use modern_tools
    use model_tools
    use covariance
    use observations
    use model_params
    !
implicit none
!
contains
    !
    !> @brief run in a subwindow and return the increment
    !! ( all subwindows are run in parallel )
    !!
    !<
    subroutine compute_4dvar_inc( mp, Niter_max, obs_std2, tol&
                                 &, u_bg, u_inc, obs&
                                 &, H, ic_cov, me_cov )
        !
        type(model_param), intent(in out) :: mp
        integer, intent(in) :: Niter_max
        real(real_kind) :: obs_std2, tol
        real(real_kind), dimension(:, 0:), intent(in)  :: u_bg 
        real(real_kind), dimension(:, 0:), intent(out) :: u_inc
        real(real_kind), dimension(:), intent(in)      :: obs
        type(obs_op), intent(in) :: H
        type(cov_type), intent(in) :: ic_cov
        type(cov_type), optional, intent(in) :: me_cov
        !
        ! Local variables
        real (kind=real_kind),dimension(size(obs)) :: Hx, beta
        !
        write(*,"(A)"), '---------------------------------------------------'
        write(*,"(A)"), "Entering compute_4dvar_inc"
        !
        ! innovation
        call H%apply( u_bg, Hx )
        Hx = obs - Hx
        write(*,"('min, max Hx ', 2ES10.2E2)"), minval(Hx), maxval(Hx)
        ! Solve the DA problem
        call psas_solve( Hx, beta, mp, Niter_max, obs_std2, tol&
                        &, H, ic_cov, me_cov )
        ! Compute the increment
        call MBGT( mp, beta, u_inc, H, ic_cov, me_cov )
        write(*,"(A)"), "Exiting compute_4dvar_inc"
        write(*,"(A)"), '---------------------------------------------------'
        !
    end subroutine compute_4dvar_inc
    !
    !> @brief Solve the DA problem using RBCG
    !! @param [in] innov innovation of one member of the ensemble
    !! @param [in,out] beta1 solution of the system
    !! @param [in] dloc
    !! @param [in] tloc
    !! @param [in] Niter_max maximum number of iterations
    !! @param [in] tol stopping criteria based on the B-norm of the gradient
    !!
    !<
    subroutine psas_solve(innov, beta1, mp, Niter_max, obs_std2, tol&
                          &, H, ic_cov, me_cov)
        real(real_kind), dimension(:), intent(in) :: innov
        real(real_kind), dimension(:), intent(in out) ::beta1
        type(model_param), intent(in out) :: mp
        integer, intent(in) :: Niter_max
        real(real_kind) :: obs_std2, tol
        type(obs_op), intent(in) :: H
        type(cov_type), intent(in) :: ic_cov
        type(cov_type), optional, intent(in) :: me_cov
        ! local variables
        ! RBCG variables
        real(real_kind), dimension(size(innov)) :: cg_r
        ! real(real_kind), dimension(size(innov)) :: cg_rtmp
        real(real_kind), dimension(size(innov)) :: cg_lambda
        real(real_kind), dimension(size(innov)) :: cg_t
        real(real_kind)                :: cg_beta
        real(real_kind), dimension(size(innov)) :: cg_p
        real(real_kind)                :: cg_alpha
        real(real_kind), dimension(size(innov)) :: cg_w
        real(real_kind), dimension(size(innov)) :: cg_q
        ! others
        real(real_kind) :: r0_norm ! B-norm of the residue at the zeroth iteration
        real(real_kind) :: ri_norm ! B-norm of the residue at the current iteration
        real(real_kind) :: cg_wTr !dot product w by r (i)
        real(real_kind) :: cg_qTt !dot product q by t
        real(real_kind) :: cg_wTrp1 !dot product w by r (i+1)
        ! real(real_kind) :: cg_gNorm !b-norm of the gradient
        real(real_kind) :: cg_omega, cg_epsilon
        integer :: cg_itr
        ! integer :: j
        logical :: converged
        !
        write(*,"('RBGC - nItrMax = ',i4.1, ', tol =',ES10.2E2)")Niter_max, tol
        !write(*,*)"min, max, norm innov", minval(innov), maxval(innov), norm2(innov)
        cg_epsilon = epsilon(1.0_real_kind)
        cg_itr = 0
        cg_lambda = 0.0
        cg_r = innov/obs_std2
        cg_p = cg_r
        ! computing HA(HA)^Tr
        ! call matrixmult(nRow, nCol, HA, cg_r, cg_w)
        !
        !
        !write(*,*)"RBCG - initialization, calling GBGT"
        call GBGT(mp, cg_r, cg_w, H, ic_cov, me_cov )
        !write(*,*)"RBCG - initialization, after GBGT"
        cg_t = cg_w
        !
        cg_wTr = dot_product(cg_w, cg_r)
        ri_norm = sqrt(cg_wTr)
        r0_norm    = ri_norm
        cg_omega= ri_norm/r0_norm
        cg_itr = 1
        write(*,"('RBGC - iter',i4.1,' of',i4.1, ', r0_norm, wTr =',2ES10.2E2)")&
        & cg_itr, Niter_max, r0_norm, cg_wTr
        !write(*,*)"min, max, norm cg_w", minval(cg_w), maxval(cg_w), norm2(cg_w)
        converged = .false.
        !
        do while( (cg_itr<= Niter_max).and.(.not.converged) )
            !
            !write(*,*)"RBCG - cg_itr=", cg_itr
            cg_q = cg_t/obs_std2 + cg_p
            ! The following two quantities are used as denominator
            ! So if one is zero, either, there is convergence or
            ! something is wrong
            cg_wTr = dot_product(cg_w, cg_r)
            cg_qTt = dot_product(cg_q, cg_t)
            if( (cg_qTt>cg_epsilon).and.(cg_wTr>cg_epsilon) )then
                cg_alpha = cg_wTr/cg_qTt
                cg_lambda = cg_lambda + cg_alpha*cg_p
                cg_r = cg_r - cg_alpha*cg_q
                !re-conjugate r
                !call matrixmult( nRow, nCol, HA, cg_r, cg_w )
                ! -------------------------------------------
                ! Writing beta to file
                call write_beta_nodir(cg_lambda, mp%beta_fpref, cg_itr)

                call GBGT(mp, cg_r, cg_w, H, ic_cov, me_cov )
                cg_wTrp1 = dot_product(cg_w, cg_r)
                cg_beta = cg_wTrp1/cg_wTr
                if( cg_beta<cg_epsilon )then
                    converged = .true.
                end if
                cg_p = cg_r + cg_beta*cg_p
                cg_t = cg_w + cg_beta*cg_t
                ri_norm = sqrt(cg_wTrp1)
                cg_omega= ri_norm/r0_norm

                !Check the convergence
                if (cg_omega > tol) then
                    write(*,"('RBGC - iter',i4.1,' of',i4.1, ', cg_omega =',ES10.2E2)")&
                    & cg_itr, Niter_max, cg_omega
                else
                    converged = .true.
                endif
                !
                !Moving to the next iteration
                !
                cg_itr=cg_itr+1
            else
                ! the denominator of alpha is zero or the denominator
                ! of beta is zero
                converged = .true.
            end if
        end do
        !
        if( .not.converged ) then
            write(*,"(A)") "RBGC - Maximum number of iterations reached", Niter_max
        else
            write(*,"('RBCG converged, total nb iterations ',i4.1)"), cg_itr
            write(*,"('Solution ||cg_lambda|| ',ES10.2E2)"), sum(abs(cg_lambda))
            if(cg_qTt<=cg_epsilon)then
                write(*,"('q*t is too small ',ES10.2E2)"), cg_qTt
            end if
            if(cg_wTr<=cg_epsilon)then
                write(*,"('w*r is too small ',ES10.2E2)"), cg_wTr
            end if
            if(cg_beta<=cg_epsilon)then
                write(*,"('beta is too small ',ES10.2E2)"), cg_beta
            end if
        end if
        !
        beta1 = cg_lambda
    end subroutine psas_solve
    !
end module psas_solver
!