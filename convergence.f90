
program convergence

    !c/////////////////////////////////////////////////////////////////////////
    !c//
    !c//     Main Driver
    !c//
    !c/////////////////////////////////////////////////////////////////////////

    use constants

implicit none

    integer (kind=int_kind)::i,j,k,i1,j1,iter,icycle,istep,i2

    real (kind=real_kind),dimension(NX,0:LT)::utr,xdata
    real (kind=real_kind),dimension(NX,0:NT)::u,rp_xt,diff,up

    real(kind=real_kind):: rms(100),rms_diff(100)
    logical res

    character*1 word(0:9)
    character*2 cycl(0:99)
    character*10  file_b(0:999)
    character*11  file_a(0:999)

    data word /'0','1','2','3','4','5','6','7','8','9'/

    istep = 0
    do j = 0 , 9
        do i = 0 , 9
            if(istep .le. NCYCLE) then
                cycl(istep) = word(j)//word(i)
                istep = istep + 1
                !c			write(*,*) 'cycle loop ',cycl(istep-1)
            endif
        enddo
    enddo

    !c		   STOP
        
    !C**       read in the data file for the measurement locations and times !!!

    open(1,file='TRUTH/xt_true.dat',status='unknown',form='unformatted')
        read(1) xdata
    close(1)

    do icycle = 4 , 4 ! M

        do k = 0, 9
            do j = 0, 9
                do i = 0, 9
                    file_a(100*k+10*j+i)='final01.'//word(k)//word(j)//word(i)
                enddo
            enddo
        enddo
        open(2,file='DATA/rms_04.dat',status='unknown',form='formatted')

        rms_diff(1) = 2.0
    
        iter = 100		  
        open(1,file=file_a(iter),status='unknown',form='formatted')
        do i = (icycle-1)*NT,icycle*NT
            i2 = i-(icycle-1)*NT
            do j = 1 , NX		 
                read(1,*)i1,j1, up(j,i2)
            enddo
        enddo
        do iter = 1 , 99

            !c		   inquire(file=file_a(iter),exist=res)

            !c		   if (res) then
            
            open(1,file=file_a(iter),status='unknown',form='formatted')
            do i = (icycle-1)*NT,icycle*NT
                    i2 = i-(icycle-1)*NT
                do j = 1 , NX		 
                    read(1,*)i1,j1, u(j,i2)
                enddo
            enddo
            diff(:,:)=(up(:,:)-u(:,:))**2.0

            !c		 if (iter .eq. 1) then
            !c		 up = u
            !c		 else
            !c		 diff = (u -up)**2.0
            !c		 up = u
            !c		 endif
            close(1)
            !c	   endif
            rms(iter) = sqrt(sum(diff)/NX/(NT+1))
            !c		if (iter .gt. 1) rms_diff(iter) = abs(rms(iter)-rms(iter-1))

            write(*,*) 'RMS per iter ==', iter, rms(iter)
            write(2,*) iter, rms(iter)! ,rms_diff(i)

        enddo ! (iter)

        close(2)
    
    enddo !  icycle

end program convergence
