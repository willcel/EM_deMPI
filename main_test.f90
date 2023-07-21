
    !***************************************************
    ! 地面电磁法的DE反演算法
    !-------------------------------------------------------------------------------------
    !                obj : The user provided file for evlauting the objective function.
    !                      subroutine obj(xc,fitness)
    !                      where "xc" is the real decision parameter vector.(input)
    !                            "fitness" is the fitness value.(output)
    !             Dim_XC : Dimension of the real decision parameters.
    !      XCmin(Dim_XC) : The lower bound of the real decision parameters.
    !      XCmax(Dim_XC) : The upper bound of the real decision parameters.
    !                VTR : The expected fitness value to reach.
    !                 NP : Population size.
    !            itermax : The maximum number of iteration.
    !               F_XC : Mutation scaling factor for real decision parameters.
    !              CR_XC : Crossover factor for real decision parameters.
    !           strategy : The strategy of the mutation operations is used in HDE.
    !            refresh : The intermediate output will be produced after "refresh"
    !                      iterations. No intermediate output will be produced if
    !                      "refresh < 1".
    !             iwrite : The unit specfier for writing to an external data file.
    ! bestmen_XC(Dim_XC) : The best real decision parameters.
    !              bestval : The best objective function.
    !             nfeval : The number of function call.
    !         method(1) = 0, Fixed mutation scaling factors (F_XC)
    !                   = 1, Random mutation scaling factors F_XC=[0, 1]
    !                   = 2, Random mutation scaling factors F_XC=[-1, 1]
    !         method(2) = 1, Random combined factor (F_CR) used for strategy = 6
    !                        in the mutation operation
    !                   = other, fixed combined factor provided by the user
    !         method(3) = 1, Saving results in a data file.
    !                   = other, displaying results only.
    !***************************************************
    Module ran_mod
    Implicit None
    ! ran return a uniform random number between 0-1
    ! normal return a normal distribution
    contains
    function ran()   !returns random number between 0 - 1
    implicit none
    integer , save :: flag = 0
    double precision :: ran
    if(flag==0) then
        call random_seed()
        flag = 1
    endif
    call random_number(ran)     ! built in fortran 90 random number function
    end function ran

    function normal(mean,sigma)
    implicit none
    integer :: flag
    double precision, parameter :: pi = 3.141592653589793239
    double precision :: u1, u2, y1, y2, normal, mean, sigma
    save flag
    data flag /0/
    u1 = ran(); u2 = ran()
    if (flag.eq.0) then
        y1 = sqrt(-2.0d0*log(u1))*cos(2.0d0*pi*u2)
        normal = mean + sigma*y1
        flag = 1
    else
        y2 = sqrt(-2.0d0*log(u1))*sin(2.0d0*pi*u2)
        normal = mean + sigma*y2
        flag = 0
    endif
    end function normal
    !The above codes are made in Fortran 90 language, if you have any question, you may write to sealin2008@hotmail.com
    End Module ran_mod

    module data_type
    implicit none
    integer(kind=4), parameter :: IB=4, RPD=8
    end module data_type

    program main
    use ran_mod
    use data_type

    implicit none

    integer(kind=IB), parameter :: NP=60, itermax=2, strategy=6, &
        refresh=10, iwrite=15
    integer(kind=IB), dimension(3), parameter :: method=(/0, 1, 1/)
    real(kind=8), parameter :: VTR=1.0e-30_RPD, CR_XC=0.5_RPD

    integer:: Dim_XC

    integer:: ntc,nolayer,npara,ns    ! number of time channel / layers / parameter per point / points
    integer:: i,j,k  ! index
    integer:: tot_para,tot_ntc
    integer:: size0,size1,size2
    integer:: l,ka

    real *8 sigma

    real *8,allocatable:: rho_true(:), hh_true(:)
    real *8,allocatable:: rho_iter(:), hh_iter(:)

    real *8 pi,snr


    real *8 sigma1,sigma2,sigma3,sigma4  ! 大地 混凝土 空气 水

    real *8,allocatable::hz1_iter(:),hz1obs(:)      ! dBz/dt HCP

    real *8,allocatable::time(:)
    real *8,allocatable::deltadobs(:,:),deltam(:,:)
    real *8,allocatable::jacobi(:,:),jacobi1(:,:)
    real *8,allocatable::EYE(:,:),b0(:,:)
    real *8,allocatable::Btemp(:,:), invBtemp(:,:)
    real *8,allocatable::Atemp(:,:),height(:,:),depth(:,:)
    real *8,allocatable::m_pre(:,:),m_next(:,:),Vobs(:,:)

    real *8,allocatable::s(:),e(:),work(:),u(:,:),v(:,:)

    real *8,allocatable::m_true(:,:), m_app(:)

    real *8,allocatable::para_mat(:,:)
    real *8,allocatable::label(:,:), h_mat(:,:)
    real *8,allocatable::true_model(:,:)



    real *8 signal_limit

    Data pi/3.1415926D0/


    real(kind=8) :: F_XC=0.5_RPD, F_CR=0.5_RPD
    integer(kind=IB) :: nfeval
    real(kind=8) :: bestval

    real(kind=8), allocatable :: XCmin(:), XCmax(:)
    real(kind=8), allocatable :: bestmem_XC(:)
    ! integer :: iTimes1, iTimes2, iTimes3, iTimes4, rate
    real :: iTimes1, iTimes2, iTimes3, iTimes4, rate


    external FTN
    open(iwrite,file='devolution.txt')


    snr = 40.0

    ntc = 28           ! 时间道
    nolayer = 5

    npara = 2*nolayer-1  ! 9

    ns = 1             ! 测点个数

    tot_para = npara*ns
    Dim_XC = tot_para
    tot_ntc = ntc*ns
    size0 = (ns-1)*(nolayer-1)
    size1 = (ns-1)*npara

    size2 = tot_ntc                ! 无约束

    ka=max(size2,tot_para)+1

    sigma1 = 0.01            ! 0 大地
    sigma2 = 133              ! 3 混凝土
    sigma3 = 1d-05           ! 2 空气
    sigma4 = 5d-02           ! 1 水


    ! 接收线圈的半径时0.4m
    signal_limit = 1.0d-06/500/pi/0.2/0.2/1600

    allocate(XCmin(Dim_XC),XCmax(Dim_XC),bestmem_XC(Dim_XC))

    allocate(hz1_iter(ntc),hz1obs(ntc),time(ntc))
    allocate(deltadobs(tot_ntc,1),deltam(tot_para,1))
    allocate(jacobi(tot_ntc,tot_para),jacobi1(ntc,npara))
    allocate(EYE(tot_para,tot_para),b0(size2,1))
    allocate(Btemp(tot_para,tot_para), invBtemp(tot_para,tot_para))
    allocate(Atemp(size2,tot_para),height(ns,nolayer-1),depth(ns,nolayer-1))
    allocate(m_pre(tot_para,1),m_next(tot_para,1),Vobs(tot_ntc,1))

    allocate(s(ka))
    allocate(e(ka))
    allocate(work(ka))
    allocate(u(size2,size2))
    allocate(v(tot_para,tot_para))

    allocate(m_true(tot_para,1), m_app(tot_para))
    allocate(para_mat(ns,nolayer),label(ns,nolayer),h_mat(ns,nolayer-1),true_model(ns,npara))

    allocate(rho_true(nolayer), hh_true(nolayer))
    allocate(rho_iter(nolayer), hh_iter(nolayer))

    print*,"size0 = ",size0
    print*,"size1 = ",size1
    print*,"size2 = ",size2
    print*,"tot_ntc = ",tot_ntc
    print*,"tot_para = ",tot_para

    !******************************************************************
    Open (16, File='res2d.dat', Status='unknown')
    Open (7, File='label.dat', Status='unknown')



    ! =================== 水平放置 =========================
    do j=1,nolayer     !计算参数矩阵
        do i=1,ns
            if(i .ge. 1 .and. i .le. 6)then
                if(j .eq. 2 .or. j .eq. 4)then
                    para_mat(i,j)=1.0/sigma2
                    label(i,j)=dlog10(1.0/sigma2)
                else
                    if(j .eq. 3)then
                        para_mat(i,j)=1.0/sigma3
                        label(i,j)=dlog10(1.0/sigma3)
                    else
                        para_mat(i,j)=1.0/sigma1
                        label(i,j)=dlog10(1.0/sigma1)
                    end if
                end if
            else
                para_mat(i,j)=1.0/sigma1
                label(i,j)=dlog10(1.0/sigma1)
            end if
        end do
    end do


    do i=1,ns
        h_mat(i,1) = 50.0
        h_mat(i,2) = 1.0
        h_mat(i,3) = 3.0
        h_mat(i,4) = 1.0

    end do

    true_model = 0.d0
    do j=1,npara     !计算参数矩阵
        do i=1,ns
            if(j .le. nolayer)then
                true_model(i,j) = label(i,j)
            else
                true_model(i,j) = h_mat(i,j-nolayer)
            end if

        end do
    end do

    write(7,20)((true_model(I,J),J=1,npara),I=1,ns)


    !********************* true value of the parameter *****************
    do i=1,ns  ! 测点 从0m出开始
        do j=1,nolayer
            rho_true(j)=para_mat(i,j)

            if(j .lt. nolayer)then
                hh_true(j) = h_mat(i,j)
            end if
        end do

        call CPU_TIME(iTimes1)
        ! forwardprocess(rho, hh, hz1_pls, nlayer,nt,t1,npls)
        call forwardprocess(rho_true,hh_true,hz1obs,nolayer,ntc,time,1) ! observe value
        call CPU_TIME(iTimes2)

        print*, 'Forward Time cost: ', iTimes2-iTimes1

        ! 添加噪声
        do j=1,ntc
            sigma = hz1obs(j)/(10**(snr/20))       ! sigma
            hz1obs(j)=hz1obs(j)+normal(0.d0,1.d0)*sigma
        end do


        !        ! 考虑分辨率
        !        do j=1,ntc
        !            if(abs(hz1obs(j)) .lt. 1.d-06)then
        !                sigma = 1.d-06/(10**(snr/20))       ! sigma
        !                Vobs(j+(i-1)*ntc,1) = 1.d-06+normal(0.d0,1.d0)*sigma
        !            else
        !                sigma = hz1obs(j)/(10**(snr/20))
        !                Vobs(j+(i-1)*ntc,1) = hz1obs(j)+normal(0.d0,1.d0)*sigma
        !            end if
        !        end do

        do j=1,ntc
            Vobs(j+(i-1)*ntc,1) = hz1obs(j)*1.d0
        end do

    end do

    Write(16,10)Vobs

    ! print*,time

    do i=1,nolayer
        hh_iter(i) = 5.0
        rho_iter(i) = 30.0
    end do


    ! apparent rho
    do k=1,ns
        ! 视电阻率成像作为反演的初始值
        if(k .ge. 1 .and. k .le. 6)then
            rho_iter(1) = 6.043
            rho_iter(2) = 0.2
            rho_iter(3) = 0.2
            rho_iter(4) = 0.1
        else
            rho_iter(1) = 30.0
            rho_iter(2) = 30.0
            rho_iter(3) = 30.0
            rho_iter(4) = 30.0
        end if


        ! 当前时刻的参数：电导率
        do i=1,npara
            if(i .le. nolayer)then  ! 电导率的初始值和范围
                m_app(i+(k-1)*npara) = dlog(rho_iter(i))

                XCmin(i+(k-1)*npara) = dlog(rho_iter(i))-4
                XCmax(i+(k-1)*npara) = dlog(rho_iter(i))+4

            else    ! 厚度的初始值和范围
                m_app(i+(k-1)*npara) = dlog(hh_iter(i-nolayer))

                XCmin(i+(k-1)*npara) = dlog(1.d-01)
                XCmax(i+(k-1)*npara) = dlog(1.d02)

            end if
        end do
    end do


    ! DE 反演程序
    ! call DE_Fortran90(FTN, Dim_XC, XCmin, XCmax, VTR, NP, itermax, F_XC,&
    !     CR_XC, strategy, refresh, iwrite, bestmem_XC, bestval, nfeval,&
    !      F_CR, method,tot_ntc,ntc,nolayer,npara,ns,m_app,Vobs)

    ! write(iwrite,205) NP, nfeval, method(1:3)
    ! write(iwrite,FMT=201) F_XC, CR_XC, F_CR
    ! write(iwrite,FMT=200) bestval

    ! do i=1,Dim_XC
    !     write(iwrite,FMT=202) i,dexp(bestmem_XC(i))
    ! end do

200 format(/2x, 'Bestval=', ES14.7)
201 format(2x, 'F_XC =',F6.3, 2x, 'CR_XC =', F6.3, 2x, 'F_CR =', F6.3)
202 format(2x, 'best_XC(',I3,') =',ES14.7)
205 format(2x, 'NP=', I4, 4x, 'No. function call =', I9, &
        /2x, 'mehtod(1:3) =',3I3)


10  Format(28E14.6)

20  Format(9E14.6)   ! 每个测点的参数个数 2*nolayers-1

    end program


    ! ---------------------- objective function -----------------------
    ! X -- parameter
    ! tot_ntc --
    ! Dim_XC -- dimension of parameter
    ! ntc -- number of time
    ! m_ap -- apparent rho
    ! p1,p2,p3,p4
    ! objval
    subroutine FTN(X,tot_ntc,Dim_XC,ntc,nolayer,npara,ns,m_ap,Vobs,objval)
    use data_type

    implicit none
    integer :: ntc,nolayer,npara,ns,tot_ntc,Dim_XC
    real *8 m_ap(Dim_XC)

    integer :: i,k

    real(kind=8), intent(in) :: X(Dim_XC)
    real(kind=8), intent(out) :: objval

    real *8 rho_iter(npara), hh_iter(npara), hz1_iter(ntc), time(ntc)

    real *8 temp0(1,1), temp1, deltam(Dim_XC)

    real *8 deltadobs(tot_ntc,1), Vobs(tot_ntc,1)


    do k=1,ns
        do i=1,npara
            if(i .le. nolayer)then
                rho_iter(i) = dexp(X(i+(k-1)*npara))
            else
                hh_iter(i-nolayer) = dexp(X(i+(k-1)*npara))
            end if
        end do

        call forwardprocess(rho_iter,hh_iter,hz1_iter,nolayer,ntc,time,1)

        ! 计算当前测点的误差
        do i=1,ntc
            deltadobs(i+(k-1)*ntc,1) = Vobs(i+(k-1)*ntc,1)-hz1_iter(i)*1.d0
        end do
    end do

    deltam = X-m_ap


    temp0 = matmul(transpose(deltadobs),deltadobs)
!    temp1 = dot_product(deltam,deltam)
!
!
!    print*,'temp0 = ',1.d15*temp0
!    print*,'temp1 = ',1.d-02*temp1


    objval = temp0(1,1)

    end subroutine FTN



    subroutine DE_Fortran90(obj, Dim_XC, XCmin, XCmax, VTR, NP, itermax, F_XC, &
        CR_XC, strategy, refresh, iwrite, bestmem_XC, bestval, nfeval, &
        F_CR, method,tot_ntc,ntc,nolayer,npara,ns,m_ap,Vobs1)
    !.......................................................................
    !
    ! Differential Evolution for Optimal Control Problems
    !
    !.......................................................................
    !  This Fortran 90 program translates from the original MATLAB
    !  version of differential evolution (DE). This FORTRAN 90 code
    !  has been tested on Compaq Visual Fortran v6.1.
    !  Any users new to the DE are encouraged to read the article of Storn and Price.
    !
    !  Refences:
    !  Storn, R., and Price, K.V., (1996). Minimizing the real function of the
    !    ICEC'96 contest by differential evolution. IEEE conf. on Evolutionary
    !    Comutation, 842-844.
    !
    !  This Fortran 90 program written by Dr. Feng-Sheng Wang
    !  Department of Chemical Engineering, National Chung Cheng University,
    !  Chia-Yi 621, Taiwan, e-mail: chmfsw@ccunix.ccu.edu.tw
    !.........................................................................
    !                obj : The user provided file for evlauting the objective function.
    !                      subroutine obj(xc,fitness)
    !                      where "xc" is the real decision parameter vector.(input)
    !                            "fitness" is the fitness value.(output)
    !             Dim_XC : Dimension of the real decision parameters.
    !      XCmin(Dim_XC) : The lower bound of the real decision parameters.
    !      XCmax(Dim_XC) : The upper bound of the real decision parameters.
    !                VTR : The expected fitness value to reach.
    !                 NP : Population size.
    !            itermax : The maximum number of iteration.
    !               F_XC : Mutation scaling factor for real decision parameters.
    !              CR_XC : Crossover factor for real decision parameters.
    !           strategy : The strategy of the mutation operations is used in HDE.
    !            refresh : The intermediate output will be produced after "refresh"
    !                      iterations. No intermediate output will be produced if
    !                      "refresh < 1".
    !             iwrite : The unit specfier for writing to an external data file.
    ! bestmen_XC(Dim_XC) : The best real decision parameters.
    !              bestval : The best objective function.
    !             nfeval : The number of function call.
    !         method(1) = 0, Fixed mutation scaling factors (F_XC)
    !                   = 1, Random mutation scaling factors F_XC=[0, 1]
    !                   = 2, Random mutation scaling factors F_XC=[-1, 1]
    !         method(2) = 1, Random combined factor (F_CR) used for strategy = 6
    !                        in the mutation operation
    !                   = other, fixed combined factor provided by the user
    !         method(3) = 1, Saving results in a data file.
    !                   = other, displaying results only.

    use data_type, only : IB, RPD
    implicit none

    integer :: ntc,nolayer,npara,ns,tot_ntc,Dim_XC
    real *8 m_ap(Dim_XC),Vobs1(tot_ntc,1)

    integer(kind=IB), intent(in) :: NP, itermax, strategy,   &
        iwrite, refresh
    real(kind=8), intent(in) :: VTR, CR_XC
    real(kind=8) :: F_XC, F_CR
    real(kind=8), dimension(Dim_XC), intent(in) :: XCmin, XCmax
    real(kind=8), dimension(Dim_XC), intent(inout) :: bestmem_XC
    real(kind=8), intent(out) :: bestval
    integer(kind=IB), intent(out) :: nfeval
    real(kind=8), dimension(NP,Dim_XC) :: pop_XC, bm_XC, mui_XC, mpo_XC,   &
        popold_XC, rand_XC, ui_XC
    integer(kind=IB) :: i, ibest, iter
    integer(kind=IB), dimension(NP) :: rot, a1, a2, a3, a4, a5, rt
    integer(kind=IB), dimension(4) :: ind
    real(kind=8) :: tempval
    real(kind=8), dimension(NP) :: val
    real(kind=8), dimension(Dim_XC) :: bestmemit_XC
    real(kind=8), dimension(Dim_XC) :: rand_C1
    integer(kind=IB), dimension(3), intent(in) :: method
    external  obj
    intrinsic max, min, random_number, mod, abs, any, all, maxloc
    integer(kind=IB) :: n,number,y,z




    !!-----Initialize a population --------------------------------------------!!

    pop_XC=0.0_RPD
    do i=1,NP
        call random_number(rand_C1)  ! 0~1
        pop_XC(i,:)=XCmin+rand_C1*(XCmax-XCmin)
    end do

    !!--------------------------------------------------------------------------!!

    !!------Evaluate fitness functions and find the best member-----------------!!
    val=0.0_RPD   ! value of objective function
    nfeval=0
    ibest=1
    call obj(pop_XC(1,:),tot_ntc,Dim_XC,ntc,nolayer,npara,ns,m_ap,Vobs1,val(1))  ! 1st population
    bestval=val(1)
    nfeval=nfeval+1
    do i=2,NP
        call obj(pop_XC(i,:),tot_ntc,Dim_XC,ntc,nolayer,npara,ns,m_ap,Vobs1,val(i))
        nfeval=nfeval+1
        if (val(i) < bestval) then
            ibest=i  ! the best ith population
            bestval=val(i)
        end if
    end do
    bestmemit_XC=pop_XC(ibest,:)
    bestmem_XC=bestmemit_XC

    !!--------------------------------------------------------------------------!!

    bm_XC=0.0_RPD
    rot=(/(i,i=0,NP-1)/)
    iter=1
    !!------Perform evolutionary computation------------------------------------!!

    do while (iter <= itermax)

        print*, iter


        popold_XC=pop_XC

        !!------Mutation operation--------------------------------------------------!!
        n=4
        do z=1,n
            call randperm(z,number)
            ind(z)=number
        enddo

        do y=1,NP
            call randperm(y,number)
            a1(y)=number
        enddo

        ! 随机生成索引
        rt=mod(rot+ind(1),NP)
        a2=a1(rt+1)
        rt=mod(rot+ind(2),NP)
        a3=a2(rt+1)
        rt=mod(rot+ind(3),NP)
        a4=a3(rt+1)
        rt=mod(rot+ind(4),NP)
        a5=a4(rt+1)
        bm_XC=spread(bestmemit_XC, DIM=1, NCOPIES=NP)  ! 按照维度dim复制bestmemit_XC矩阵NP次

        !----- Generating a random sacling factor--------------------------------!
        select case (method(1))
        case(1)
            call random_number(F_XC)  ! 0~1
        case(2)
            call random_number(F_XC)
            F_XC=2.0_RPD*F_XC-1.0_RPD   !-1~1
        end select

        !---- select a mutation strategy-----------------------------------------!
        select case (strategy)
        case (1)
            ui_XC=bm_XC+F_XC*(popold_XC(a1,:)-popold_XC(a2,:))

            case default
            ui_XC=popold_XC(a3,:)+F_XC*(popold_XC(a1,:)-popold_XC(a2,:))

        case (3)
            ui_XC=popold_XC+F_XC*(bm_XC-popold_XC+popold_XC(a1,:)-popold_XC(a2,:))

        case (4)
            ui_XC=bm_XC+F_XC*(popold_XC(a1,:)-popold_XC(a2,:)+popold_XC(a3,:)-popold_XC(a4,:))

        case (5)
            ui_XC=popold_XC(a5,:)+F_XC*(popold_XC(a1,:)-popold_XC(a2,:)+popold_XC(a3,:) &
                -popold_XC(a4,:))
        case (6) ! A linear crossover combination of bm_XC and popold_XC
            if (method(2) == 1) call random_number(F_CR)
            ui_XC=popold_XC+F_CR*(bm_XC-popold_XC)+F_XC*(popold_XC(a1,:)-popold_XC(a2,:))

        end select
        !!--------------------------------------------------------------------------!!
        !!------Crossover operation-------------------------------------------------!!
        call random_number(rand_XC)
        mui_XC=0.0_RPD
        mpo_XC=0.0_RPD
        where (rand_XC < CR_XC)
            mui_XC=1.0_RPD
            !           mpo_XC=0.0_RPD
        elsewhere
            !           mui_XC=0.0_RPD
            mpo_XC=1.0_RPD
        end where

        ui_XC=popold_XC*mpo_XC+ui_XC*mui_XC
        !!--------------------------------------------------------------------------!!
        !!------Evaluate fitness functions and find the best member-----------------!!
        do i=1,NP
            !!------Confine each of feasible individuals in the lower-upper bound-------!!
            ui_XC(i,:)=max(min(ui_XC(i,:),XCmax),XCmin)
            call obj(ui_XC(i,:),tot_ntc,Dim_XC,ntc,nolayer,npara,ns,m_ap,Vobs1,tempval)
            nfeval=nfeval+1
            if (tempval < val(i)) then
                pop_XC(i,:)=ui_XC(i,:)
                val(i)=tempval
                if (tempval < bestval) then
                    bestval=tempval
                    bestmem_XC=ui_XC(i,:)
                end if
            end if
        end do
        bestmemit_XC=bestmem_XC
        if( (refresh > 0) .and. (mod(iter,refresh)==0)) then
!            if (method(3)==1) write(unit=iwrite,FMT=203) iter
            write(unit=*, FMT=203) iter
            do i=1,Dim_XC
!                if (method(3)==1) write(unit=iwrite, FMT=202) i, bestmem_XC(i)
                write(*,FMT=202) i,dexp(bestmem_XC(i))
            end do
            if (method(3)==1) write(unit=iwrite, FMT=201) bestval
            write(unit=*, FMT=201) bestval
        end if
        iter=iter+1
        if ( bestval <= VTR .and. refresh > 0) then
            write(unit=iwrite, FMT=*) ' The best fitness is smaller than VTR'
            write(unit=*, FMT=*) 'The best fitness is smaller than VTR'
            exit
        endif

        print*,bestval
    end do
    !!------end the evolutionary computation------------------------------!!
201 format(2x, 'bestval =', ES14.7, /)
202 format(5x, 'bestmem_XC(', I3, ') =', ES12.5)
203 format(2x, 'No. of iteration  =', I8)
    end subroutine DE_Fortran90


    subroutine randperm(num,number)
    use data_type, only : IB, RPD
    implicit none
    integer(kind=IB) :: num
    integer(kind=IB) :: number, i, j, k
    real(kind=8), dimension(num) :: rand2
    intrinsic random_number
    call random_number(rand2)
    do i=1,num
        number=1
        do j=1,num
            if (rand2(i) > rand2(j)) then
                number=number+1
            end if
        end do
        do k=1,i-1
            if (rand2(i) <= rand2(k) .and. rand2(i) >= rand2(k)) then
                number=number+1
            end if
        end do
    end do
    return
    end



    Subroutine forwardprocess(rho, hh, hz1_pls, nlayer,nt,t1,npls)
    integer::nlayer, nt,npls
    Real *8 t_st,t_ed,delta_t
    Real *8 t(npls*nt), hz1(npls*nt),tlog(npls*nt)

    Real *8 hz1_pls(nt),t1(nt)

    Real *8 s2, ss2
    Real *8 r, rplus, tt1, tt2, tt3, ti0, pi, tm
    Real *8 rho(nlayer), hh(nlayer), frq(67)

    Real *8 xr, yr, rt, ht, hr, zplus, zminus
    Complex *16 func(5, 67)

    Data pi/3.1415926D0/
    Data ngau/1000/
    Common /para/r

    Common /funn/frq, func

    !**************************************************************
    Call filter
    !**************************************************************
    ht = 1.
    hr = 1.
    xr = 0.58
    yr = 0.

    zplus = ht - hr
    zminus = ht + hr
    r = dsqrt(xr*xr+yr*yr)
    rplus = dsqrt(r*r+zplus*zplus)

    !**************************************************************
    ! 发射线圈
    rt = 0.5          ! 线圈的半径
    nturn = 3.         ! 线圈的匝数
    tm = 200
    ti0 = tm/(4*rt*rt*nturn)  ! 在达到磁矩时的电流


    !**************************************************************

    t_st = 2.04e-3
    t_ed = 20e-3

    delta_t = (dlog10(t_ed)-dlog10(t_st))/(nt-1)

    do k=1,npls
        do i=(nt*(k-1)+1),(nt*k)
            if(k==1)then                         ! 第一个脉冲
                tlog(i)=dlog10(t_st)+(i-1)*delta_t
                t(i)=10**tlog(i)
            else
                t(i)=t(i-nt)+t_ed            ! 这里加的t_ed是假设t_ed此时是周期
            end if
        end do
    end do

    t1 = t(1:nt)


    ic = 2


    If (ic==2) Then   ! 单极性方波
        tt1 = 2e-3
        tt2 = t_ed-tt1
    End If



    !********************* impulse and step wave ************************
    If (ic==0 .Or. ic==1) Then
        ik = 0
        Do i = 1,size(t)

            Call frt(rho, hh, t(i), hz1(i), 2, zplus, zminus, ic, ik, nlayer)

            ik = ik + 1

            ! Transformation of cylindrical coordinate system
            hz1(i) = tm*hz1(i)

        End Do
        !******************* for rectangle wave(multi-pulse)***********************
    Else If (ic==2) Then
        Do i = 1, size(t)

            ss2 = 0.D0

            ik = 0
            Do ip = 1, npls
                If (t(i)>(ip-1)*(tt1+tt2) .And. t(i)<ip*(tt1+tt2)) kpls = ip   ! 判断抽道时间位于哪一个脉冲段
            End Do

            Do ip = 1,kpls-1

                Call frt(rho, hh, t(i)-(ip-1)*(tt1+tt2), s2, 2, zplus, zminus, 0, ik, nlayer)  ! 脉冲响应

                ss2 = ss2 + s2      ! positive delta


                ik = ik + 1

                Call frt(rho, hh, t(i)-(ip-1)*(tt1+tt2)-tt1, s2, 2, zplus, zminus, 0, ik, nlayer)

                ss2 = ss2 + (-1)*s2        ! negative delta

            End Do

            !****************** contribution from resting pulses**************************
            If (t(i)>(kpls-1)*(tt1+tt2))Then

                Call frt(rho, hh, t(i)-(kpls-1)*(tt1+tt2), s2, 2, zplus, zminus, 0, ik, nlayer)

                ss2 = ss2 + s2

                ik = ik + 1

                If (t(i)>(kpls-1)*(tt1+tt2)+tt1)Then   ! second delta

                    Call frt(rho, hh, t(i)-(kpls-1)*(tt1+tt2)-tt1, s2, 2, zplus, zminus, 0, ik, nlayer)

                    ss2 = ss2 + (-1)*s2

                End If
            End If

            hz1(i) = -tm*ss2

        End Do
    End If


    ! 抽道叠加
    hz1_pls = 0.d0

    do i=1,nt
        do k=1,npls
            hz1_pls(i) = hz1_pls(i)+hz1(i+(k-1)*nt)/npls
        end do
    end do


    End Subroutine forwardprocess








    Subroutine gauleg(x1, x2, x, w, n)
    Implicit Real *8(A-H, O-Z)
    Real *8 x1, x2, xm, xl, x(n), w(n)
    Parameter (eps=3.D-14)
    m = (n+1)/2
    xm = 0.5D0*(x2+x1)
    xl = 0.5D0*(x2-x1)
    Do i = 1, m
        z = cos(3.141592654D0*(i-0.25D0)/(n+0.5D0))
1       Continue
        p1 = 1.D0
        p2 = 0.D0
        Do j = 1, n
            p3 = p2
            p2 = p1
            p1 = ((2.D0*j-1.D0)*z*p2-(j-1.D0)*p3)/j
        End Do
        pp = n*(z*p1-p2)/(z*z-1.D0)
        z1 = z
        z = z1 - p1/pp
        If (abs(z-z1)>eps) Goto 1
        x(i) = xm - xl*z
        x(n+1-i) = xm + xl*z
        w(i) = 2.D0*xl/((1.D0-z*z)*pp*pp)
        w(n+1-i) = w(i)
    End Do
    Return
    End Subroutine gauleg

    !        spl(nfrq, nc, frq, funr0, f, funr1)
    Subroutine spl(nx, n2, x, fx, x2, fx2)
    Real *8 x(nx), fx(nx), c(3, nx), x2(n2), fx2(n2), xint, xl1, xl2
    Call splin1(nx, fx, c)
    xl1 = dlog10(x(1))
    xl2 = dlog10(x(nx))
    Do ix = 1, n2
        xint = dlog10(x2(ix))
        Call splin2(nx, xint, xl1, xl2, c, fx2(ix))
    End Do
    End Subroutine spl

    Subroutine splin1(n, y, c)
    Real *8 y(n), c(3, n), p
    n1 = n - 1
    Do i = 2, n1
        c(1, i) = y(i+1) - 2.*y(i) + y(i-1)
    End Do
    c(2, 1) = 0.D0
    c(3, 1) = 0.D0
    Do i = 2, n1
        p = 4. + c(2, i-1)
        c(2, i) = -1.D0/p
        c(3, i) = (c(1,i)-c(3,i-1))/p
    End Do
    c(1, n) = 0.D0
    Do ii = 2, n1
        i = n + 1 - ii
        c(1, i) = c(2, i)*c(1, i+1) + c(3, i)
    End Do
    c(1, 1) = 0.
    Do i = 1, n1
        c(2, i) = y(i+1) - y(i) - c(1, i+1) + c(1, i)
        c(3, i) = y(i) - c(1, i)
    End Do
    c(3, n) = y(n)
    Return
    End Subroutine splin1

    Subroutine splin2(n, xint, x1, x2, c, yint)
    Real *8 c(3, n), xint, x1, x2, yint, h, u, p, q
    h = (x2-x1)/dble(float(n-1))
    If (xint<x1) Goto 10
    If (xint>=x2) Goto 20
    u = (xint-x1)/h
    i = 1 + int(u)
    p = u - i + 1
    q = 1.D0 - p
    yint = c(1, i)*q**3 + c(1, i+1)*p**3 + c(2, i)*p + c(3, i)
    Return
10  p = (xint-x1)/h
    yint = c(2, 1)*p + c(3, 1)
    Return
20  p = (xint-x2)/h
    yint = c(2, n-1)*p + c(3, n)
    Return
    End Subroutine splin2

    Subroutine frt(rho1, hh1, t, ft, item, zplus, zminus, ic, ik, nlayer)
    Complex *16 fun, iomega, func(5, 67)
    Real *8 t, ft, zplus, zminus, pi, q, rho1(nlayer), hh1(nlayer)
    Real *8 frq(67), funr0(67), funi0(67)
    Real *8 f(160), omega(160), funr1(160), funi1(160), h(200)

    ! integer :: iTimes1, iTimes2, iTimes3, iTimes4, rate
    real :: iTimes1, iTimes2, iTimes3, iTimes4, rate

    Common /funn/frq, func
    Data pi, q/3.141592654D0, 1.258925412D0/
    Data ncnull, nc, ndec, (h(i), i=1, 160)/80,160,10,&
        & 2.59511139938829d-13,3.66568771323555d-13,5.17792876616242d-13,&
        & 7.31400730405791d-13,1.03313281156235d-12,1.45933600088387d-12,&
        & 2.06137146234699d-12,2.91175733962418d-12,4.11297804457870d-12,&
        & 5.80971771117984d-12,8.20647323099742d-12,1.15919058389365d-11,&
        & 1.63740746547780d-11,2.31288803930431d-11,3.26705938902288d-11,&
        & 4.61481520721098d-11,6.51864545047052d-11,9.20775899532545d-11,&
        & 1.30064200980219d-10,1.83718747396255d-10,2.59512512377884d-10,&
        & 3.66566596154242d-10,5.17796324027279d-10,7.31395266627501d-10,&
        & 1.03314147106736d-09,1.45932227649333d-09,2.06139321404013d-09,&
        & 2.91172286551380d-09,4.11303268236158d-09,5.80963111612975d-09,&
        & 8.20661047490285d-09,1.15916883220051d-08,1.63744193958818d-08,&
        & 2.31283340152144d-08,3.26714598407299d-08,4.61467796330556d-08,&
        & 6.84744728867720d-08,5.46574677490374d-08,1.13319898777493d-07,&
        & 2.16529974157527d-07,2.88629942214140d-07,3.42872728051125d-07,&
        & 4.79119488706262d-07,7.42089418889752d-07,1.07736520535271d-06,&
        & 1.46383231306575d-06,2.01727682134668d-06,2.89058197617431d-06,&
        & 4.15237808867022d-06,5.84448989361742d-06,8.18029430348419d-06,&
        & 1.15420854481494d-05,1.63897017145322d-05,2.31769096113890d-05,&
        & 3.26872676331330d-05,4.60786866701851d-05,6.51827321351636d-05,&
        & 9.20862589540037d-05,1.30169142615951d-04,1.83587481111627d-04,&
        & 2.59595544393723d-04,3.66324383719323d-04,5.18210697462501d-04,&
        & 7.30729969562531d-04,1.03385239132389d-03,1.45738764044730d-03,&
        & 2.06298256402732d-03,2.90606401578959d-03,4.11467957883740d-03,&
        & 5.79034253321120d-03,8.20005721235220d-03,1.15193892333104d-02,&
        & 1.63039398900789d-02,2.28256810984487d-02,3.22248555163692d-02,&
        & 4.47865101670011d-02,6.27330674874545d-02,8.57058672847471d-02,&
        & 1.17418179407605d-01,1.53632645832305d-01,1.97718111895102d-01,&
        & 2.28849924263247d-01,2.40310905012422d-01,1.65409071929404d-01,&
        & 2.84709685167114d-03,-2.88015846269687d-01,-3.69097391853225d-01,&
        & -2.50109865922601d-02,5.71811109500426d-01,-3.92261390212769d-01,&
        & 7.63282774297327d-02,5.16233692927851d-02,-6.48015160576432d-02,&
        & 4.89045522502552d-02,-3.26934307794750d-02,2.10542570949745d-02,&
        & -1.33862848934736d-02,8.47098801479259d-03,-5.35134515919751d-03,&
        & 3.37814023806349d-03,-2.13157364002470d-03,1.34506352474558d-03,&
        & -8.48929743771803d-04,5.35521822356713d-04,-3.37744799986382d-04,&
        & 2.13268792633204d-04,-1.34629969723156d-04,8.47737416679279d-05,&
        & -5.34940635827096d-05,3.3904416298191d-05,-2.13315638358794d-05,&
        & 1.33440911625019d-05,-8.51629073825634d-06,5.44362672273211d-06,&
        & -3.32112278417896d-06,2.07147190852386d-06,-1.42009412555511d-06,&
        & 8.78247754998004d-07,-4.5566280473703d-07,3.38598103040009d-07,&
        & -2.87407830772251d-07,1.07866150545699d-07,-2.4724024185358d-08,&
        & 5.35535110396030d-08,-3.3789981131378d-08,2.13200367531820d-08,&
        & -1.34520337740075d-08,8.48765950790546d-09,-5.35535110396018d-09,&
        & 3.37899811131383d-09,-2.13200367531819d-09,1.34520337740075d-09,&
        & -8.48765950790576d-10,5.35535110396015d-10,-3.37899811131382d-10,&
        & 2.13200367531811d-10,-1.34520337740079d-10,8.48765950790572d-11,&
        & -5.35535110396034d-11,3.37899811131381d-11,-2.13200367531818d-11,&
        & 1.34520337740074d-11,-8.48765950790571d-12,5.35535110396031d-12,&
        & -3.37899811131379d-12,2.13200367531817d-12,-1.34520337740073d-12,&
        & 8.48765950790567d-13,-5.35535110396029d-13,3.37899811131377d-13,&
        & -2.13200367531816d-13,1.34520337740078d-13,-8.48765950790596d-14,&
        & 5.35535110396007d-14,-3.37899811131377d-14,2.13200367531816d-14,&
        & -1.34520337740083d-14,8.4876550790558d-15,-5.35535110396025d-15,&
        & 3.37899811131389d-15/
    Data nfrq, (frq(i), i=1, 67)/67,&
        & 0.10000000d-02,0.14677993d-02,0.21544347d-02,0.31622777d-02,&
        & 0.46415888d-02,0.68129207d-02,0.10000000d-01,0.14677993d-01,&
        & 0.21544347d-01,0.31622777d-01,0.46415888d-01,0.68129207d-01,&
        & 0.10000000d+00,0.14677993d+00,0.21544347d+00,0.31622777d+00,&
        & 0.46415888d+00,0.68129207d+00,0.10000000d+01,0.14677993d+01,&
        & 0.21544347d+01,0.31622777d+01,0.46415888d+01,0.68129207d+01,&
        & 0.10000000d+02,0.14677993d+02,0.21544347d+02,0.31622777d+02,&
        & 0.46415888d+02,0.68129207d+02,0.10000000d+03,0.14677993d+03,&
        & 0.21544347d+03,0.31622777d+03,0.46415888d+03,0.68129207d+03,&
        & 0.10000000d+04,0.14677993d+04,0.21544347d+04,0.31622777d+04,&
        & 0.46415888d+04,0.68129207d+04,0.10000000d+05,0.14677993d+05,&
        & 0.21544347d+05,0.31622777d+05,0.46415888d+05,0.68129207d+05,&
        & 0.10000000d+06,0.14677993d+06,0.21544347d+06,0.31622777d+06,&
        & 0.46415888d+06,0.68129207d+06,0.10000000d+07,0.14677993d+07,&
        & 0.21544347d+07,0.31622777d+07,0.46415888d+07,0.68129207d+07,&
        & 0.10000000d+08,0.14677993d+08,0.21544347d+08,0.31622777d+08,&
        & 0.46415888d+08,0.68129207d+08,0.1000000d+09/
    
    call SYSTEM_CLOCK(count_rate=rate)
    
    call CPU_TIME(iTimes1)

    If (ik==0) Then
        Do i = 1, nfrq
            Call forward(rho1, hh1, frq(i), func(item,i), item, zplus, zminus, nlayer) ! mu0*H = B
        End Do
    End If

    call CPU_TIME(iTimes2)

    Do nn = 1, nc
        n = -nc + ncnull + nn
        omega(nn) = q**(-(n-1))/t
        f(nn) = omega(nn)/(2.D0*pi)
    End Do

    Do i = 1, nfrq
        funr0(i) = dreal(func(item,i))
        funi0(i) = dimag(func(item,i))
    End Do
    Call spl(nfrq, nc, frq, funr0, f, funr1)  !interpolation
    Call spl(nfrq, nc, frq, funi0, f, funi1)
    ft = 0.D0
    Do nn = 1, nc
        If (ic==0) Then
        iomega = (1.D0, 0.D0)
    Else If (ic==1) Then
        iomega = 1.D0/((0.,-1.D0)*omega(nn)) ! -1/iw
    End If
    fun = dcmplx(funr1(nn), funi1(nn))*iomega
    ita = max0(1, nn-nc+1)
    ite = min0(1, nn)
    Do it = ita, ite
        itn = nc - nn + it
        ft = ft + dimag(fun)*dsqrt(omega(nn))*h(itn) ! primary field st
    End Do
    End Do
    ft = -ft*dsqrt(2.D0/pi/t)

    call CPU_TIME(iTimes3)

    ! print*, '========================='
    ! print*, 'Forward Time cost: ', real(iTimes2-iTimes1)/real(rate)
    ! print*, 'Other Time cost: ', real(iTimes3-iTimes2)/real(rate)

    ! print*, 'Forward Time cost: ', iTimes2-iTimes1
    ! print*, 'Other Time cost: ', iTimes3-iTimes2
    ! print*, '========================='

    Return

    End Subroutine frt

    Subroutine forward(rho2, hh2,  f, fun, item, zplus, zminus, nlayer)
    Complex *16 t3, t5, t6, hf, fun
    Real *8 pi, f, zplus, zminus
    Real *8 r, rho2(nlayer), hh2(nlayer)
    Common /para/r
    pi = 3.1415926D0
    If (item==1)Then
        hf = t6(rho2, hh2,f, zminus, nlayer)/(4.D0*pi)
    Else If (item==2)Then
        hf = -t3(rho2, hh2,f, zminus, nlayer)/(4.D0*pi)
    Else If (item==3)Then
        hf = (-t3(rho2, hh2,f,zminus, nlayer)+t5(rho2, hh2,f,zminus, nlayer)/r)/(4.D0*pi)
        If (zplus<0.D0) hf = -hf
    Else If (item==4)Then
        hf = t5(rho2, hh2,f, zminus, nlayer)/(4.D0*pi*r)
        If (zplus<0.D0) hf = -hf
    Else If (item==5) Then
        hf = -t6(rho2, hh2,f, zminus, nlayer)/(4.D0*pi)
        If (zplus<0.D0) hf = -hf
    Else
        Print *, 'item must be between 1 and 5'
    End If
    fun = hf*4.D-7*pi
    Return
    End Subroutine forward


    Complex *16 Function t3(rho3, hh3, f, z, nlayer)
    Complex *16 s, s1, b
    Real *8 r, rho3(nlayer), hh3(nlayer)
    Real *8 h0, h1, f, z, u, fac, expc
    Common /para/r
    Common /hankel/nc, ncnull, h0(100), h1(100)
    fac = 0.1D0*dlog(10.D0)
    s = (0.D0, 0.D0)
    Do nn = 1, nc
        nu = nn
        mn = nc - nn + 1
        nnn = ncnull - nc + nu
        u = expc(-(nnn-1)*fac)/r
        s1 = (b(rho3, hh3,f,u, nlayer)-u)/(b(rho3, hh3,f,u, nlayer)+u)*expc(-u*z)*u*u
        s = s + s1*h0(mn)
    End Do
    t3 = s/r
    Return
    End Function t3

    Complex *16 Function t5(rho5, hh5, f, z, nlayer)
    Complex *16 s, s1, b
    Real *8 r, rho5(nlayer), hh5(nlayer)
    Real *8 h0, h1, f, z, u, fac, expc
    Common /para/r
    Common /hankel/nc, ncnull, h0(100), h1(100)
    fac = 0.1D0*dlog(10.D0)
    s = (0.D0, 0.D0)
    Do nn = 1, nc
        nu = nn
        mn = nc - nn + 1
        nnn = ncnull - nc + nu
        u = expc(-(nnn-1)*fac)/r
        s1 = (b(rho5, hh5,f,u, nlayer)-u)/(b(rho5, hh5,f,u, nlayer)+u)*expc(-u*z)*u
        s = s + s1*h1(mn)
    End Do
    t5 = s/r
    Return
    End Function t5

    Complex *16 Function t6(rho6, hh6, f, z, nlayer)
    Complex *16 s, s1, b
    Real *8 r, rho6(nlayer), hh6(nlayer)
    Real *8 h0, h1, f, z, u, fac, expc
    Common /para/r
    Common /hankel/nc, ncnull, h0(100), h1(100)
    fac = 0.1D0*dlog(10.D0)
    s = (0.D0, 0.D0)
    Do nn = 1, nc
        nu = nn
        mn = nc - nn + 1
        nnn = ncnull - nc + nu
        u = expc(-(nnn-1)*fac)/r
        s1 = (b(rho6, hh6, f,u, nlayer)-u)/(b(rho6, hh6, f,u, nlayer)+u)*expc(-u*z)*u*u
        s = s + s1*h1(mn)
    End Do
    t6 = s/r
    Return
    End Function t6

    ! output: BE
    ! input: f - frequency
    ! input: u - k
    ! input: mu - 磁导率
    ! nlayer -  no. of layers
    Complex *16 Function b(rho4, hh4, f, u, nlayer)
    Complex *16 alpha, s1, s2
    Real *8 f, u, pi,gam
    Real *8 r, rho4(nlayer), hh4(nlayer)
    real *8 mu0
    Common /para/r

    pi = 3.1415926D0
    mu0 = pi*4d-07

    b = cdsqrt(u*u+(0.D0,1.D0)*mu0*2*pi*f/rho4(nlayer))  ! alpha nlayer omega = 2*pi*f
    If (nlayer==1) Return
    Do i = nlayer-1, 1, -1
        gam = 1
        alpha = cdsqrt(u*u+(0.D0,1.D0)*mu0*2*pi*f/rho4(i))
        s1 = (0.D0, 0.D0)
        If (dreal(2.D0*alpha*hh4(i))<400.D0) s1 = cdexp(-2.D0*alpha*hh4(i))
        s2 = (1.D0-s1)/(1.D0+s1)  ! tanh
        b = alpha*(b+alpha*s2*gam)/(alpha*gam+b*s2)
    End Do
    b = mu0/mu0*b
    End Function b

    Real *8 Function expc(x)
    Real *8 x, x1
    x1 = x
    If (dabs(x1)>650.D0) x1 = dsign(650.D0, x1)
    expc = dexp(x1)
    Return
    End Function expc

    Subroutine filter
    Common /hankel/nc, ncnull, h0(100), h1(100)
    Real *8 h0, h1

    Data(h0(i),i=1,48)/&
        & 2.89878288d-07,3.64935144d-07,4.59426126d-07,5.78383226d-07,&
        & 7.28141338d-07,9.16675639d-07,1.15402625d-06,1.45283298d-06,&
        & 1.82900834d-06,2.30258511d-06,2.89878286d-06,3.64935148d-06,&
        & 4.59426119d-06,5.78383236d-06,7.28141322d-06,9.16675664d-06,&
        & 1.15402621d-05,1.45283305d-05,1.82900824d-05,2.30258527d-05,&
        & 2.89878259d-05,3.64935186d-05,4.59426051d-05,5.78383329d-05,&
        & 7.28141144d-05,9.16675882d-05,1.15402573d-04,1.45283354d-04,&
        & 1.82900694d-04,2.30258630d-04,2.89877891d-04,3.64935362d-04,&
        & 4.59424960d-04,5.78383437d-04,7.28137738d-04,9.16674828d-04,&
        & 1.15401453d-03,1.45282561d-03,1.82896826d-03,2.30254535d-03,&
        & 2.89863979d-03,3.64916703d-03,4.59373308d-03,5.78303238d-03,&
        & 7.27941497d-03,9.16340705d-03,1.15325691d-02,1.45145832d-02/

    Data(h0(i),i=49,100)/&
        & 1.82601199d-02,2.29701042d-02,2.88702619d-02,3.62691810d-02,&
        & 4.54794031d-02,5.69408192d-02,7.09873072d-02,8.80995426d-02,&
        & 1.08223889d-01,1.31250483d-01,1.55055715d-01,1.76371506d-01,&
        & 1.85627738d-01,1.69778044d-01,1.03405245d-01,-3.02583233d-02,&
        & -2.27574393d-01,-3.62173217d-01,-2.05500446d-01,3.37394873d-01,&
        & 3.17689897d-01,-5.13762160d-01,3.09130264d-01,-1.26757592d-01,&
        & 4.61967890d-02,-1.80968674d-02,8.35426050d-03,-4.47368304d-03,&
        & 2.61974783d-03,-1.60171357d-03,9.97717882d-04,-6.26275815d-04,&
        & 3.94338818d-04,-2.48606354d-04,1.56808604d-04,-9.89266288d-05,&
        & 6.24152398d-05,-3.93805393d-05,2.48472358d-05,-1.56774945d-05,&
        & 9.89181741d-06,-6.24131160d-06,3.93800058d-06,-2.48471018d-06,&
        & 1.56774609d-06,-9.89180896d-07,6.24130948d-07,-3.93800005d-07,&
        & 2.48471005d-07,-1.56774605d-07,9.89180888d-08,-6.24130946d-08/

    Data(h1(i),i=1,48)/&
        & 1.84909557d-13,2.85321327d-13,4.64471808d-13,7.16694771d-13,&
        & 1.16670043d-12,1.80025587d-12,2.93061898d-12,4.52203829d-12,&
        & 7.36138206d-12,1.13588466d-11,1.84909557d-11,2.85321327d-11,&
        & 4.64471808d-11,7.16694771d-11,1.16670043d-10,1.80025587d-10,&
        & 2.93061898d-10,4.52203829d-10,7.36138206d-10,1.13588466d-09,&
        & 1.84909557d-09,2.85321326d-09,4.64471806d-09,7.16694765d-09,&
        & 1.16670042d-08,1.80025583d-08,2.93061889d-08,4.52203807d-08,&
        & 7.36138149d-08,1.13588452d-07,1.84909521d-07,2.85321237d-07,&
        & 4.64471580d-07,7.16694198d-07,1.16669899d-06,1.80025226d-06,&
        & 2.93060990d-06,4.52201549d-06,7.36132477d-06,1.13587027d-05,&
        & 1.84905942d-05,2.85312247d-05,4.64449000d-05,7.16637480d-05,&
        & 1.16655653d-04,1.79989440d-04,2.92971106d-04,4.51975783d-04/

    Data(h1(i),i=49,100)/&
        & 7.35565435d-04,1.13444615d-03,1.84548306d-03,2.84414257d-03,&
        & 4.62194743d-03,7.10980590d-03,1.15236911d-02,1.76434485d-02,&
        & 2.84076233d-02,4.29770596d-02,6.80332569d-02,9.97845929d-02,&
        & 1.51070544d-01,2.03540581d-01,2.71235377d-01,2.76073871d-01,&
        & 2.16691977d-01,-7.83723737d-02,-3.40675627d-01,-3.60693673d-01,&
        & 5.13024526d-01,-5.94724729d-02,-1.95117123d-01,1.99235600d-01,&
        & -1.38521553d-01,8.79320859d-02,-5.50697146d-02,3.45637848d-02,&
        & -2.17527180d-02,1.37100291d-02,-8.64656417d-03,5.45462758d-03,&
        & -3.44138864d-03,2.17130686d-03,-1.36998628d-03,8.64398952d-04,&
        & -5.45397874d-04,3.44122545d-04,-2.17126585d-04,1.36997597d-04,&
        & -8.64396364d-05,5.45397224d-05,-3.44122382d-05,2.17126544d-05,&
        & -1.36997587d-05,8.64396338d-06,-5.45397218d-06,3.44122380d-06,&
        & -2.17126543d-06,1.36997587d-06,-8.64396337d-07,5.45397218d-07/

    nc = 100
    ncnull = 60
    Return
    End Subroutine filter
