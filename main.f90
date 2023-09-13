program main
    use mod_autodiff
    use mod_condition
    use mod_cubic_eq_d
    use mod_cubic_eq
    use mod_fugacity
    use mod_gauss
    use mod_ini_chemical
    use mod_ini_flash
    use mod_ini_N
    use mod_input
    use mod_phase_stability_analysis

    implicit none

    !!初期計算用
    integer::i,j,iteration,time,phase_judge0,phase0
    real(8),allocatable,dimension(:,:)::amat,cmat,emat
    real(8),allocatable,dimension(:)::bmat,z0,k0,lnk0,x0,y0,w0,alpha0,dmat,fmat,theta0
    real(8)::V0,L0,P0,P0old,error,wt,lumda,z_factor0,fai0
    type(diffs),allocatable::fxs(:)

    !!
    real(8),dimension(com_2phase+com_ion)::Nc0,Nc0old
    real(8),dimension(com_mine)::Nm0,Nm0old
    real(8)::Nt,Sw0,lnk1,lnk2,lnk3,lnk4
    real(8),allocatable,dimension(:,:)::chemi_mat


    call file_open

    allocate(amat(com_2phase,com_2phase),bmat(com_2phase),z0(com_2phase+com_ion),k0(com_2phase),lnk0(com_2phase)&
        ,x0(com_2phase),y0(com_2phase+com_ion),w0(com_2phase),alpha0(com_2phase),chemi_mat(chemi+mine,com_all)&
        ,theta0(chemi+mine))


    !!化学反応の係数のインプット
    do i=1,chemi+mine
        read(1,*) (chemi_mat(i,j),j=1,com_all)
    end do
    
    !初期設定
    
    
    P0 =iniPressure
    
    
    !平衡定数の初期予測
    lnk0(1) = -6.6318705003566869!-5.0d0
    lnk0(2) = 3.4598654370079327  !3.0d0
    lnk0(3) = 5.0d0
    lnk0(4) = 3.0d0
    k0(1) = exp(lnk0(1))!1.0d0/Pc(1)*exp(5.37*(1.0d0+acentric(1)*(1.0d0-1.0d0/tc(1))))!exp(lnk0(1))!=1.0d0/pr(i)*exp(5.37*(1.0d0+omega(i)*(1.0d0-1.0d0/tr(i)))
    k0(2) = exp(lnk0(2))!1.0d0/Pc(2)*exp(5.37*(1.0d0+acentric(2)*(1.0d0-1.0d0/tc(2))))!exp(lnk0(2))
    k0(3) = exp(lnk0(3))
    k0(4) = exp(lnk0(4))
    
    
    !初期の全モル分率
    !z0(1) = 0.95d0 !H2O
    !z0(2) = 0.01d0 !CO2
    !z0(3) = 0.0d0 !CH4
    !z0(4) = 0.01d0 !H2S
    do i =com_2phase+1,com_2phase+com_ion
    !    z0(i) = (1.0d0 - (z0(1)+z0(2)+z0(3)+z0(4)))/(com_2phase+com_ion-com_2phase) !そのほかのイオン
    end do
    
    !!!検証 これ
    z0=0.0000d0
    z0(1)=0.982!0.98295661537557144!0.98294928708382090!0.98294  !H2O
    z0(2)=2.34d0*10.0d0**(-7.0d0)!3.4817505806941896E-007!1.8992194875568811E-005!0.00001 !CO2
    z0(3)=1.70d0*10.0d0**(-2.0d0)!1.6990282021069454E-002!1.6990001958143271E-002!0.01699 !CH4
    z0(5)=6.57d0*10.0d0**(-11.0d0)!6.1817158796360923E-011!6.3674605894641180E-007!1E-5 !H+
    z0(6)=1.76d0*10.0d0**(-5.0d0)!2.8331115169197601E-005!1.1005716001412855E-006!1E-5 !HCO3-
    z0(7)=2.75d0*10.0d0**(-7.0d0)!1.3212068937893931E-006!9.9072355493762180E-006!1E-5 !CO32-
    z0(8)=2.60d0*10.0d0**(-7.0d0)!1.0885035523759696E-006!9.9072356593756316E-006!1E-5 !OH-
    z0(10)=4.23d0*10.0d0**(-6.0d0)!1.1030901971414173E-005!1.0092765812147074E-005!1E-5 !Mg2+
    z0(11)=1.0d0-z0(1)-z0(2)-z0(3)-z0(5)-z0(6)-z0(7)-z0(8)-z0(10) !SiO2
    
    !!学部ケースこれ
    !z0=0.0d0
    !z0(1)=0.990029871
    !z0(2)=1.68415E-05
    !z0(3)=0.008000107!0.0d0
    !z0(4)=0.0d0
    !z0(5)=9.43763E-10
    !z0(6)=0.000118031
    !z0(7)=1.3007E-07
    !z0(8)=1.86134E-08
    !z0(10)=9.13985E-06
    !z0(11)=1.0d0-z0(1)-z0(2)-z0(3)-z0(5)-z0(6)-z0(7)-z0(8)-z0(10) !SiO2

    !?相安定性解析ケース
    z0=0.0d0
    z0(1) = 0.999d0 !H2O
    z0(2) = 2.38d-7 !CO2
    z0(5) = 6.68d-11 !H+
    z0(6) = 1.79d-5 !HCO3-
    z0(7) = 2.79d-7 !CO32-
    z0(8) = 2.68d-7 !OH-
    z0(10) = 4.3d-6 !Mg2+
    z0(11) = 1.0d0 -z0(1)-z0(2)-z0(3)-z0(5)-z0(6)-z0(7)-z0(8)-z0(10) !SiO2

    !z0=0.0d0
     
    !z0(2)=4.6029673472613991E-004!1.0d0-z0(1)
    !z0(1)=1.0d0-z0(2)
    !chemi_mat(:,:) = 0
    
    !相安定解析の主要変数の初期値
    if (z0(1) >= z0(2)+z0(3)+z0(4)) then !液相から気相が出てくるパターン
        do i =1,com_2phase
            w0(i) = z0(i)*k0(i)
        end do
    else !気相から液相が出てくるパターン
        do i=1,com_2phase
            w0(i) = z0(i)/k0(i)
        end do
    end if
    V0 = 1.0d0 - z0(1)
    L0 = 1.0d0 -V0
    do i=1,com_2phase
        alpha0(i) = 2.0d0 *sqrt(w0(i))
    end do
    
    !write(*,*) lnk0
    
    
    do iteration =1,100
        if ((z0(1) >= z0(2)+z0(3)+z0(4))) then
            call phase_stability_liquid(alpha0,P0,z0,fxs)
        else
            call phase_stability_vapor(alpha0,P0,z0,fxs)
        end if
        
        call jacobian_mat(fxs,amat)
        call outxs(fxs,bmat)
        bmat = -bmat
        do i=1,com_2phase
            if (z0(i) == 0.0d0) then !存在しない成分は計算しないよ
                do j =1,com_2phase
                    amat(i,j) = 0.0d0
                    amat(j,i) = 0.0d0
                end do
                amat(i,i) = 1.0d0
                bmat(i) =0.0d0
            end if
        end do
            
        call pvgauss(com_2phase,amat,bmat)
        
        do i=1,com_2phase
            alpha0(i) =alpha0(i) +bmat(i)
            w0(i) = (alpha0(i) / 2.0d0) ** 2.0d0
        end do
        
        error = sqrt(dot_product(bmat,bmat))
        !write(*,*) error,iteration
        if (error < 10.0d0**(-5.0d0)) then
            exit
        end if
    end do
    
    
    wt =0.0d0
    do i=1,com_2phase
        wt =wt+w0(i)
    end do
    
    !write(*,*)(log(wt))
    lumda =1.0d0 -log(wt)
    if (lumda >= 1.0d0) then
        phase_judge0 = 1
        if (z0(1) > z0(2) +z0(3)+z0(4)) then
            V0 = 0.0d0
            L0 = 1.0d0 - V0
            phase0 = 1 !液相のみ
        else
            V0 = 1.0d0
            L0 = 1.0d0 - V0
            phase0 = 3 !気相のみ
            !write(*,*) V0
        end if
    else
        phase_judge0 = 2
        phase0 = 2 !2相
        v0 = z0(2) + z0(3) + z0(4)
        L0 = 1.0d0 - v0
        if (z0(1) >= (z0(2)+z0(3)+z0(4))) then !液相から気相が出てくるパターン
            do i=1,com_2phase
                if (z0(i) == 0d0) then
                    k0(i) = 1.0d0
                    lnk0(i) = log(k0(I))
                else
                    k0(i) = w0(i) / z0(i)
                    lnk0(i) = log(k0(i))
                end if
            end do
        else !気相から液相が出てくるパターン
            do i=1,com_2phase
                if (z0(i)== 0.0d0) then
                    k0(i) = 5.0
                    lnk0(i) = log(k0(i))
                else
                    k0(i) = z0(i) / w0(i)
                    lnk0(i) = log(k0(i))
                end if
            end do
        end if
    end if
    deallocate(amat,bmat)
    write(*,*) lnk0,lumda
    do i=1,com_2phase
    !    write(*,*) z0(i)
    end do
    !write(*,*) phase_judge0,'phase'
    
    
    allocate(cmat(com_2phase+1,com_2phase+1),dmat(com_2phase+1))
    if (phase_judge0 == 2) then
        do iteration = 1,100
            call ini_fla(P0,z0,lnk0,v0,z_factor0,fxs)
            call jacobian_mat(fxs,cmat)
            call outxs(fxs,dmat)
            dmat = -dmat
            do i=1,com_2phase
                if (z0(i) == 0.0d0) then !存在しない成分は計算しないよ
                    do j =1,com_2phase
                        cmat(i,j) = 0.0d0
                        cmat(j,i) = 0.0d0
                    end do
                    cmat(i,i) = 1.0d0
                    dmat(i) =0.0d0
                end if
            end do
            !do i=1,com_2phase+1
                !write(*,*) (cmat(i,j),j=1,com_2phase+1)
            !end do
            call pvgauss(com_2phase+1,cmat,dmat)
            do i=1,com_2phase
                lnk0(i) = lnk0(i) + dmat(i)
                k0(i) = exp(lnk0(i))
            end do
            v0 = v0 +dmat(com_2phase+1)
            
            error = sqrt(dot_product(dmat,dmat))
            if (error < 10.0d0**(-5.0d0)) then
                L0 = 1.0d0 - V0
                exit
            end if
        end do
    end if
    deallocate(cmat,dmat)
    !lnk1=lnk0(1) !chemicalの前に謎にlnkの値が0になるバグが10回に1回くらいおこるので、一旦保存して、あとで代入する。まじでなぞ
    !lnk2=lnk0(2)
    !lnk3=lnk0(3)
    !lnk4=lnk0(4)
    
    if (phase_judge0 == 2) then !2相
        call ini_N_2phase(z_factor0,V0,z0,x0,y0,k0,Nc0,Nt,Sw0)
    else if (phase0 == 1) then !液相のみ
        call ini_N_liquid(V0,z0,x0,y0,Nc0,Nt,Sw0)
    else !気相のみ
        call ini_N_vapor(V0,z0,x0,y0,Nc0,Nt,Sw0)
    end if
    
    !!岩石の初期モル数
    Nm0(1)=Nm1_ini
    Nm0(2)=Nm2_ini
    Nm0(3)=Nm3_ini
    Nm0(4)=Nm4_ini
    Nm0(5)=Nm5_ini
    do i=1,com_2phase+com_ion
        Nc0old(i) =Nc0(i)
    end do
    do i=1,com_mine
        Nm0old(i) =Nm0(i)
    end do
    P0old =p0
    
    
    !lnk0(1)=lnk1
    !lnk0(2)=lnk2
    !lnk0(3)=lnk3
    !lnk0(4)=lnk4
    
    !write(*,*) Nc0
    !write(*,*) Nc0(10)
    !液相がある場合、電離
    if (phase_judge0 == 2 .or. phase0 == 1) then
        allocate(emat(eq,eq),fmat(eq))
        do time =1,100!00!00 !time iteration
            !write(11,*) 'day',time
            
            do iteration =1,100
                !write(11,*) iteration
                call ini_chemi(P0,P0old,Nc0,Nc0old,Nm0,Nm0old,lnk0,V0,Sw0,fai0,chemi_mat,phase_judge0,z0,theta0,z_factor0,fxs)
                call jacobian_mat(fxs,emat)
                call outxs(fxs,fmat)
                fmat = -fmat
                !do i=1,eq
                    !write(10,*) (emat(i,j),j=1,eq)
                !end do
                do i=1,com_2phase
                    if (Nc0(i) == 0.0d0) then !存在しない成分は計算しないよ
            !            write(10,*) 'a'
                        do j =1,eq
                            emat(i,j) = 0.0d0
                            emat(j,i) = 0.0d0
                        end do
                        emat(i,i) = 1.0d0
                        fmat(i) =0.0d0
                    end if
                end do
                do i=1,com_2phase+com_ion
                    if (Nc0(i) == 0.0d0) then !存在しない成分は計算しないよ
            !            write(10,*) 'i'
                        do j =1,eq
                            emat(i+chemi,j) = 0.0d0
                            emat(j,i+chemi) = 0.0d0
                        end do
                        emat(i+chemi,i+chemi) = 1.0d0
                        fmat(i+chemi) =0.0d0
                    end if
                end do
                if (phase0 == 1) then !水相だけの時はrachford解かないよ
            !        write(10,*) 'u'
                    do i =1,eq
                        emat(i,eq) = 0.0d0
                        emat(eq,i) = 0.0d0
                    end do
                    emat(eq,eq) = 1.0d0
                    fmat(eq) = 0.0d0
                end if
                !do i=1,eq
                    !write(10,*) (emat(i,j),j=1,eq)
                !end do
            
                call pvgauss(eq,emat,fmat)
                
                
                do i=1,com_2phase
                    lnk0(i) = lnk0(i) + fmat(i)
                    k0(i) = exp(lnk0(i))
                end do
                do i=1,com_2phase+com_ion
                    Nc0(i) = Nc0(i) + fmat(i+com_2phase)
                end do
                do i=1,com_mine
                    Nm0(i) = Nm0(i) + fmat(i+com_2phase+com_2phase+com_ion)
                end do
                P0 = P0 + fmat(eq-1)
                v0 = v0 + fmat(eq)
                error = sqrt(dot_product(fmat,fmat))
                !write(10,*) time,iteration,error,'---------------------------------------------------'
                 !write(*,*) 'V:',V0
                 !write(*,*) p0
                !do i=1,eq
                    !write(10,*) (emat(i,j),j=1,eq)
                !end do
                if (error < 10.0d0**(-5.0d0)) then
                    goto 10
                end if
                
            end do
10          &
            do j=1,com_2phase+com_ion
                Nc0old(j) = Nc0(j)
            end do
            do j=1,com_mine
                Nm0old(j) = Nm0(j)
            end do
            P0 = inipressure
            P0old = P0
            
            
            if (iteration == 100) then
                write(*,*) '収束せず'
            end if

            
            
        end do
        deallocate(emat,fmat)
    end if
    
    
    !write(*,*) Nc0
    write(14,*) phase_judge0,'phases'
    write(14,*) P0,'pressure'
    write(14,*) 'z_factor:',z_factor0
    write(14,*) 'V0',V0,'Sw',Sw0
    write(14,*) 'Nc','z'
    do i=1,com_2phase+com_ion
        write(14,*) Nc0(i),z0(i)
    end do
    do i=1,com_mine
        write(14,*) Nm0(i)
    end do
    write(14,*) 'theta'
    do i=1,chemi+mine
        write(14,*) theta0(i)
    end do

    write(11,*) V0
    do j=1,com_2phase+com_ion
        write(11,*) Nc0(j),z0(j)
    end do 
    do j=1,com_mine
        write(11,*) Nm0(j)
    end do
    do j=1,com_2phase
        write(11,*) lnk0(j)
    end do
    write(11,*) fai0
    write(11,*) P0
    write(11,*) phase_judge0
    write(11,*) phase0

    

endprogram