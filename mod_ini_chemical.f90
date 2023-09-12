module mod_ini_chemical
    use mod_autodiff
    use mod_condition
    use mod_input
    use mod_fugacity
    
    implicit none
    
    integer,private::i,j
    
    contains
    
    subroutine ini_chemi(P0,P0old,Nc0,Nc0old,Nm0,Nm0old,lnk0,V0,Sw0,fai0,chemi_mat,phase_judge0,z0,theta0,z_factor0,g)
        implicit none
        integer,intent(in)::phase_judge0
        real(8),intent(inout)::P0,P0old,V0,Sw0,z_factor0,fai0
        real(8),intent(inout),dimension(com_2phase)::lnk0
        real(8),intent(inout),dimension(com_2phase+com_ion)::Nc0,Nc0old
        real(8),intent(out),allocatable,dimension(:)::z0
        real(8),intent(inout),dimension(com_mine)::Nm0,Nm0old
        real(8),intent(in),dimension(chemi+mine,com_all)::chemi_mat
        real(8),intent(out),allocatable,dimension(:)::theta0
        type(diffs),allocatable,intent(out)::g(:)
        type(diffs),allocatable,target::xd(:)
        type(diffs),pointer::P,Nc(:),lnk(:),V,Nm(:)
        real(8),allocatable,dimension(:)::x0
        
        real(8)::k0(com_2phase),v_L2,v_L3,v_L4,Sw00,Sg00,faimineold,faiold
        real(8),allocatable,dimension(:)::kakuninn,w0
        
        
        real(8),dimension(com_2phase)::Tc,Pc,acentric,a,b
        real(8),dimension(com_2phase,com_2phase)::bic
        type(diffs)::a_mix_V,b_mix_V,A_V,B_V,a_coef(3),C_V,D_V(com_2phase),E_V(com_2phase),G_V,sigma
        
        real(8)::A_ant,B_ant,C_ant,P1_sat,fai1_sat_L,A_rc,B_rc,C_rc,Tr1,kaku
        real(8),dimension(com_2phase)::A_hen,B_hen,C_hen,lnH_sat
        
        
        
        type(diffs)::Nt,z(com_2phase+com_ion),fra,k(com_2phase),x(com_2phase+com_ion),y(com_2phase),lnfai_V(com_2phase)&
                    ,lnfai_L(com_2phase),z_factor
        type(diffs)::MD_V,MV_V,v_L1,v_L(com_2phase+com_ion),MV_L,MD_L,Sw,Sg,faimine,fai,w(com_2phase+com_ion),Ac(chemi)&
                    ,Q(chemi+mine),theta(chemi+mine),min
        type(diffs)::rs(chemi+mine),rs_sum(com_all)
        
        
        real(8)::ks(chemi+mine),ke(chemi),km(mine),min0
        
        allocate(x0(eq),theta0(chemi+mine))
        
        
        !!自動微分の下準備
        do i=1,com_2phase
            k0(i) = lnk0(i)
            x0(i) = lnk0(i)
        end do
        do i=1,com_2phase+com_ion
            x0(com_2phase+i) = Nc0(i)
        end do
        do i=1,com_mine
            x0(com_2phase+com_2phase+com_ion+i) = Nm0(i)
        end do
        x0(eq-1) = P0
        x0(eq) = V0
        
        call diffsset1(x0,xd)
        call sizeset(x0,g)
        
        lnk => xd(1:com_2phase)
        Nc => xd(com_2phase+1:com_2phase+com_2phase+com_ion)
        Nm => xd(com_2phase+com_2phase+com_ion+1:com_2phase+com_2phase+com_ion+com_mine)
        P => xd(eq-1)
        V => xd(eq)
        
        
        !!モル分率
        call residualvectorset3(eq,Nt)
        do i=1,com_2phase +com_ion
            Nt = Nt +Nc(i)
        end do
        do i=1,com_2phase+com_ion
            z(i) = Nc(i)/Nt
        end do
        call outxs(z,z0)
        !write(*,*) lnk0
        !
        !!!rachford-rice
        if (phase_judge0 == 2) then
            call residualvectorset3(eq,fra)
            do i=1,com_2phase
                fra=fra+(1.0d0-exp(lnk(i)))*z(i)/(1.0d0-V+V*exp(lnk(i)))
                !call out_diffsx(fra,kaku)
                !write(*,*) kaku
            end do
            do i=com_2phase+1,com_2phase+com_ion
                fra = fra+z(i)/(1.0d0-V)   
            end do
        else
            call residualvectorset3(eq,fra)
            call residualvectorset3(eq,V)
        end if
        
        
        
        !モル分率
        if (phase_judge0 == 2) then
            do i=1,com_2phase
                if(z0(i) == 0.0d0) then
                    call residualvectorset3(eq,k(i))
                else
                    k(i) = exp(lnk(i))
                end if
                x(i) = z(i) /(1.0d0-V+V*k(i))
                y(i) = x(i) * k(i)
            end do
            do i=com_2phase+1,com_2phase+com_ion
                x(i) = z(i)/(1.0d0-V)
            end do
        else
            do i=1,com_2phase
                call residualvectorset4(k0(i),eq,k(i))
                call residualvectorset3(eq,y(i))
            end do
            x = z
        end if
        
        Tc(1) = 647.30d0
        Tc(2) = 304.2d0
        Tc(3) = 190.6d0
        Tc(4) = 373.2d0
        Pc(1) = 217.6d0 * 101325.d0
        Pc(2) = 72.8d0 * 101325.d0
        Pc(3) = 45.4d0 * 101325.d0
        Pc(4) = 88.2d0 * 101325.d0
        acentric(1) = 0.344d0
        acentric(2) = 0.225d0
        acentric(3) = 0.008d0
        acentric(4) = 0.1d0
        bic(1,1)=0.0d0
        bic(1,2)=0.2d0
        bic(1,3)=0.4907d0
        bic(1,4)=0.12d0
        bic(2,1)=bic(1,2)
        bic(2,2)=0.0d0
        bic(2,3)=0.105d0
        bic(2,4)=0.135d0
        bic(3,1)=bic(1,3)
        bic(3,2)=bic(2,3)
        bic(3,3)=0.0d0
        bic(3,4)=0.07d0
        bic(4,1)=bic(1,4)
        bic(4,2)=bic(2,4)
        bic(4,3)=bic(3,4)
        bic(4,4)=0.0d0
        
        !call outxs(y,kakuninn)
        !write(*,*) kakuninn
        if (phase_judge0 == 2) then
            !call vapor_fugacity_chemi(y,lnfai_V,P,z_factor0,z_factor)
            
        
            call residualvectorset3(eq,a_mix_V)
            !call residualvectorset4(0.0d0,eq,a_mix_V)
            call residualvectorset3(eq,b_mix_V)
            do j=1,com_2phase
                a(j) = (0.457235528921d0*(R**2.0d0)*(tc(j)**2.0d0)/pc(j)) * ((1.0d0+(0.37464d0+1.54226d0*acentric(j)-0.26992*&
                (acentric(j)**2.0d0))*(1.0d0-sqrt(temp/tc(j))) ))**2.0d0
                b(j) = 0.0777960739039d0*R*tc(j)/pc(j)
                b_mix_V=b_mix_V+y(j)*b(j)
            end do
            !write(*,*) a(4)
        
        
            do j=1,com_2phase
                do i=1,com_2phase
                    a_mix_V=a_mix_V+y(j)*y(i)*(1.0d0-bic(j,i))*sqrt(a(j)*a(i))
                end do
            end do
            !call outxs(,kakuninn)
            !write(*,*) kakuninn
        
            A_V          = (p/((R*temp)**2.0d0)) * a_mix_V
            B_V          = (p/(R*temp))        * b_mix_V
            a_coef(1)       =(B_V-1.0d0)
            a_coef(2)       = A_V - 2.0d0*B_V - 3.0d0*B_V**2.0d0
            a_coef(3)       =(B_V**2.0d0 + B_V**3.0d0-A_V*B_V)
        
            call cubic_eq_d(a_coef(1),a_coef(2),a_coef(3),z_factor)
            call out_diffsx(z_factor,z_factor0)
            !?---------------------------------------------------------
            !write(*,*) z_factor0
        
        
            !気相のフガシティ係数
            C_V = z_factor-B_V
            
            do j=1,com_2phase
                call residualvectorset3(eq,sigma)
                D_V(j)=b(j)*(z_factor-1.0d0)/b_mix_V
                do i=1,com_2phase
                    sigma=sigma+y(i)*(1.0d0-bic(j,i))*dsqrt(a(j)*a(i))
                end do
                E_V(j)=(A_V/(2.0d0*dsqrt(2.0d0)*B_V))*((2.0d0/a_mix_V)*sigma-b(j)/b_mix_V)
            end do
            G_V=(z_factor+B_V*(1.0d0+dsqrt(2.0d0)))/(z_factor+B_V*(1.0d0-dsqrt(2.0d0)))
            do j=1,com_2phase
                lnfai_V(j)=D_V(j)-log(C_V)-E_V(j)*log(G_V)
            end do
            !lnfai_V = y
            !call residualvectorset3(eq,lnfai_V)
            !call outxs(lnfai_V,kakuninn)
            !write(*,*) kakuninn(1),kakuninn(4)
            !write(*,*) z_factor0
        else
            do i=1,com_2phase
                call residualvectorset3(eq,lnfai_V(i))
            end do
        end if
        !
        !!call outxs(lnfai_V,kakuninn)
        !!write(*,*) kakuninn
        !
        !!call liquid_fugacity_chemi(lnfai_L,P,v_L1,v_L2,v_L3,v_L4)
        !H2O
        A_ant=8.02754d0
        B_ant=1705.616d0
        C_ant=-41.745d0
        P1_sat=133.322368d0*10.0d0**(A_ant-B_ant/(Temp+C_ant))
        
        if (Temp_F <= 90) then
            fai1_sat_L=1.0d0
        else
            fai1_sat_L=0.9958+9.68330d0*Temp_F/(10.0d0**5.d0)-6.17*Temp_F**2.0d0/(10.0d0**7.d0)-&
            3.08333*Temp_F**3.0d0/(10.0d0**10.0d0)
        end if
        
        A_rc=5.916365d0-0.01035794d0*temp+0.9270048d0*(10.0d0**(-5.0d0))*(temp**2.0d0)-1127.522d0/temp+100674.1d0/(temp**2.0d0)
        B_rc=0.5204914d0*(10.0d0**(-2.0d0))-0.10482101d0*(10.0d0**(-4.0d0))*temp+0.8328532d0*(10.0d0**(-8.0d0))*(temp**2.0d0)&
            -1.1702939d0/temp+102.2783d0/(temp**2.0d0)
        C_rc=0.118547d0*(10.0d0**(-7.0d0))-0.6599143d0*(10.0d0**(-10.0d0))*temp
        
        !液相のモル体積
        v_L1=(M1/10.0d0**3.0d0)*(A_rc-B_rc*(P/(9.80665d0*10.0d0**4.0d0))-C_rc*(P/(9.80665d0*10.0d0**4.0d0))**2.0d0)
        v_L2=(-47.7518+4.336154/10.0d0*temp-5.945771/10.0d0**4.0d0*temp**2.0d0)/10.0d0**6.0d0
        v_L3=exp(3.541d0+1.23d-3*(temp-273.15))*10**(-6.d0)
        v_L4=(160.5567d0-5.538776d0*10.0d0**(-1.0d0)*temp)*10**(-6.d0)
        
        
        lnfai_L(1)=log((P1_sat*fai1_sat_L*exp(M1/(10.0d0**3.0d0*R*Temp)*(A_rc*(p-P1_sat)-B_rc*(P**2.0d0-P1_sat**2.0d0)/&
                    (2.0d0*9.80665d0*10.0d0**4.d0)&
                    -C_rc*(P**3.0d0-P1_sat**3.0d0)/(3.0d0*(9.80665d0*10.0d0**4.0d0)**2.0d0))))/P)
        !______________________________________________________________________________________________________
        
        
        !!CO2の液相のフガシティ(Henry)----------------------------------------
        A_hen(2)=-9.4234d0
        B_hen(2)=4.0087d0
        C_hen(2)=10.3199d0
        !ここはH2Oなので注意!
        Tr1=Temp/Tc(1)
        
        lnH_sat(2)=log(P1_sat)+A_hen(2)/Tr1+B_hen(2)*((1.0d0-Tr1)**0.355d0)/Tr1+C_hen(2)*dexp(1.0d0-Tr1)/Tr1**0.41d0
        lnfai_L(2)=lnH_sat(2)+v_L2*(P-P1_sat)/(R*temp)-log(P)
        !?------------------------------------------------------------------
        !!CH4の液相のフガシティ(Henry)----------------------------------------
        A_hen(3)=-11.0094d0
        B_hen(3)=4.8362d0
        C_hen(3)=12.5220d0
        !ここはH2Oなので注意!
        Tr1=Temp/Tc(1)
        
        lnH_sat(3)=log(P1_sat)+A_hen(3)/Tr1+B_hen(3)*((1.0d0-Tr1)**0.355d0)/Tr1+C_hen(3)*dexp(1.0d0-Tr1)/Tr1**0.41d0
        lnfai_L(3)=lnH_sat(3)+v_L3*(P-P1_sat)/(R*temp)-log(P)
        !?------------------------------------------------------------------
        !!H2Sの液相のフガシティ(Henry)----------------------------------------
        A_hen(4)=-5.7131d0
        B_hen(4)=5.3727d0
        C_hen(4)=5.4227d0
        !ここはH2Oなので注意!
        Tr1=Temp/Tc(1)
        
        lnH_sat(4)=log(P1_sat)+A_hen(4)/Tr1+B_hen(4)*((1.0d0-Tr1)**0.355d0)/Tr1+C_hen(4)*dexp(1.0d0-Tr1)/Tr1**0.41d0
        lnfai_L(4)=lnH_sat(4)+v_L4*(P-P1_sat)/(R*temp)-log(P)
        !?------------------------------------------------------------------
        !
        !
        !!気相について
        if (phase_judge0 == 2) then
            MD_V = P/z_factor/R/temp
            MV_V = 1.0d0/MD_V
        else
            call residualvectorset3(eq,MD_V)
            call residualvectorset3(eq,MV_V)
        end if
        
        !!液相について
        v_L(1) = v_L1![m^3/mol]
        !write(*,*) v_L2,v_L3,v_L4
        
        call residualvectorset4(v_L2,eq,v_L(2))
        call residualvectorset4(v_L3,eq,v_L(3))
        call residualvectorset4(v_L4,eq,v_L(4))
        do i=1,com_ion
            call residualvectorset3(eq,v_L(i+com_2phase))
            !call residualvectorset4(v_L2/10,eq,v_L(i+com_2phase))
            !v_L(i+com_2phase) = v_L(1)/10.0d0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!要検討!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        end do
        !
        call residualvectorset3(eq,MV_L)
        do i=1,com_2phase+com_ion
            MV_L = MV_L + x(i)*v_L(i)
        end do
        MD_L =1.0d0/MV_L
        
        !!飽和率
        Sw = Nt * (1.0d0-V) * MV_L
        Sg = Nt * V * MV_V
        call out_diffsx(Sw,Sw00)
        call out_diffsx(Sg,Sg00)
        !write(*,*) 'Sw,Sg',Sw00,Sg00
        
        !!孔隙率
        faimine = faiini - (Nm(1)-Nm1_ini)/Nm1_MD - (Nm(2)-Nm2_ini)/Nm2_MD - (Nm(3)-Nm3_ini)/Nm3_MD&
                         - (Nm(4)-Nm4_ini)/Nm4_MD - (Nm(5)-Nm5_ini)/Nm5_MD
        faimineold = faiini - (Nm0old(1)-Nm1_ini)/Nm1_MD - (Nm0old(2)-Nm2_ini)/Nm2_MD - (Nm0old(3)-Nm3_ini)/Nm3_MD&
                            - (Nm0old(4)-Nm4_ini)/Nm4_MD - (Nm0old(5)-Nm5_ini)/Nm5_MD
        
        fai = faimine * (1.0d0+Cr*(P-iniPressure))
        faiold = faimineold*(1.0d0+Cr*(P0old-iniPressure))
        
        call out_diffsx(fai,fai0)
        !write(*,*) fai0
        
        !!!!!!chemical!!!!!
        
        
        !!速度定数[mol/m^3/s]------------ここでは初期条件で平衡にするためなので適当
        do i=1,chemi+mine
            ks(i) =1.0d0*10.0d0**(-3.0d0)
        end do
        
        
        !!平衡定数
        ke(1) = ke1
        ke(2) = ke2
        ke(3) = ke3
        ke(4) = ke4
        ke(5) = ke5
        km(1) = km1
        km(2) = km2
        km(3) = km3
        km(4) = km4
        km(5) = km5
        
        
        !!活量
        do i=1,com_2phase+com_ion
            !w(i)=MD_L*x(i)/(MD_L*x(1)*18.015d0*10.0d0**(-3.0d0))
            w(i)=MD_L*x(i)/(MD_L*x(1)*18.015d0*10.0d0**(-3.0d0)) 
        end do
        call residualvectorset4(1.0d0,eq,w(1))
        call outxs(w,w0)
        !write(12,*) w0
        
        !!比表面積
        Ac(1)=AA1*Nm(1)/Nm1_ini
        Ac(2)=AA2*Nm(2)/Nm2_ini
        Ac(3)=AA3*Nm(3)/Nm3_ini
        Ac(4)=AA4*Nm(4)/Nm4_ini
        Ac(5)=AA5*Nm(5)/Nm5_ini
        
        !!活量積
        do i=1,chemi
            call residualvectorset4(1.0d0,eq,Q(i))
            do j=1,com_2phase+com_ion
                if (Nc0(j) /= 0 ) then
                    Q(i)=Q(i)*w(j)**chemi_mat(i,j)
                end if
            end do
            call out_diffsx(Q(i),kaku)
            !write(*,*) i,kaku
            theta(i)=Q(i)/ke(i)
        end do
        do i=1,mine
            call residualvectorset4(1.0d0,eq,Q(i+chemi))
            do j=1,com_2phase+com_ion
                if (Nc0(j) /= 0 ) then
                    Q(i+chemi)=Q(i+chemi)*w(j)**chemi_mat(i+chemi,j)
                end if
            end do
            call out_diffsx(Q(i+chemi),kaku)
            !write(*,*) i+chemi,kaku
            theta(i+chemi)=Q(i+chemi)/km(i)
        end do
        !call outxs(w,theta0)
        !write(11,*) theta0
        !call outxs(Q,theta0)
        !write(11,*) theta0
        call outxs(theta,theta0)
        !write(11,*) theta0
        !    
        !!化学反応の反応速度
        
        do i=1,chemi
            min0=100000000000.0d0
            call residualvectorset4(min0,eq,min)
            if (theta0(i) < 1.0d0) then !順反応
                do j=2,com_2phase+com_ion
                    if (chemi_mat(i,j) < 0.0d0) then !減少成分
                        if (min0 > w0(j)) then !最小濃度成分
                            min0=w0(j)
                            min=w(j)
                        end if
                    end if
                end do
            else !逆反応
                do j=2,com_2phase+com_ion
                    if (chemi_mat(i,j) > 0.0d0) then !減少成分
                        if (min0 > w0(j)) then !最小濃度成分
                            min0=w0(j)
                            min=w(j)
                        end if
                    end if
                end do
            end if
            !if (abs(theta0(i)-1.0d0) < 10.0d0) then
            !    call residualvectorset3(eq,rs(i))
            !else
                rs(i)=ks(i)*(1.0d0-theta(i))*fai*Sw*min
            !end if
        end do
        
        !!鉱物反応
        do i=1,mine
            rs(i+chemi)=ks(i+chemi)*(1.0d0-theta(i+chemi))*fai*Sw
        end do
        
        !!反応速度の合計
        do i=1,com_all
            call residualvectorset3(eq,rs_sum(i))
            do j=1,chemi+mine
                rs_sum(i)=rs_sum(i)+rs(j)*chemi_mat(j,i)
            end do
        end do
        call outxs(Q,kakuninn)
        !write(*,*) kakuninn
        
        
        
        !!熱力学平衡条件式
        do i=1,com_2phase
            g(i) = lnk(i) + lnfai_V(i) - lnfai_L(i)
        end do
        
        !!物質収支式
        do i=1,com_2phase+com_ion
            g(i+com_2phase) = (-1.0d0)*(Nc(i)*fai-Nc0old(i)*faiold)/dt + rs_sum(i)
        end do
        do i=1,com_mine
            g(i+com_2phase+com_2phase+com_ion) = (-1.0d0)*(Nm(i)*(1.0d0-fai)-Nm0old(i)*(1.0d0-faiold))/dt + rs_sum(i)
        end do
        
        g(eq-1) = Sw+Sg-1.0d0
        g(eq) = fra

        
        
    end subroutine
    
    end module
    
    
    