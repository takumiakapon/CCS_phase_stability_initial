module mod_fugacity
    use mod_autodiff
    use mod_condition
    use mod_input
    use mod_cubic_eq_d
    
    implicit none
    
    integer,private::i,j,ii
    
    contains
    
    subroutine vapor_fugacity(y,lnfai_V,P)
        implicit none
        type(diffs),intent(in)::y(com_2phase)
        real(8),intent(in)::P
        type(diffs),intent(out)::lnfai_V(com_2phase)
        
        real(8),dimension(com_2phase)::Tc,Pc,acentric,a,b
        real(8),dimension(com_2phase,com_2phase)::bic
        type(diffs)::a_mix_V,b_mix_V,A_V,B_V,z_factor,C_V,sigma,G_V
        type(diffs)::a_coef(3)
        type(diffs),dimension(com_2phase)::D_V,E_V
        real(8)::z_factor0
        real(8),allocatable,dimension(:)::kakuninn
        real(8)::kaku
        
        
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
        
        call residualvectorset3(com_2phase,a_mix_V)
        call residualvectorset3(com_2phase,b_mix_V)
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
        
        
        !!気相のフガシティ係数
        C_V = z_factor-B_V
        !write(*,*) 'e'
        do j=1,com_2phase
        !    write(*,*) 'i'
            call residualvectorset3(com_2phase,sigma)
        !    write(*,*) 'u'
            D_V(j)=b(j)*(z_factor-1.0d0)/b_mix_V
            !call outxs(D_V,kakuninn)
            
        !    write(*,*) 'k'
            do i=1,com_2phase
                sigma=sigma+y(i)*(1.0d0-bic(j,i))*sqrt(a(j)*a(i))
                !call out_diffsx(sigma,kaku)
                !write(*,*) kaku
            end do
        !    write(*,*) 's'
            E_V(j)=(A_V/(2.0d0*sqrt(2.0d0)*B_V))*((2.0d0/a_mix_V)*sigma-b(j)/b_mix_V)!?ここ？2成分目が計算できていない？
            !#TODOよくわからないからサブルーチンじゃなくて普通に書く？
        !    write(*,*) 'n'
            call outxs(E_V,kakuninn)
        !    write(*,*) kakuninn
        end do
        
        G_V=(z_factor+B_V*(1.0d0+sqrt(2.0d0)))/(z_factor+B_V*(1.0d0-sqrt(2.0d0)))
        do j=1,com_2phase
            lnfai_V(j)=D_V(j)-log(C_V)-E_V(j)*log(G_V)
        end do
        
        call outxs(lnfai_V,kakuninn)
        !write(*,*) kakuninn(1),kakuninn(4)
        !write(*,*) z_factor0
        !write(*,*) 'o'
        
    end subroutine
    
    
    
    subroutine liquid_fugacity(lnfai_L,P)
        implicit none
        real(8),intent(in)::P
        type(diffs),intent(out)::lnfai_L(com_2phase)
        
        real(8),dimension(com_2phase)::Tc,Pc,acentric,A_hen,B_hen,C_hen,lnH_sat,lnfai_L0
        real(8),dimension(com_2phase,com_2phase)::bic
        
        real(8)::A_ant,B_ant,C_ant,P1_sat,fai1_sat_L,A_rc,B_rc,C_rc,v_L1,v_L2,v_L3,v_L4,Tr1
        
        
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
        !write(*,*) v_L4

        lnfai_L0(1)=log((P1_sat*fai1_sat_L*exp(M1/(10.0d0**3.0d0*R*Temp)*(A_rc*(p-P1_sat)-B_rc*(P**2.0d0-P1_sat**2.0d0)/&
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
        lnfai_L0(2)=lnH_sat(2)+v_L2*(P-P1_sat)/(R*temp)-log(P)
        !?------------------------------------------------------------------
        !!CH4の液相のフガシティ(Henry)----------------------------------------
        A_hen(3)=-11.0094d0
        B_hen(3)=4.8362d0
        C_hen(3)=12.5220d0
        !ここはH2Oなので注意!
        Tr1=Temp/Tc(1)

        lnH_sat(3)=log(P1_sat)+A_hen(3)/Tr1+B_hen(3)*((1.0d0-Tr1)**0.355d0)/Tr1+C_hen(3)*dexp(1.0d0-Tr1)/Tr1**0.41d0
        lnfai_L0(3)=lnH_sat(3)+v_L3*(P-P1_sat)/(R*temp)-log(P)
        !?------------------------------------------------------------------
        !!H2Sの液相のフガシティ(Henry)----------------------------------------
        A_hen(4)=-5.7131d0
        B_hen(4)=5.3727d0
        C_hen(4)=5.4227d0
        !ここはH2Oなので注意!
        Tr1=Temp/Tc(1)

        lnH_sat(4)=log(P1_sat)+A_hen(4)/Tr1+B_hen(4)*((1.0d0-Tr1)**0.355d0)/Tr1+C_hen(4)*dexp(1.0d0-Tr1)/Tr1**0.41d0
        lnfai_L0(4)=lnH_sat(4)+v_L4*(P-P1_sat)/(R*temp)-log(P)
        !?------------------------------------------------------------------
        
        do i=1,com_2phase
            call residualvectorset4(lnfai_L0(i),com_2phase,lnfai_L(i))
        end do
        
        !write(*,*) lnfai_L0(1),lnfai_L0(4
        
        
    end subroutine
    
    
    
    subroutine vapor_fugacity_ini(y,lnfai_V,P,z_factor0)
        implicit none
        type(diffs),intent(in)::y(com_2phase)
        real(8),intent(inout)::P,z_factor0
        type(diffs),intent(out)::lnfai_V(com_2phase)
        
        real(8),dimension(com_2phase)::Tc,Pc,acentric,a,b
        real(8),dimension(com_2phase,com_2phase)::bic
        type(diffs)::a_mix_V,b_mix_V,A_V,B_V,z_factor,C_V,sigma,G_V
        type(diffs)::a_coef(3)
        type(diffs),dimension(com_2phase)::D_V,E_V
        real(8),allocatable,dimension(:)::kakuninn
        real(8)::kaku
        
        
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
        
        call residualvectorset3(com_2phase+1,a_mix_V)
        call residualvectorset3(com_2phase+1,b_mix_V)
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
        
        
        !!気相のフガシティ係数
        C_V = z_factor-B_V
        
        do j=1,com_2phase
            call residualvectorset3(com_2phase+1,sigma)
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
        call outxs(lnfai_V,kakuninn)
        !write(*,*) kakuninn(1),kakuninn(4)
        !write(*,*) z_factor0
        
    end subroutine
    
    subroutine liquid_fugacity_ini(lnfai_L,P)
        implicit none
        real(8),intent(in)::P
        type(diffs),intent(out)::lnfai_L(com_2phase)
        
        real(8),dimension(com_2phase)::Tc,Pc,acentric,A_hen,B_hen,C_hen,lnH_sat,lnfai_L0
        real(8),dimension(com_2phase,com_2phase)::bic
        
        real(8)::A_ant,B_ant,C_ant,P1_sat,fai1_sat_L,A_rc,B_rc,C_rc,v_L1,v_L2,v_L3,v_L4,Tr1
        
        
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
        

        lnfai_L0(1)=log((P1_sat*fai1_sat_L*exp(M1/(10.0d0**3.0d0*R*Temp)*(A_rc*(p-P1_sat)-B_rc*(P**2.0d0-P1_sat**2.0d0)/&
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
        lnfai_L0(2)=lnH_sat(2)+v_L2*(P-P1_sat)/(R*temp)-log(P)
        !?------------------------------------------------------------------
        !!CH4の液相のフガシティ(Henry)----------------------------------------
        A_hen(3)=-11.0094d0
        B_hen(3)=4.8362d0
        C_hen(3)=12.5220d0
        !ここはH2Oなので注意!
        Tr1=Temp/Tc(1)

        lnH_sat(3)=log(P1_sat)+A_hen(3)/Tr1+B_hen(3)*((1.0d0-Tr1)**0.355d0)/Tr1+C_hen(3)*dexp(1.0d0-Tr1)/Tr1**0.41d0
        lnfai_L0(3)=lnH_sat(3)+v_L3*(P-P1_sat)/(R*temp)-log(P)
        !?------------------------------------------------------------------
        !!H2Sの液相のフガシティ(Henry)----------------------------------------
        A_hen(4)=-5.7131d0
        B_hen(4)=5.3727d0
        C_hen(4)=5.4227d0
        !ここはH2Oなので注意!
        Tr1=Temp/Tc(1)

        lnH_sat(4)=log(P1_sat)+A_hen(4)/Tr1+B_hen(4)*((1.0d0-Tr1)**0.355d0)/Tr1+C_hen(4)*dexp(1.0d0-Tr1)/Tr1**0.41d0
        lnfai_L0(4)=lnH_sat(4)+v_L4*(P-P1_sat)/(R*temp)-log(P)
        !?------------------------------------------------------------------
        
        do i=1,com_2phase
            call residualvectorset4(lnfai_L0(i),com_2phase+1,lnfai_L(i))
        end do
        
        !write(*,*) lnfai_L0(1),lnfai_L0(4)
        
        
        
    end subroutine
    
    
    subroutine vapor_fugacity_chemi(y,lnfai_V,P,z_factor0,z_factor)
        implicit none
        type(diffs),intent(inout)::y(com_2phase),P,z_factor
        real(8),intent(inout)::z_factor0
        type(diffs),intent(out)::lnfai_V(com_2phase)
        
        real(8),dimension(com_2phase)::Tc,Pc,acentric,a,b
        real(8),dimension(com_2phase,com_2phase)::bic
        type(diffs)::a_mix_V,b_mix_V,A_V,B_V,C_V,sigma,G_V
        type(diffs)::a_coef(3)
        type(diffs),dimension(com_2phase)::D_V,E_V
        real(8),allocatable,dimension(:)::kakuninn
        real(8)::kaku
        
        
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
        
        call residualvectorset3(eq,a_mix_V)
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
        call outxs(lnfai_V,kakuninn)
        !write(*,*) kakuninn(1),kakuninn(4)
        !write(*,*) z_factor0
        
    end subroutine
    
    
    subroutine liquid_fugacity_chemi(lnfai_L,P,v_L1,v_L2,v_L3,v_L4)
        implicit none
        type(diffs),intent(in)::P
        type(diffs),intent(out)::lnfai_L(com_2phase),v_L1
        real(8),intent(out)::v_L2,v_L3,v_L4
        
        real(8),dimension(com_2phase)::Tc,Pc,acentric,A_hen,B_hen,C_hen,lnH_sat,lnfai_L0
        real(8),dimension(com_2phase,com_2phase)::bic
        
        real(8)::A_ant,B_ant,C_ant,P1_sat,fai1_sat_L,A_rc,B_rc,C_rc,Tr1
               
        
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
        
        !do i=1,com_2phase
            !call residualvectorset4(lnfai_L0(i),eq,lnfai_L(i))
        !end do
        
        !write(*,*) lnfai_L0(1),lnfai_L0(4)
        
        
        
    end subroutine

    subroutine vapor_fugacity_main(y,lnfai_V,P,z_factor0,z_factor,q_judge,phase_judge)
        implicit none
        type(diffs),intent(inout)::y(com_2phase,n),P(n),z_factor(n)
        real(8),intent(inout)::z_factor0(n)
        integer,intent(inout)::q_judge,phase_judge(n)
        type(diffs),intent(out)::lnfai_V(com_2phase,n)
        
        real(8),dimension(com_2phase)::Tc,Pc,acentric,a,b
        real(8),dimension(com_2phase,com_2phase)::bic
        type(diffs)::a_mix_V(n),b_mix_V(n),A_V(n),B_V(n),C_V(n),sigma(com_2phase,n),G_V(n)
        type(diffs)::a_coef(3)
        type(diffs),dimension(com_2phase,n)::D_V,E_V
        real(8),allocatable,dimension(:)::kakuninn
        real(8)::kaku
        
        
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
        
        do i=1,n
            if (phase_judge(i) == 2) then
                call residualvectorset3(n*eq+q_judge,a_mix_V(i))
                call residualvectorset3(n*eq+q_judge,b_mix_V(i))
                do j=1,com_2phase
                    a(j) = (0.457235528921d0*(R**2.0d0)*(tc(j)**2.0d0)/pc(j)) * ((1.0d0+(0.37464d0+1.54226d0*acentric(j)-0.26992*&
                    (acentric(j)**2.0d0))*(1.0d0-sqrt(temp/tc(j))) ))**2.0d0
                    b(j) = 0.0777960739039d0*R*tc(j)/pc(j)
                    b_mix_V(i)=b_mix_V(i)+y(j,i)*b(j)
                end do
                !write(*,*) a(4)
                do j=1,com_2phase
                    do ii=1,com_2phase
                        a_mix_V(i)=a_mix_V(i)+y(j,i)*y(ii,i)*(1.0d0-bic(j,ii))*sqrt(a(j)*a(ii))
                    end do
                end do
                !call outxs(,kakuninn)
                !write(*,*) kakuninn
        
                A_V(i)          = (p(i)/((R*temp)**2.0d0)) * a_mix_V(i)
                B_V(i)          = (p(i)/(R*temp))        * b_mix_V(i)
                a_coef(1)       =(B_V(i)-1.0d0)
                a_coef(2)       = A_V(i) - 2.0d0*B_V(i) - 3.0d0*B_V(i)**2.0d0
                a_coef(3)       =(B_V(i)**2.0d0 + B_V(i)**3.0d0-A_V(i)*B_V(i))
        
                call cubic_eq_d(a_coef(1),a_coef(2),a_coef(3),z_factor(i))
                call out_diffsx(z_factor(i),z_factor0(i))
            end if
        end do
        !?---------------------------------------------------------
        !write(*,*) z_factor0
        
        
        !気相のフガシティ係数
        do i=1,n
            if (phase_judge(i) == 2) then
                C_V(i) = z_factor(i)-B_V(i)
        
                
                do j=1,com_2phase
                    call residualvectorset3(n*eq+q_judge,sigma(j,i))
                    D_V(j,i)=b(j)*(z_factor(i)-1.0d0)/b_mix_V(i)
                    do ii=1,com_2phase
                        sigma(j,i)=sigma(j,i)+y(ii,i)*(1.0d0-bic(j,ii))*dsqrt(a(j)*a(ii))
                    end do
                    E_V(j,i)=(A_V(i)/(2.0d0*dsqrt(2.0d0)*B_V(i)))*((2.0d0/a_mix_V(i))*sigma(j,i)-b(j)/b_mix_V(i))
                end do
                G_V(i)=(z_factor(i)+B_V(i)*(1.0d0+dsqrt(2.0d0)))/(z_factor(i)+B_V(i)*(1.0d0-dsqrt(2.0d0)))
                do j=1,com_2phase
                    lnfai_V(j,i)=D_V(j,i)-log(C_V(i))-E_V(j,i)*log(G_V(i))
                end do
            else
                do j=1,com_2phase
                    call residualvectorset3(n*eq+q_judge,lnfai_V(j,i))
                end do
            end if
        end do
        !lnfai_V = y
        !call residualvectorset3(eq,lnfai_V)
        !call outxs(lnfai_V,kakuninn)
        !write(*,*) kakuninn(1),kakuninn(4)
        !write(*,*) z_factor0
        
    end subroutine
    
    
    subroutine liquid_fugacity_main(lnfai_L,P,v_L1,v_L2,v_L3,v_L4)
        implicit none
        type(diffs),intent(in)::P(n)
        type(diffs),intent(out)::lnfai_L(com_2phase,n),v_L1(n)
        real(8),intent(out)::v_L2(n),v_L3(n),v_L4(n)
        
        real(8),dimension(com_2phase)::Tc,Pc,acentric,A_hen,B_hen,C_hen,lnH_sat,lnfai_L0
        real(8),dimension(com_2phase,com_2phase)::bic
        
        real(8)::A_ant,B_ant,C_ant,P1_sat,fai1_sat_L,A_rc,B_rc,C_rc,Tr1
               
        
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
        do i=1,n
            v_L1(i)=(M1/10.0d0**3.0d0)*(A_rc-B_rc*(P(i)/(9.80665d0*10.0d0**4.0d0))&
                    -C_rc*(P(i)/(9.80665d0*10.0d0**4.0d0))**2.0d0)
            v_L2(i)=(-47.7518+4.336154/10.0d0*temp-5.945771/10.0d0**4.0d0*temp**2.0d0)/10.0d0**6.0d0
            v_L3(i)=exp(3.541d0+1.23d-3*(temp-273.15))*10**(-6.d0)
            v_L4(i)=(160.5567d0-5.538776d0*10.0d0**(-1.0d0)*temp)*10**(-6.d0)

            lnfai_L(1,i)=log((P1_sat*fai1_sat_L*exp(M1/(10.0d0**3.0d0*R*Temp)*(A_rc*(p(i)-P1_sat)&
                    -B_rc*(P(i)**2.0d0-P1_sat**2.0d0)/(2.0d0*9.80665d0*10.0d0**4.d0)&
                    -C_rc*(P(i)**3.0d0-P1_sat**3.0d0)/(3.0d0*(9.80665d0*10.0d0**4.0d0)**2.0d0))))/P(i))
        end do
        

        
        !______________________________________________________________________________________________________
        
        
        !!CO2の液相のフガシティ(Henry)----------------------------------------
        A_hen(2)=-9.4234d0
        B_hen(2)=4.0087d0
        C_hen(2)=10.3199d0
        !ここはH2Oなので注意!
        Tr1=Temp/Tc(1)

        lnH_sat(2)=log(P1_sat)+A_hen(2)/Tr1+B_hen(2)*((1.0d0-Tr1)**0.355d0)/Tr1+C_hen(2)*dexp(1.0d0-Tr1)/Tr1**0.41d0
        do i=1,n
            lnfai_L(2,i)=lnH_sat(2)+v_L2(i)*(P(i)-P1_sat)/(R*temp)-log(P(i))
        end do
        !?------------------------------------------------------------------
        !!CH4の液相のフガシティ(Henry)----------------------------------------
        A_hen(3)=-11.0094d0
        B_hen(3)=4.8362d0
        C_hen(3)=12.5220d0
        !ここはH2Oなので注意!
        Tr1=Temp/Tc(1)

        lnH_sat(3)=log(P1_sat)+A_hen(3)/Tr1+B_hen(3)*((1.0d0-Tr1)**0.355d0)/Tr1+C_hen(3)*dexp(1.0d0-Tr1)/Tr1**0.41d0
        do i=1,n
            lnfai_L(3,i)=lnH_sat(3)+v_L3(i)*(P(i)-P1_sat)/(R*temp)-log(P(i))
        end do
        !?------------------------------------------------------------------
        !!H2Sの液相のフガシティ(Henry)----------------------------------------
        A_hen(4)=-5.7131d0
        B_hen(4)=5.3727d0
        C_hen(4)=5.4227d0
        !ここはH2Oなので注意!
        Tr1=Temp/Tc(1)

        lnH_sat(4)=log(P1_sat)+A_hen(4)/Tr1+B_hen(4)*((1.0d0-Tr1)**0.355d0)/Tr1+C_hen(4)*dexp(1.0d0-Tr1)/Tr1**0.41d0
        do i=1,n
            lnfai_L(4,i)=lnH_sat(4)+v_L4(i)*(P(i)-P1_sat)/(R*temp)-log(P(i))
        end do
        !?------------------------------------------------------------------
        
        !do i=1,com_2phase
            !call residualvectorset4(lnfai_L0(i),eq,lnfai_L(i))
        !end do
        
        !write(*,*) lnfai_L0(1),lnfai_L0(4)
        
        
        
    end subroutine
    end module
    
        