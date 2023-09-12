module mod_ini_N
    use mod_condition
    use mod_input
    use mod_cubic_eq
    
    implicit none
    contains
    
    subroutine ini_N_2phase(z_factor0,V0,z0,x0,y0,k0,Nc0,Nt,Sw)
        implicit none
        real(8),intent(in)::z_factor0,V0
        real(8),intent(inout),dimension(com_2phase+com_ion)::z0,x0,Nc0
        real(8),intent(inout),dimension(com_2phase)::k0,y0
        real(8),intent(out)::Nt,Sw
        
        
        integer::i
        real(8)::MD_V,MV_V,MD_L,MV_L,A_rc,B_rc,C_rc,P0
        real(8),dimension(com_2phase+com_ion)::v_L
        
        P0 = iniPressure
        
        !モル分率
        do i=1,com_2phase
            x0(i) = z0(i)/(1.0d0-V0+V0*k0(i))
            y0(i) = k0(i) *x0(i)
        end do
        do i=com_2phase +1, com_2phase + com_ion
            x0(i) = z0(i) /(1.0d0-V0)
        end do
        
        
        !!気相について
        MD_V = P0/z_factor0/R/temp !モル密度
        MV_V = 1.0d0 / MD_V !モル体積
        
        !!液相について
        A_rc=5.916365d0-0.01035794d0*temp+0.9270048d0*(10.0d0**(-5.0d0))*(temp**2.0d0)-1127.522d0/temp+100674.1d0/(temp**2.0d0)!!
        B_rc=0.5204914d0*(10.0d0**(-2.0d0))-0.10482101d0*(10.0d0**(-4.0d0))*temp+0.8328532d0*(10.0d0**(-8.0d0))*(temp**2.0d0)&!!
            -1.1702939d0/temp+102.2783d0/(temp**2.0d0)!!
        C_rc=0.118547d0*(10.0d0**(-7.0d0))-0.6599143d0*(10.0d0**(-10.0d0))*temp!!
        v_L(1)=(M1/10.0d0**3.0d0)*(A_rc-B_rc*(iniPressure/(9.80665d0*10.0d0**4.0d0))-&
                C_rc*(iniPressure/(9.80665d0*10.0d0**4.0d0))**2.0d0)!!
        v_L(2)=(-47.7518+4.336154/10.0d0*temp-5.945771/10.0d0**4.0d0*temp**2.0d0)/10.0d0**6.0d0
        v_L(3)=exp(3.541d0+1.23d-3*(temp-273.15))*10**(-6.d0) !CH4
        v_L(4)=(160.5567d0-5.538776d0*10.0d0**(-1.0d0)*temp)*10**(-6.d0)
        do i=com_2phase + 1, com_2phase + com_ion
            v_L(i) = 0.0d0!v_L(1)/10.0d0
        end do
        
        MV_L = 0.0d0
        do i=1,com_2phase +com_ion
            MV_L = MV_L +x0(i) *v_L(i)
        end do
        MD_L = 1.0d0/MV_L
        
        !!モル数
        Nt = 1.0d0/(V0*MV_V+(1.0d0-V0)*MV_L)
        do i= 1,com_2phase +com_ion
            Nc0(i) = Nt*z0(i)
        end do
        
        !!飽和率
        Sw = Nt*(1.0d0-V0)*MV_L
        
        !write(*,*) '成分','液相', '気相','モル数'
        do i =1,com_2phase
        !    write(*,*) i,x0(i),y0(i),Nc0(i)
        end do
        do i=com_2phase +1, com_2phase+com_ion
        !    write(*,*) i, x0(i),Nc0(i)
        end do
        
        !write(*,*) y0
    end subroutine
    
    subroutine ini_N_liquid(V0,z0,x0,y0,Nc0,Nt,Sw)
        implicit none
        real(8),intent(in)::V0
        real(8),intent(inout),dimension(com_2phase+com_ion)::z0,Nc0,x0
        real(8),intent(inout),dimension(com_2phase)::y0
        real(8),intent(out)::Nt,Sw
        
        integer::i
        real(8)::MD_L,MV_L,A_rc,B_rc,C_rc,L0,Sg
        real(8),dimension(com_2phase+com_ion)::v_L
        
        L0 = 1.0d0 - V0
        
        
        !!液相について
        A_rc=5.916365d0-0.01035794d0*temp+0.9270048d0*(10.0d0**(-5.0d0))*(temp**2.0d0)-1127.522d0/temp+100674.1d0/(temp**2.0d0)!!
        B_rc=0.5204914d0*(10.0d0**(-2.0d0))-0.10482101d0*(10.0d0**(-4.0d0))*temp+0.8328532d0*(10.0d0**(-8.0d0))*(temp**2.0d0)&!!
            -1.1702939d0/temp+102.2783d0/(temp**2.0d0)!!
        C_rc=0.118547d0*(10.0d0**(-7.0d0))-0.6599143d0*(10.0d0**(-10.0d0))*temp!!
        v_L(1)=(M1/10.0d0**3.0d0)*(A_rc-B_rc*(iniPressure/(9.80665d0*10.0d0**4.0d0))-&
                C_rc*(iniPressure/(9.80665d0*10.0d0**4.0d0))**2.0d0)!!
        v_L(2)=(-47.7518+4.336154/10.0d0*temp-5.945771/10.0d0**4.0d0*temp**2.0d0)/10.0d0**6.0d0
        v_L(3)=exp(3.541d0+1.23d-3*(temp-273.15))*10**(-6.d0) !CH4
        v_L(4)=(160.5567d0-5.538776d0*10.0d0**(-1.0d0)*temp)*10**(-6.d0)
        do i=com_2phase + 1, com_2phase + com_ion
            v_L(i) = 0.0d0!v_L(2)/10.0d0
        end do
        
        
        
        MV_L = 0.0d0
        do i=1,com_2phase +com_ion
            MV_L = MV_L +z0(i) *v_L(i)
        end do
        MD_L = 1.0d0/MV_L
        
        x0=z0
        y0 = 0.0d0
        !!モル数--------
        Nt=1.0d0/(L0*MV_L)
        do i=1,com_2phase+com_ion
            Nc0(i)=Nt*x0(i)
        end do
        
        !!飽和率-------
        Sw=Nt*L0*MV_L
        Sg =1.0d0-Sw
        
        !write(*,*) '成分','液相','気相','モル数'
        do i =1,com_2phase
        !    write(*,*) i,x0(i),y0(i),Nc0(i)
        end do
        do i=com_2phase +1, com_2phase+com_ion
        !    write(*,*) i, x0(i),'なし',Nc0(i)
        end do
        
    end subroutine
    
    subroutine ini_N_vapor(V0,z0,x0,y0,Nc0,Nt,Sw)
        implicit none
        real(8),intent(in)::V0
        real(8),intent(inout),dimension(com_2phase+com_ion)::z0,x0,Nc0
        real(8),intent(inout),dimension(com_2phase)::y0
        real(8),intent(inout)::Nt,Sw
        
        real(8),dimension(com_2phase)::Tc,Pc,acentric,a,b
        real(8),dimension(com_2phase,com_2phase)::bic
        real(8)::L0,a_mix_V,b_mix_V,A_V,B_V,a_coef(3),z_factor,p0,MD_V,MV_V,Sg,p
        
        integer::i,j
        L0 =1.0d0-V0
        
        p=inipressure
        x0=0.0d0
        do i=1,com_2phase
            y0(i) = z0(i)
        end do
        
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
        
        b_mix_V = 0.0d0
        a_mix_V= 0.0d0
        do j=1,com_2phase
            a(j) = (0.457235528921d0*(R**2.0d0)*(tc(j)**2.0d0)/pc(j)) * ((1.0d0+(0.37464d0+1.54226d0*acentric(j)-0.26992*&
            (acentric(j)**2.0d0))*(1.0d0-sqrt(temp/tc(j))) ))**2.0d0
            b(j) = 0.0777960739039d0*R*tc(j)/pc(j)
            b_mix_V=b_mix_V+y0(j)*b(j)
        end do
        !write(*,*) a(4)
        
        
        do j=1,com_2phase
            do i=1,com_2phase
                a_mix_V=a_mix_V+y0(j)*y0(i)*(1.0d0-bic(j,i))*sqrt(a(j)*a(i))
            end do
        end do
        !call outxs(,kakuninn)
        !write(*,*) kakuninn
        
        A_V          = (p/((R*temp)**2.0d0)) * a_mix_V
        B_V          = (p/(R*temp))        * b_mix_V
        a_coef(1)       =(B_V-1.0d0)
        a_coef(2)       = A_V - 2.0d0*B_V - 3.0d0*B_V**2.0d0
        a_coef(3)       =(B_V**2.0d0 + B_V**3.0d0-A_V*B_V)
        
        call cubic_eq(a_coef(1),a_coef(2),a_coef(3),z_factor)
        
        p0 =inipressure    
        
        MD_V = P0/z_factor/R/temp
        MV_V = 1.0d0/MD_V
        
        Nt = 1.0d0/(V0*MV_V)
        do i=1,com_2phase+com_ion
            Nc0(i) = Nt*z0(i)
        end do
        
        
    
        
        
         Sg=Nt*V0*MV_V
         Sw =1.0d0-Sg
         
        ! write(*,*) '成分','液相','気相','モル数'
        do i =1,com_2phase
        !    write(*,*) i,x0(i),y0(i),Nc0(i)
        end do
        do i=com_2phase +1, com_2phase+com_ion
        !    write(*,*) i, x0(i),'なし',Nc0(i)
        end do
         
         end subroutine
    
end module
    
        
        