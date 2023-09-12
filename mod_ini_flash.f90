module mod_ini_flash
    use mod_autodiff
    use mod_condition
    use mod_input
    use mod_fugacity
    use mod_cubic_eq_d

    implicit none
    integer,private::i,j

    contains

    subroutine ini_fla(P,z0,lnk0,v0,z_factor0,g)
        implicit none
        real(8),intent(inout)::P,V0,z_factor0
        real(8),intent(inout),dimension(com_2phase)::lnk0
        real(8),intent(inout),dimension(com_2phase+com_ion)::z0
        type(diffs),allocatable,intent(out)::g(:)
        type(diffs),allocatable,target::xd(:)
        type(diffs),pointer::lnk(:),V
        real(8),allocatable,dimension(:)::x0
        
        
        type(diffs)::fra
        type(diffs),dimension(com_2phase)::k,y,lnfai_V,lnfai_L
        type(diffs),dimension(com_2phase+com_ion)::z,x
        
        real(8),allocatable,dimension(:)::kakuninn
        real(8)::kaku

        real(8)::Tc(com_2phase),Pc(com_2phase),acentric(com_2phase),bic(com_2phase,com_2phase),a(com_2phase),b(com_2phase)
        type(diffs)::a_mix_V,b_mix_V,A_V,B_V,a_coef(3),z_factor,C_V,D_V(com_2phase),E_V(com_2phase),sigma,G_V
        
        
        allocate(x0(com_2phase+1))
        
        !!自動微分の下準備----------------------------
        do i=1,com_2phase
            x0(i)=lnk0(i)
        end do
        x0(com_2phase+1) = v0

        call diffsset1(x0,xd)
        call sizeset(x0,g)

        lnk => xd(1:com_2phase)
        V => xd(com_2phase+1)
        !?-----------------------------------------
        call outxs(lnk,kakuninn)
        !write(*,*) kakuninn
        
        
        !!rachford-rice
        do i=1,com_2phase+com_ion
            call residualvectorset4(z0(i),com_2phase+1,z(i))
        end do
        call residualvectorset3(com_2phase+1,fra)
        do i=1,com_2phase
            fra=fra+(1.0d0-exp(lnk(i)))*z(i)/(1.0d0-V+V*exp(lnk(i)))
            call out_diffsx(fra,kaku)
            !write(*,*) kaku
        end do
        do i=com_2phase+1,com_2phase+com_ion
            fra = fra+z(i)/(1.0d0-V)   
        end do
        
        
            
        do i=1,com_2phase
            k(i)=exp(lnk(i))
            x(i)=z(i)/(1.d0-V+V*k(i))
            y(i)=x(i)*k(i)
        end do
        do i=com_2phase+1,com_2phase+com_ion
            x(i)=z(i)/(1.0d0-V)
        end do
        
        !call vapor_fugacity_ini(y,lnfai_V,P,z_factor0)
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
        
        call liquid_fugacity_ini(lnfai_L,P)
        
        
        do i=1,com_2phase
            g(i) = lnk(i)+lnfai_V(i)-lnfai_L(i)
        end do 
        g(com_2phase+1) = fra
        call out_diffsx(fra,kaku)
        !write(*,*) kaku
        
        
        !!相安定解析と初期フラッシュでは主要変数の数が違うから、それを「mod_fugacity」では考慮できていない→主要変数の数も引数にしちゃう？
        
        
        
        
    end subroutine
    
    
    end module
    