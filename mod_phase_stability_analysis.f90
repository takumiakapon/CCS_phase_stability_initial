module mod_phase_stability_analysis
    use mod_autodiff
    use mod_condition
    use mod_input
    use mod_fugacity
    use mod_cubic_eq_d
    
    implicit none
    
    integer,private::i,j
    
    contains
    
    subroutine phase_stability_liquid(alpha0,P0,z0,g)
        !!設定用
        implicit none
        real(8),intent(inout),dimension(com_2phase)::alpha0
        real(8),intent(in)::P0,z0(com_2phase+com_ion)
        type(diffs),allocatable,intent(out)::g(:)
        type(diffs),allocatable,target::xd(:)
        type(diffs),pointer::alpha(:)
        real(8),allocatable,dimension(:)::x0
        !!計算用
        type(diffs),dimension(com_2phase)::x,y,w,lnfai_V,lnfai_L
        type(diffs)::wt
        
        real(8),allocatable,dimension(:)::kakuninn
        
        
        allocate(x0(com_2phase))
        
        !!自動微分の下準備
        do i=1,com_2phase
            x0(i) = alpha0(i)
        end do
        
        call diffsset1(x0,xd)
        call sizeset(x0,g)
        
        alpha => xd(1:com_2phase)
        call outxs(alpha,kakuninn)
        !write(*,*) kakuninn
        
        do i=1,com_2phase
            w(i) = (alpha(i) / 2.0d0) ** 2.0d0
        end do
        
        
        call outxs(w,kakuninn)
        !write(*,*) kakuninn
        !!気相のモル分率
        call residualvectorset3(com_2phase,wt)
        do i=1,com_2phase
            wt = wt + w(i)
        end do
         call outxs(alpha,kakuninn)
        !write(*,*) kakuninn
        do i=1,com_2phase
            y(i) = w(i) /wt
        end do
        
        !!液相のモル分率
        do i=1,com_2phase
            call residualvectorset4(z0(i),com_2phase,x(i))
        end do
       
        call vapor_fugacity(y,lnfai_V,P0)
        
        call liquid_fugacity(lnfai_L,P0)
        call outxs(lnfai_V,kakuninn)
        !write(*,*) kakuninn
        
        
        
        do i=1,com_2phase
            g(i) = sqrt(w(i)) * (log(w(i)) + lnfai_V(i) -log(x(i)) - lnfai_L(i))
        end do
        
    end subroutine
    
    
    subroutine phase_stability_vapor(alpha0,P0,z0,g)
        implicit none
        !!設定用
        real(8),intent(inout),dimension(com_2phase)::alpha0
        real(8),intent(in)::P0,z0(com_2phase+com_ion)
        type(diffs),allocatable,intent(out)::g(:)
        type(diffs),allocatable,target::xd(:)
        type(diffs),pointer::alpha(:)
        real(8),allocatable,dimension(:)::x0
        
        !!計算用
        type(diffs),dimension(com_2phase)::x,y,w,lnfai_L,lnfai_V
        type(diffs)::wt
        real(8),allocatable,dimension(:)::kakuninn
        
        
        allocate(x0(com_2phase))
        
        !!自動微分の下準備
        do i=1,com_2phase
            x0(i) = alpha0(i)
        end do
        
        call diffsset1(x0,xd)
        call sizeset(x0,g)
        
        alpha => xd(1:com_2phase)
        
        
        do i=1,com_2phase
            if (alpha0(i) == 0.0d0) then
                call residualvectorset3(com_2phase,w(i))
            else
                w(i) = (alpha(i) / 2.0d0) ** 2.0d0
            end if
        end do
        
        call residualvectorset3(com_2phase,wt)
        do i=1,com_2phase
            wt = wt + w(i)
        end do
        !!液相のモル分率
        do i =1,com_2phase
            x(i) = w(i) / wt
        end do
        
        !気相のモル分率
        do i=1,com_2phase
            call residualvectorset4(z0(i),com_2phase,y(i))
        end do
        
        call vapor_fugacity(y,lnfai_V,P0)
        call liquid_fugacity(lnfai_L,P0)
        
        call outxs(lnfai_L,kakuninn)
        !write(*,*) kakuninn
        
        
        
        do i=1,com_2phase
            g(i) = sqrt(w(i)) * (log(w(i)) + lnfai_L(i) -log(y(i)) - lnfai_V(i))
        end do
        
        
    end subroutine

    subroutine phase_stability_liquid2(alpha0,P0,z0,g)
        !!設定用
        implicit none
        real(8),intent(inout),dimension(com_2phase)::alpha0
        real(8),intent(in)::P0,z0(com_2phase+com_ion)
        type(diffs),allocatable,intent(out)::g(:)
        type(diffs),allocatable,target::xd(:)
        type(diffs),pointer::alpha(:)
        real(8),allocatable,dimension(:)::x0
        !!計算用
        type(diffs),dimension(com_2phase)::x,y,w,lnfai_V,lnfai_L
        type(diffs)::wt
        

        real(8)::Tc(com_2phase),Pc(com_2phase),acentric(com_2phase),bic(com_2phase,com_2phase),a(com_2phase),b(com_2phase)
        real(8)::z_factor0
        type(diffs)::a_mix_V,b_mix_V,A_V,B_V,sigma,D_V(com_2phase),E_V(com_2phase),a_coef(3),z_factor,C_V,G_V
        real(8),allocatable,dimension(:)::kakuninn
        
        
        allocate(x0(com_2phase))
        
        !!自動微分の下準備
        do i=1,com_2phase
            x0(i) = alpha0(i)
        end do
        
        
        call diffsset1(x0,xd)
        call sizeset(x0,g)
        
        alpha => xd(1:com_2phase)
        call outxs(alpha,kakuninn)
        write(*,*) kakuninn
        
        do i=1,com_2phase
            w(i) = (alpha(i) / 2.0d0) ** 2.0d0
        end do
        
        
        call outxs(w,kakuninn)
        !write(*,*) kakuninn
        !!気相のモル分率
        call residualvectorset3(com_2phase,wt)
        do i=1,com_2phase
            wt = wt + w(i)
        end do
         call outxs(alpha,kakuninn)
        !write(*,*) kakuninn
        do i=1,com_2phase
            y(i) = w(i) /wt
        end do
        call outxs(y,kakuninn)
        !write(*,*) sum(kakuninn)
        
        !!液相のモル分率
        do i=1,com_2phase
            call residualvectorset4(z0(i),com_2phase,x(i))
        end do
        !write(*,*) 'a'
        !call vapor_fugacity(y,lnfai_V,P0) 
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
        do j=1,com_2phase
            do i=1,com_2phase
                a_mix_V=a_mix_V+y(j)*y(i)*(1.0d0-bic(j,i))*sqrt(a(j)*a(i))
            end do
        end do
        !call outxs(,kakuninn)
        !write(*,*) kakuninn
        
        A_V          = (p0/((R*temp)**2.0d0)) * a_mix_V
        B_V          = (p0/(R*temp))        * b_mix_V
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
            call residualvectorset3(com_2phase,E_V(j))
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
            
        !    write(*,*) 'n'
            call outxs(E_V,kakuninn)
            !write(*,*) kakuninn
        end do
        
        G_V=(z_factor+B_V*(1.0d0+sqrt(2.0d0)))/(z_factor+B_V*(1.0d0-sqrt(2.0d0)))
        do j=1,com_2phase
            lnfai_V(j)=D_V(j)-log(C_V)-E_V(j)*log(G_V)
        end do
        
        call outxs(lnfai_V,kakuninn)
        !write(*,*) kakuninn(1),kakuninn(4)
        !write(*,*) z_factor0
        !write(*,*) 'o'


        !write(*,*) 'i'
        call liquid_fugacity(lnfai_L,P0)
        call outxs(lnfai_V,kakuninn)
        !write(*,*) kakuninn
        
        
        
        do i=1,com_2phase
            g(i) = sqrt(w(i)) * (log(w(i)) + lnfai_V(i) -log(x(i)) - lnfai_L(i))
        end do
        
    end subroutine
    
    
    subroutine phase_stability_vapor2(alpha0,P0,z0,g)
        implicit none
        !!設定用
        real(8),intent(inout),dimension(com_2phase)::alpha0
        real(8),intent(in)::P0,z0(com_2phase+com_ion)
        type(diffs),allocatable,intent(out)::g(:)
        type(diffs),allocatable,target::xd(:)
        type(diffs),pointer::alpha(:)
        real(8),allocatable,dimension(:)::x0
        
        !!計算用
        type(diffs),dimension(com_2phase)::x,y,w,lnfai_L,lnfai_V
        type(diffs)::wt
        real(8),allocatable,dimension(:)::kakuninn
        
        allocate(x0(com_2phase))
        
        !!自動微分の下準備
        do i=1,com_2phase
            x0(i) = alpha0(i)
        end do
        
        call diffsset1(x0,xd)
        call sizeset(x0,g)
        
        alpha => xd(1:com_2phase)
        
        
        do i=1,com_2phase
            if (alpha0(i) == 0.0d0) then
                call residualvectorset3(com_2phase,w(i))
            else
                w(i) = (alpha(i) / 2.0d0) ** 2.0d0
            end if
        end do
        
        call residualvectorset3(com_2phase,wt)
        do i=1,com_2phase
            wt = wt + w(i)
        end do
        !!液相のモル分率
        do i =1,com_2phase
            x(i) = w(i) / wt
        end do
        
        !気相のモル分率
        do i=1,com_2phase
            call residualvectorset4(z0(i),com_2phase,y(i))
        end do
        
        call vapor_fugacity(y,lnfai_V,P0)
        call liquid_fugacity(lnfai_L,P0)
        
        call outxs(lnfai_L,kakuninn)
        !write(*,*) kakuninn
        
        
        
        do i=1,com_2phase
            g(i) = sqrt(w(i)) * (log(w(i)) + lnfai_L(i) -log(y(i)) - lnfai_V(i))
        end do
        
        
    end subroutine
    
    
end module