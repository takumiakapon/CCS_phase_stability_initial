module mod_autodiff
    implicit none

!****************************************!
!**一変数の場合の構造体宣言**************!
!****************************************!    
    type diff
        double precision,private::x,dx
    end type diff
!****************************************!
!**多変数の場合の構造体宣言**************!
!****************************************!     
    type diffs
        double precision,private::x
        double precision,allocatable,private::dx(:)
    end type diffs
!****************************************!
!**変数と微分値の入力********************!
!****************************************! 
    interface diffset
        module procedure diffset1,diffset2,diffsset2
    end interface diffset

!--足し算の拡張--------------------------!
    interface operator(+)
        module procedure diffplusdiff,rplusdiff,diffplusr,diffsplusdiffs,rplusdiffs,&
        diffsplusr
    end interface
!--掛け算の拡張--------------------------!  
    interface operator(*)
        module procedure difftimediff,rtimediff,difftimer,diffstimediffs,rtimediffs,&
        diffstimer
    end interface
!--引き算の拡張--------------------------! 
    interface operator(-)
       module procedure diffminuisdiff,rminuisdiff,diffminuisr,&
                        diffsminuisdiffs,rminuisdiffs,diffsminuisr
    end interface
!--割り算の拡張--------------------------! 
    interface operator(/)
       module procedure diffdivdiff,rdivdiff,diffdivr,&
                        rdivdiffs,diffsdivr,diffsdivdiffs
    end interface
!--指数の拡張----------------------------!
    interface operator(**)
       module procedure diffpowdiff,diffspowdiffs,diffpowr,diffspowr,&
                        diffpowint,diffspowint,rpowdiff,rpowdiffs
    end interface
!--表示の拡張----------------------------! 
    interface disp
        module procedure dispdiff,dispdiffs,dispdiffs1,disp
    end interface disp
!--指数の拡張----------------------------!    
    interface exp
        module procedure expdiff,expdiffs
    end interface exp
!--自然対数の拡張------------------------!     
    interface log
       module procedure logdiff,logdiffs
    end interface log
!--sinの拡張----------------------------!
    interface sin
        module procedure diffsin,diffssin
    end interface sin
!--cosの拡張----------------------------!    
    interface cos
        module procedure diffcos,diffscos
    end interface cos
!--tanの拡張----------------------------!     
    interface tan
        module procedure difftan,diffstan
    end interface tan
!--平方根の拡張--------------------------!    
    interface sqrt
        module procedure diffsqrt,diffssqrt
    end interface sqrt
!--常用対数の拡張------------------------!    
    interface log10
        module procedure difflog10,diffslog10
    end interface log10
!--asinの拡張----------------------------!       
    interface asin
        module procedure diffasin,diffsasin
    end interface asin
!--acosの拡張----------------------------!        
    interface acos
        module procedure diffacos,diffsacos
    end interface acos
!--atanの拡張----------------------------!       
    interface atan
        module procedure diffatan,diffsatan
    end interface atan
!--sinhの拡張----------------------------!       
    interface sinh
        module procedure diffsinh,diffssinh
    end interface sinh
!--coshの拡張----------------------------!      
    interface cosh
        module procedure diffcosh,diffscosh
    end interface cosh
!--tanhの拡張----------------------------!      
    interface tanh
        module procedure difftanh,diffstanh
    end interface tanh
   
    interface outx
        module procedure outx
    end interface outx
 
contains
!****************************************!
!**1変数の場合***************************!
!****************************************!
    
!--主要変数の値と微分値の入力------------!    
    function diffset1(c) result(f)
        double precision,intent(in)::c
        type(diff)::f
        f%x=c
        f%dx=1d0
    end function diffset1
!--変数の値と微分値の入力----------------!   
    function diffset2(c1,c2) result(f)
        double precision,intent(in)::c1,c2
        type(diff)::f
        f%x=c1
        f%dx=c2
    end function diffset2
!--f(x)+g(x)-----------------------------! 
    function diffplusdiff(f1,f2) result(g)
        type(diff),intent(in)::f1,f2
        type(diff)::g
        g=diffset(f1%x+f2%x,f1%dx+f2%dx)
    end function diffplusdiff
!--変数の値と微分値の表示----------------!    
    subroutine dispdiff(f)
        type(diff),intent(in)::f
        write(*,*)'f%x:',f%x,'f%dx:',f%dx
    end subroutine dispdiff
!--c+f(x)--------------------------------!     
    function rplusdiff(c,f) result(g)
        double precision,intent(in)::c
        type(diff),intent(in)::f
        type(diff)::g
        g=diffset(c+f%x,f%dx)
    end function rplusdiff
!--f(x)+c--------------------------------!      
    function diffplusr(f,c) result(g)
        double precision,intent(in)::c
        type(diff),intent(in)::f
        type(diff)::g
        g=diffset(f%x+c,f%dx)
    end function diffplusr
!--f(x)*g(x)-----------------------------!     
    function difftimediff(f1,f2) result(g)
        type(diff),intent(in)::f1,f2
        type(diff)::g
        g=diffset(f1%x*f2%x,f1%dx*f2%x+f1%x*f2%dx)
    end function difftimediff
 !--c*f(x)--------------------------------!     
    function rtimediff(c,f) result(g)
        double precision,intent(in)::c
        type(diff),intent(in)::f
        type(diff)::g
        g=diffset(c*f%x,c*f%dx)
    end function rtimediff
 !--f(x)*c--------------------------------!   
    function difftimer(f,c) result(g)
        double precision,intent(in)::c
        type(diff),intent(in)::f
        type(diff)::g
        g=diffset(f%x*c,f%dx*c)
    end function difftimer
!--exp(f(x))------------------------------!    
    function expdiff(f) result(g)
        type(diff),intent(in)::f
        type(diff)::g
        g=diffset(exp(f%x),f%dx*exp(f%x))
    end function expdiff
!--変数の値(x)を返す----------------------!    
    function outx(f) result(x)
        type(diff),intent(in)::f
        double precision::x
        x=f%x
    end function outx
!--微分値(dx)を返す-----------------------!     
    function outdx(f) result(dx)
        type(diff),intent(in)::f
        double precision::dx
        dx=f%dx
    end function outdx
!--f(x)-g(x)-----------------------------!    
    function diffminuisdiff(f1,f2) result(g)
        type(diff),intent(in)::f1,f2
        type(diff)::g
        g=diffset(f1%x-f2%x,f1%dx-f2%dx)
    end function diffminuisdiff
!--c-f(x)--------------------------------!      
    function rminuisdiff(c,f) result(g)
        double precision,intent(in)::c
        type(diff),intent(in)::f
        type(diff)::g
        g=diffset(c-f%x,-f%dx)
    end function rminuisdiff
!--f(x)-c--------------------------------!     
    function diffminuisr(f,c) result(g)
        double precision,intent(in)::c
        type(diff),intent(in)::f
        type(diff)::g
        g=diffset(f%x-c,f%dx)
    end function diffminuisr
!--c/f(x)---------------------------------!      
    function rdivdiff(c,f) result(g)
      double precision,intent(in)::c
      type(diff),intent(in)::f
      type(diff)::g
      g=diffset(c/f%x,-1d0*c*f%dx/(f%x*f%x))
    end function rdivdiff
!--f(x)/c---------------------------------!     
    function diffdivr(f,c) result(g)
      double precision,intent(in)::c
      type(diff),intent(in)::f
      type(diff)::g
      g=diffset(f%x/c,f%dx/c)
    end function diffdivr
!--f(x)/g(x)------------------------------!     
    function diffdivdiff(f1,f2) result(g)
      type(diff),intent(in)::f1,f2
      type(diff)::g
      g=diffset(f1%x/f2%x,(f1%dx*f2%x-f1%x*f2%dx)/(f2%x*f2%x))
    end function diffdivdiff
!--ln(f(x))-------------------------------!     
    function logdiff(f) result(g)
      type(diff),intent(in)::f
      type(diff)::g
      g=diffset(log(f%x),f%dx/f%x)
    end function logdiff
!--f(x)**g(x)-----------------------------!     
    function diffpowdiff(f1,f2) result(g)
      type(diff),intent(in)::f1,f2
      type(diff)::g
      g=diffset(f1%x**f2%x,(f1%x**f2%x)*(f2%dx*log(f1%x)+f2%x*f1%dx/f1%x))
    end function diffpowdiff
!--f(x)**c（実数）------------------------!       
    function diffpowr(f,c) result(g)
      type(diff),intent(in)::f
      double precision,intent(in)::c
      type(diff)::g
      g=diffset(f%x**c,c*f%dx*f%x**(c-1d0))
    end function diffpowr
!--f(x)**c（整数）------------------------!     
    function diffpowint(f,c) result(g)
      type(diff),intent(in)::f
      integer,intent(in)::c
      type(diff)::g
      g=diffset(f%x**c,dble(c)*f%dx*f%x**(c-1))
    end function diffpowint
!--c**f(x)--------------------------------!         
    function rpowdiff(c,f) result(g)
      type(diff),intent(in)::f
      double precision,intent(in)::c
      type(diff)::g
      g=exp(f*log(c))
     !g=diffset(c**f%x,c**f%x*f%dx*log(c))
    end function rpowdiff
!--sin(f(x))------------------------------! 
    function diffsin(f) result(g)
        type(diff),intent(in)::f
        type(diff)::g
        g=diffset(sin(f%x),f%dx*cos(f%x))
    end function diffsin
!--cos(f(x))------------------------------!     
    function diffcos(f) result(g)
        type(diff),intent(in)::f
        type(diff)::g
        g=diffset(cos(f%x),-1d0*f%dx*sin(f%x))
    end function diffcos
!--tan(f(x))------------------------------!     
    function difftan(f) result(g)
        type(diff),intent(in)::f
        type(diff)::g
        g=sin(f)/cos(f)
       !g=diffset(sin(f%x)/cos(f%x),f%dx/cos(f%x)**(2d0))
    end function difftan
!--sqrt(f(x))------------------------------!        
    function diffsqrt(f) result(g)
        type(diff),intent(in)::f
        type(diff)::g
        g=f**0.5d0
       !g=diffset(f%x**0.5d0,0.5d0*f%dx*f%x**(-0.5d0))
    end function diffsqrt
!--log10(f(x))------------------------------!        
    function difflog10(f) result(g)
        type(diff),intent(in)::f
        type(diff)::g
        g=log(f)/log(10d0)
    end function difflog10
!--asin(f(x))------------------------------!    
    function diffasin(f) result(g)
        type(diff),intent(in)::f
        type(diff)::g
        g=diffset(asin(f%x),f%dx/sqrt(1d0-f%x**2))
    end function diffasin
!--acos(f(x))------------------------------!        
    function diffacos(f) result(g)
        type(diff),intent(in)::f
        type(diff)::g
        g=diffset(acos(f%x),-f%dx/sqrt(1d0-f%x**2))
    end function diffacos
!--atan(f(x))------------------------------!    
    function diffatan(f) result(g)
        type(diff),intent(in)::f
        type(diff)::g
        g=diffset(atan(f%x),f%dx/(1d0+f%x**2))
    end function diffatan
!--sinh(f(x))------------------------------!    
    function diffsinh(f) result(g)
        type(diff),intent(in)::f
        type(diff)::g
        g=diffset(sinh(f%x),f%dx*cosh(f%x))
    end function diffsinh
!--cosh(f(x))------------------------------!      
    function diffcosh(f) result(g)
        type(diff),intent(in)::f
        type(diff)::g
        g=diffset(cosh(f%x),f%dx*sinh(f%x))
    end function diffcosh
!--tanh(f(x))------------------------------!     
    function difftanh(f) result(g)
        type(diff),intent(in)::f
        type(diff)::g
        g=sinh(f)/cosh(f)
    end function difftanh

!****************************************!
!**多変数の場合**************************!
!****************************************!

 ! 多変数関数を数値的に解く際のresidual vectorの大きさを決定するためのサブルーチン
    subroutine residualvectorset(c,f)
        double precision,intent(in)::c(:)
        type(diffs),allocatable,intent(out)::f(:)
        integer::n,i
        n=size(c)
        if(allocated(f))then
            deallocate(f)
        end if
        allocate(f(n))
        do i=1,n
            if(allocated(f(i)%dx))then
                deallocate(f(i)%dx)
            end if
            allocate(f(i)%dx(n))
        end do
    end subroutine residualvectorset

    subroutine residualvectorset2(n,f)
        integer,intent(in)::n
        type(diffs),intent(out)::f

            if(allocated(f%dx))then
                deallocate(f%dx)
            end if
            allocate(f%dx(n))
    end subroutine residualvectorset2

    subroutine residualvectorset3(n,f)
        integer,intent(in)::n
        type(diffs),intent(out)::f
        integer :: i
            f%x = 0d0
            if(allocated(f%dx))then
                deallocate(f%dx)
            end if
            allocate(f%dx(n))
            do i=1,n
                f%dx(i)     =0d0
            end do

    end subroutine residualvectorset3

    subroutine residualvectorset4(c,n,f)
        double precision,intent(in)::c
        integer,intent(in)::n
        type(diffs),intent(out)::f
        integer :: i
            if(allocated(f%dx))then
                deallocate(f%dx)
            end if
            allocate(f%dx(n))
            do i=1,n
                f%dx(i)     =0d0
            end do
            f%x = c

    end subroutine residualvectorset4

    subroutine residualvectorset5(c,f)
        double precision,intent(in)::c
        type(diffs),intent(inout)::f

            f%x = c

    end subroutine residualvectorset5

    subroutine residualvectorset6(c,g,f)
        double precision,intent(in)::c
        type(diffs),intent(in)::g
        type(diffs),intent(out)::f
        integer :: i,n,m
            n=size(g%dx)
            if(allocated(f%dx))then
                deallocate(f%dx)
            end if
            allocate(f%dx(n))
            m=1
            do i=1,n
            if (g%dx(i)>0) exit
            m=m+1
            end do
            do i=1,n
                f%dx(i)     =0d0
            end do
            f%dx(m)     =1d0
            f%x = c

    end subroutine residualvectorset6

!--主要変数の値と微分値の入力------------!  
    subroutine diffsset1(c,f)
        double precision,intent(in)::c(:)
        type(diffs),allocatable,intent(out)::f(:)
        integer::n,i
        n=size(c)!nは主要変数の数!
        allocate(f(n))
        do i=1,n
            allocate(f(i)%dx(n))
            f(i)%dx     =0d0
            f(i)%dx(i)  =1d0
            f(i)%x      =c(i)
        end do
    end subroutine diffsset1

    subroutine diffsset3(n,c,f)
        integer,intent(in)::n
        double precision,intent(in)::c
        type(diffs),intent(inout)::f
        integer::i
        
        f%dx     =0d0
        f%dx(n)  =1d0
        f%x      =c
    end subroutine diffsset3
!--変数の値と微分値の入力----------------!      
    function diffsset2(c1,c2) result(f)
        double precision,intent(in)::c1,c2(:)
        type(diffs)::f
        integer::n
        f%x=c1
        n=size(c2)
        if (allocated(f%dx)) then
            deallocate(f%dx)
        end if
        allocate(f%dx(n))
        f%dx=c2
    end function diffsset2
!--変数の値と微分値(0)の入力----------------!      
!    function diffsset4(c1) result(f)
!        double precision,intent(in)::c1
!        type(diffs)::f
!        integer::n,i
!        n=size(f%dx)
!        f%x=c1
!        do i=1,n
!        f%dx(i)=0d0
!        end do
!    end function diffsset4
!--主要変数の数に対応した新たな変数のサイズ決定--!   
    subroutine sizeset(f,g)
        double precision,intent(in)::f(:)
        type(diffs),intent(out),allocatable::g(:)
        integer::n,i
        n=size(f)
        if (allocated(g)) then
            deallocate(g)
        end if
        allocate(g(n))
        do i=1,n
            if (allocated(g(i)%dx)) then
                deallocate(g(i)%dx)
            end if
            allocate(g(i)%dx(n))
        end do
    end subroutine sizeset
!--主要変数の値と微分値の表示------------!
    subroutine dispdiffs(f)
        type(diffs),intent(in)::f(:)
        integer::n,i,j,m
        n=size(f)

        do i=1,n
            write(*,'(a2,i3,a10,f22.15)')'f(',i,')%x     :',f(i)%x
            m=size(f(i)%dx)
            do j=1,m
                write(*,'(a2,i3,a5,i3,a2,f22.15)')'f(',i,')%dx(',j,'):',f(i)%dx(j)
            end do
        end do
    end subroutine dispdiffs
!--変数の値と微分値の表示----------------!        
    subroutine dispdiffs1(f)
        type(diffs),intent(in)::f
        integer::n,j
        n=size(f%dx)

        write(*,'(a2,1x,a9,f22.15)')'f','%x      :',f%x
        do j=1,n
            write(*,'(a2,a5,i3,a2,f22.15)')'f','%dx(',j,'):',f%dx(j)
        end do
    end subroutine dispdiffs1
!--f(x,y)+g(x,y)-----------------------------! 
    function diffsplusdiffs(f1,f2) result(g)
        type(diffs),intent(in)::f1,f2
        type(diffs)::g
        g=diffset(f1%x+f2%x,f1%dx+f2%dx)
    end function diffsplusdiffs
!--c+f(x,y)--------------------------------!    
    function rplusdiffs(c,f) result(g)
        double precision,intent(in)::c
        type(diffs),intent(in)::f
        type(diffs)::g
        g=diffset(c+f%x,f%dx)
    end function rplusdiffs
!--f(x,y)+c--------------------------------!     
    function diffsplusr(f,c) result(g)
        double precision,intent(in)::c
        type(diffs),intent(in)::f
        type(diffs)::g
        g=diffset(f%x+c,f%dx)
    end function diffsplusr
!--f(x,y)*g(x,y)-----------------------------!    
    function diffstimediffs(f1,f2) result(g)
        type(diffs),intent(in)::f1,f2
        type(diffs)::g
        g=diffset(f1%x*f2%x,f1%dx*f2%x+f1%x*f2%dx)
    end function diffstimediffs
!--c*f(x,y)--------------------------------!        
    function rtimediffs(c,f) result(g)
        double precision,intent(in)::c
        type(diffs),intent(in)::f
        type(diffs)::g
        g=diffset(c*f%x,c*f%dx)
    end function rtimediffs
!--f(x,y)*c--------------------------------!      
    function diffstimer(f,c) result(g)
        double precision,intent(in)::c
        type(diffs),intent(in)::f
        type(diffs)::g
        g=diffset(f%x*c,f%dx*c)
    end function diffstimer
!--exp(f(x,y))------------------------------!       
    function expdiffs(f) result(g)
        type(diffs),intent(in)::f
        type(diffs)::g
        g=diffset(exp(f%x),f%dx*exp(f%x))
    end function expdiffs
!--f(x,y)-g(x,y)-----------------------------!        
    function diffsminuisdiffs(f1,f2) result(g)
        type(diffs),intent(in)::f1,f2
        type(diffs)::g
        g=diffset(f1%x-f2%x,f1%dx-f2%dx)
    end function diffsminuisdiffs
!--c-f(x,y)--------------------------------!        
    function rminuisdiffs(c,f) result(g)
        double precision,intent(in)::c
        type(diffs),intent(in)::f
        type(diffs)::g
        g=diffset(c-f%x,-f%dx)
    end function rminuisdiffs
!--f(x,y)-c--------------------------------!       
    function diffsminuisr(f,c) result(g)
        double precision,intent(in)::c
        type(diffs),intent(in)::f
        type(diffs)::g
        g=diffset(f%x-c,f%dx)
    end function diffsminuisr
!--c/f(x,y)---------------------------------!       
    function rdivdiffs(c,f) result(g)
      double precision,intent(in)::c
      type(diffs),intent(in)::f
      type(diffs)::g
      g=diffset(c/f%x,-1d0*c*f%dx/(f%x*f%x))
    end function rdivdiffs
!--f(x,y)/c---------------------------------!      
    function diffsdivr(f,c) result(g)
      double precision,intent(in)::c
      type(diffs),intent(in)::f
      type(diffs)::g
      g=diffset(f%x/c,f%dx/c)
    end function diffsdivr
!--f(x,y)/g(x,y)------------------------------!      
    function diffsdivdiffs(f1,f2) result(g)
      type(diffs),intent(in)::f1,f2
      type(diffs)::g
      g=diffset(f1%x/f2%x,(f1%dx*f2%x-f1%x*f2%dx)/(f2%x*f2%x))
    end function diffsdivdiffs
!--ln(f(x,y))-------------------------------!       
    function logdiffs(f) result(g)
      type(diffs),intent(in)::f
      type(diffs)::g
      g=diffset(log(f%x),f%dx/f%x)
    end function logdiffs
    
    function diffspowdiffs(f1,f2) result(g)
      type(diffs),intent(in)::f1,f2
      type(diffs)::g
      g=diffset(f1%x**f2%x,(f1%x**f2%x)*(f2%dx*log(f1%x)+f2%x*f1%dx/f1%x))
    end function diffspowdiffs
!--f(x,y)**g(x,y)-----------------------------!      
    function diffspowr(f,c) result(g)
      type(diffs),intent(in)::f
      double precision,intent(in)::c
      type(diffs)::g
      g=diffset(f%x**c,c*f%dx*f%x**(c-1d0))
    end function diffspowr
!--f(x,y)**c---------------------------------!       
    function diffspowint(f,c) result(g)
      type(diffs),intent(in)::f
      integer,intent(in)::c
      type(diffs)::g
      g=diffset(f%x**c,dble(c)*f%dx*f%x**(c-1))
    end function diffspowint
!--c**f(x,y)-----------------------------!       
    function rpowdiffs(c,f) result(g)
      type(diffs),intent(in)::f
      double precision,intent(in)::c
      type(diffs)::g
      g=exp(f*log(c))
    end function rpowdiffs
!--sin(f(x,y))------------------------------!     
    function diffssin(f) result(g)
        type(diffs),intent(in)::f
        type(diffs)::g
        g=diffset(sin(f%x),f%dx*cos(f%x))
    end function diffssin
!--cos(f(x,y))------------------------------!     
    function diffscos(f) result(g)
        type(diffs),intent(in)::f
        type(diffs)::g
        g=diffset(cos(f%x),-1d0*f%dx*sin(f%x))
    end function diffscos
!--tan(f(x,y))------------------------------!     
    function diffstan(f) result(g)
        type(diffs),intent(in)::f
        type(diffs)::g
        g=sin(f)/cos(f)
    end function diffstan
!--sqrt(f(x,y))------------------------------!     
    function diffssqrt(f) result(g)
        type(diffs),intent(in)::f
        type(diffs)::g
        g=f**0.5d0
    end function diffssqrt
!--log10(f(x,y))------------------------------!      
    function diffslog10(f) result(g)
        type(diffs),intent(in)::f
        type(diffs)::g
        g=log(f)/log(10d0)
    end function diffslog10
!--asin(f(x,y))------------------------------!      
    function diffsasin(f) result(g)
        type(diffs),intent(in)::f
        type(diffs)::g
        g=diffset(asin(f%x),f%dx/sqrt(1d0-f%x**2))
    end function diffsasin
!--acos(f(x,y))------------------------------!     
    function diffsacos(f) result(g)
        type(diffs),intent(in)::f
        type(diffs)::g
        g=diffset(acos(f%x),-f%dx/sqrt(1d0-f%x**2))
    end function diffsacos
!--atan(f(x,y))------------------------------!    
    function diffsatan(f) result(g)
        type(diffs),intent(in)::f
        type(diffs)::g
        g=diffset(atan(f%x),f%dx/(1d0+f%x**2))
    end function diffsatan
!--sinh(f(x,y))------------------------------!       
    function diffssinh(f) result(g)
        type(diffs),intent(in)::f
        type(diffs)::g
        g=diffset(sinh(f%x),f%dx*cosh(f%x))
    end function diffssinh
!--cosh(f(x,y))------------------------------!         
    function diffscosh(f) result(g)
        type(diffs),intent(in)::f
        type(diffs)::g
        g=diffset(cosh(f%x),f%dx*sinh(f%x))
    end function diffscosh
!--tanh(f(x,y))------------------------------!         
    function diffstanh(f) result(g)
        type(diffs),intent(in)::f
        type(diffs)::g
        g=sinh(f)/cosh(f)
    end function diffstanh

!--構造体からヤコビ行列作成------------------!   
    subroutine jacobian_mat(f,g)
        type(diffs),intent(in)::f(:)
        double precision,allocatable,intent(out)::g(:,:)
        integer::n,i,j
        if (allocated(g)) then
            deallocate(g)
        end if
        n=size(f)
        allocate(g(n,n))
        do i=1,n
            do j=1,n
                g(i,j)=f(i)%dx(j)
            end do
        end do
    end subroutine jacobian_mat
!--変数の値(x)を返す----------------------!        
    subroutine out_diffsx(f,x)
        type(diffs),intent(in)::f
        double precision,intent(out)::x
        x=f%x
    end subroutine out_diffsx
!--主要変数の値(x)を返す------------------!    
    subroutine outxs(f,x)
        type(diffs),intent(in)::f(:)
        double precision,allocatable,intent(out)::x(:)
        integer::n
        n=size(f)
        if (allocated(x)) then
            deallocate(x)
        end if
        allocate(x(n))
        x=f%x
    end subroutine outxs
!--微分値(dx)を返す-----------------------!        
    subroutine outdxs(f,dx)
        type(diffs),intent(in)::f
        double precision,allocatable,intent(out)::dx(:)
        if (allocated(dx)) then
            deallocate(dx)
        end if
        allocate(dx(size(f%dx)))
        dx=f%dx
    end subroutine outdxs
!--微分値(dx)の総数を返す-----------------------! 
    subroutine outndx(f,n)
        type(diffs),intent(in)::f
        integer,intent(out)::n
        n = size(f%dx)
    end subroutine outndx
    !diffm型f%x取り出し
    subroutine outdiffmx(f,fx)
        type(diffs),intent(in)::f(:)
        double precision,allocatable,intent(out)::fx(:)
        integer::n
        n=size(f)
        if(allocated(fx))then
            deallocate(fx)
        end if
        allocate(fx(n))
        fx=f%x
    end subroutine outdiffmx
!--行列の書き出し-------------------------!       
    subroutine disp(a)
        double precision,intent(in)::a(:,:)
        integer::n,m,i,j
        n=size(a,1)
        m=size(a,2)
        do i=1,n
           do j=1,m-1
              write(11,'(e15.5)',advance='no')a(i,j) 
           end do
           write(11,'(e15.5)')a(i,m) 
        end do
        return
    end subroutine disp
!--構造体から行列の値を計算する---------------!      
    function jacobian_det(f) result(g)
        type(diffs),intent(in)::f(:)
        double precision::g
        double precision,allocatable::t(:,:)
        call jacobian_mat(f,t)        
        g=determinant(t)
    end function jacobian_det
!--行列の値を計算する-------------------------!           
    recursive function determinant(A) result(f)
        double precision,intent(in)::A(:,:)
        double precision::f
        integer::n,m,i
        double precision,allocatable::t(:,:)
        n=size(A,1)
        m=size(A,2)
        allocate(t(n-1,n-1))
        if (n/=m) then
           write(*,*)'Error : determinant : n/=m'
           stop
        end if
        if (n==2) then
           f=A(1,1)*A(2,2)-A(1,2)*A(2,1)
        else
           f=0
           do i=1,n
              if (i==1) then
                 t=a(2:n,2:n)
              else if (i==n) then
                 t=a(2:n,1:n-1)
              else
                 t(:,1:i-1)=A(2:n,1:i-1)
                 t(:,i:n-1)=A(2:n,i+1:n)
              end if
              f=f+A(1,i)*determinant(t)*(-1)**(i+1)
           end do
        end if
  end function determinant

      !Jacobian作成用サブルーチン
    subroutine makeJacobian(g,a)
        type(diffs),intent(in)::g(:)
        double precision,allocatable,intent(out)::a(:,:)
        integer::n,i,j
        if(allocated(a))then
            deallocate(a)
        end if
        n=size(g)
        allocate(a(n,n))
        do i=1,n
            do j=1,n
                a(i,j)=g(i)%dx(j)
            end do
        end do
    end subroutine makeJacobian
    
end module mod_autodiff
