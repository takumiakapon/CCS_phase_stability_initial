module mod_cubic_eq_d
    use mod_autodiff
    implicit none

    contains

    subroutine cubic_eq_d(a,b,c,z)
        implicit none
        type(diffs),intent(in)::a,b,c
        type(diffs),intent(out)::z

        double precision, parameter :: pai = 4d0*atan(1d0)    
        type(diffs) :: p,q,d,theta
        double precision ::d_x
        type(diffs), allocatable, dimension(:) :: x, y, LA
        double precision, allocatable, dimension(:) :: x_x
        integer :: i,ns, nsmax,nsmin
        double precision :: xmax,xmin

        p = (-1d0)*(a**2d0)/9d0 + b/3d0
        q = 2d0*(a**3) /27d0 - a*b/3d0 + c
        d = q**2 + 4d0*(p**3)

        !---   d>0 one real solution   ---
        call out_diffsx(d,d_x)
        if(d_x>0d0) then
            allocate(x(1),y(1),LA(2))
            d = sqrt(d)
            LA(1) = ((-1d0)*q + d)*0.5d0
            LA(2) = ((-1d0)*q +(-1d0)* d)*0.5d0
            do i=1,2
                LA(i) = fLA( LA(i) )
            end do
            y(1) = LA(1)+LA(2)

        !---   d=0 two real solutions   ---        
        elseif (d_x == 0d0) then
            allocate(x(2), y(2), LA(1))
            LA = (-1d0)*q*0.5d0
            LA(1) = fLA( LA(1) )

            y(1) = 2d0 * LA(1)
            y(2) = (-1d0)*LA(1)

        !---   d<0 three real solutions   ---
        else
            allocate(x(3), y(3))
            theta = acos( (-0.5d0) * q /  sqrt((-1d0)*p**3d0)  )
            y(1) = 2d0 * sqrt((-1d0)*p) * cos(  theta            / 3d0 )
            y(2) = 2d0 * sqrt((-1d0)*p) * cos( (theta + 2d0*pai) / 3d0 )
            y(3) = 2d0 * sqrt((-1d0)*p) * cos( (theta + 4d0*pai) / 3d0 )
        endif
    
        ns = size(y)
        do i=1,ns    
            x(i) = y(i) +(-1d0)*a/3d0
        end do
        call outxs(x,x_x)
    
        xmax=x_x(1)
        nsmax=1
        do i=1,ns-1
            if (x_x(i+1)>=xmax)then
                xmax=x_x(i+1)
                nsmax=i+1
            end if
        end do
        xmin=x_x(1)
        nsmin=1
        do i=1,ns-1
            if (x_x(i+1)<=xmax)then
                xmin=x_x(i+1)
                nsmin=i+1
            end if
        end do

        z=x(nsmax)
    
    
    contains
        type(diffs) function fLA(LA)
            type(diffs), intent(in) :: LA
            double precision ::LA0

            call out_diffsx(LA,LA0)
            if(LA0>0d0) then
                fLA = LA**(1d0/3d0)
            else
                fLA = (-1d0)*( ((-1d0)*LA)**(1d0/3d0) )
            endif
    
        end function fLA

    end subroutine 

end module
