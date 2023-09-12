module mod_cubic_eq
    contains
  subroutine cubic_eq(A,B,C,x)
      implicit none
      real(8),intent(in)::A,B,C
      real(8),intent(out)::x
      real(8)::f,df,dx,x1,error,epsi
      integer::k,iter
      
      iter=1000
      epsi=0.0001d0
  
  
      x=1.0d0
  
      do k=1,iter
        f=x**3.0d0+A*x**2.0d0+B*x+C
        df=3.0d0*x**2.0d0+2.0d0*A*x+B
        dx=-f/df
        x1=x+dx
        error=dabs(x1-x)
        if(error<epsi)exit
        x=x1
      end do
  
  end subroutine cubic_eq
  end module mod_cubic_eq
  
  