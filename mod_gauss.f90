module mod_gauss
    implicit none
    contains
  
  subroutine pvgauss(n,a,b)
      implicit none 
      ! Dummy arguments
      integer,intent(in)::n
      real(8),intent(inout)::a(:,:),b(:)
      
      ! Local variables
      integer:: i,j,k,jp(n),imax,jmax,iik
      real(8):: akk,aik,val,pivot
      
      ! Set the initial order of a solution vector
      do i=1,n
        jp(i)=i
      enddo
      
      do k=1,n
        ! Set a(k,k) as pivot  
        pivot=dabs(a(k,k)); imax=k; jmax=k
        
        ! Find the largest element from the lower right square, a(k-n,k-n)
        do i=k,n
          do j=k,n
            if(dabs(a(i,j))>pivot)then
              pivot=dabs(a(i,j))
              imax=i; jmax=j
            endif
          enddo
        enddo
        
        ! Re-ordering due to the interchange of k-th and jmax-th laws
        if (jmax/=k) then
          iik=jp(k); jp(k)=jp(jmax); jp(jmax)=iik
        endif
        
        ! Switch k-th and jmax-th laws
        do i=1,n
          pivot=a(i,k); a(i,k)=a(i,jmax); a(i,jmax)=pivot
        enddo
        
        ! Switch k-th and imax-th lines
        do j=k,n
          pivot=a(k,j); a(k,j)=a(imax,j); a(imax,j)=pivot
        enddo
        
        ! Switch k-th and imax-th elements in the right-hand side vector
        pivot=b(k); b(k)=b(imax); b(imax)=pivot
        
        ! Forward elimination
        akk=a(k,k)
        do j=k,n
          a(k,j)=a(k,j)/akk
        enddo
        b(k)=b(k)/akk
        
        do i=k+1,n
          aik=a(i,k)
          do j=k,n
            a(i,j)=a(i,j)-aik*a(k,j)
          enddo
          b(i)=b(i)-aik*b(k)
        enddo
      enddo
      
      ! Backward substitution
      do k=n,1,-1
        val=0.0d0
        do j=k+1,n
          val=val+a(k,j)*b(j)
        enddo
        b(k)=b(k)-val
      enddo
      
      ! Re-ordering
      b(jp)=b
      
      return
  end subroutine pvgauss
  end module