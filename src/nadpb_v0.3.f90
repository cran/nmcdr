!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                             !! 
!!    nonparametric multiple change-point detect procedure     !! 
!!    author: Changliang Zou                                   !! 
!!    license:  GPL-2                                          !! 
!!                                                             !! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine NMCDP(n,y,steps,ncp,cpp,ke,cpe,w,tmin)
!DEC$ ATTRIBUTES DLLEXPORT,C,REFERENCE,ALIAS:'nmcdp_' :: NMCDP
!! n: number of data
!! y: data
!! steps: upper bound of the number of change-points
!! ncp: number of candidate change-points returned by isp
!! cpp: candiate change-points 
!! ke: true number of change-points
!! cpe:  true change-points
!! w: store the middle value of DP algorithm
 integer n,steps,m,cpp(ncp),ncp,cpp1(0:ncp+1),h(2:steps+1,n)
!! cpp1 is simply the 'extended' cpp
!! m is just a middle variate.
!! h is used for tracing-back procedure in DP algorithm. 
 integer cpe(steps+1),cpee(steps+1),ke
 real*8 y(n),w(0:(n-1),n),x(0:n),fa(steps+1,n),bic(n),tmin
!! x : 'extended' data
!! xx : sorted data y
!! fa: store value for DP algorithm.
!! tmin: store the minimum BIC.
 real*8 tempmin
 integer r1(n),r(n),rr(n),ipe(n),r2(n),ni
!! r1: sort-index
!! r: sorted index.
!! rr: soted index for middle variable.
!! r2: only used for initialization.
!! ni: middle variate. 
!! ipe: middle variable.

 tmin = 10000000.0d0
 x(0) = 0.0d0
 x(1:n) = y
 w = 0.0d0
 cpp1(0) = 0
 cpp1(1:ncp) = cpp(1:ncp)-1
 cpp1(ncp+1) = n
 
 h = 0
 do i = 1,n
   r1(i) = i
   r2(i) = i
 enddo
 call Dquicksort(n,x(1:n),r1)
 do i = 1,n
   r(r1(i)) = i
 enddo
 do m = 1,ncp+1
   do i = 0,m-1
      w(i,m) = 0.0d0
	    ni = cpp1(m)-cpp1(i)
	    tempmin = 100.0d0
	    if(ni/=0) then
	       ipe(1:ni)=r2(1:ni)
	       rr(cpp1(i)+1:cpp1(m)) = r(cpp1(i)+1:cpp1(m))
	       call Iquicksort(ni,rr(cpp1(i)+1:cpp1(m)),ipe(1:ni))	
         l=0
	       do k = 1,r(cpp1(i)+ipe(1))-1
	          l=l+1
	          f_ik=0.0
	          if(f_ik==0.0) then
	            aa=0.0d0
              else
   		        aa=-(f_ik*log(f_ik)+(1.0-f_ik)*log(1.0-f_ik))/((k-0.5d0)*(n-k+0.5d0))
            endif
	          w(i,m)=w(i,m)+aa
	       enddo
	       do j=1,ni-1
            l=l+1
	          f_ik=(j-0.5d0)/dble(ni)
	          aa=-(f_ik*log(f_ik)+(1.0-f_ik)*log(1.0-f_ik))/((r(cpp1(i)+ipe(j))-0.5d0)*( &
&               n-r(cpp1(i)+ipe(j))+0.5d0))
	          w(i,m)=w(i,m)+aa
            do k=r(cpp1(i)+ipe(j))+1,r(cpp1(i)+ipe(j+1))-1
               l=l+1 
	             f_ik=(j)/dble(ni)
               aa=-(f_ik*log(f_ik)+(1.0-f_ik)*log(1.0-f_ik))/((k-0.5d0)*(n-k+0.5d0))
	             w(i,m)=w(i,m)+aa
            enddo
         enddo
         l=l+1
	       f_ik=(ni-0.5d0)/dble(ni)
	       aa=-(f_ik*log(f_ik)+(1.0-f_ik)*log(1.0-f_ik))/((r(cpp1(i)+ipe(ni))-0.5d0)*(&
&            n-r(cpp1(i)+ipe(ni))+0.5d0))
	       w(i,m)=w(i,m)+aa
         do k=r(cpp1(i)+ipe(ni))+1,n
            l=l+1
	          f_ik=1.0d0
            if(f_ik==1.0) then
	             aa=0.0d0
            else
   		         aa=-(f_ik*log(f_ik)+(1.0-f_ik)*log(1.0-f_ik))/((k-0.5d0)*(n-k+0.5d0))
	          endif
		           w(i,m)=w(i,m)+aa
	    enddo
      w(i,m)=1.0*ni*w(i,m)
	 endif

   enddo
 enddo
 
 fa(1,:) = w(0,:)
 do ll = 1, steps-1
   do i = 2, ll+1
     do m=1,ncp+1
     tempmin = 10000000.0d0
       do j = 1, m-1
         aa = fa(i-1,j)+w(j,m)
	     if (tempmin>aa) then 
	       tempmin = aa
	       h(i,m) = j
	     endif
       enddo
       fa(i,m) = tempmin
     enddo	 
   enddo
   cpee(ll+1) = ncp+1
   bic(ll) = fa(ll+1,ncp+1)+0.5*(log(1.0*n))**2.0*ll/n
   if(bic(ll) < tmin) then
    tmin = bic(ll)
    do i = ll, 1, -1
      cpee(i) = h(i+1,cpee(i+1))
    enddo
    cpe(2:ll+1) = cpp1(cpee(1:ll))+1
	ke = ll
   endif
  enddo
 RETURN
end subroutine NMCDP



RECURSIVE SUBROUTINE Iquicksort(n,list, order) 
! intialize order to 1:n before calling.
implicit none

INTEGER, DIMENSION (n), INTENT(INOUT)  :: list,order

! Local variable
INTEGER :: i,n

CALL quick_sort_1(1, SIZE(list))

CONTAINS

RECURSIVE SUBROUTINE quick_sort_1(left_end, right_end) 

implicit none
INTEGER, INTENT(IN) :: left_end, right_end

!     Local variables
INTEGER             :: i, j, itemp,temp, reference
INTEGER, PARAMETER  :: max_simple_sort_size = 6

IF (right_end < left_end + max_simple_sort_size) THEN
  ! Use interchange sort for small lists
  CALL interchange_sort(left_end, right_end)

ELSE
  ! Use partition ("quick") sort
  reference = list((left_end + right_end)/2)
  i = left_end - 1; j = right_end + 1

  DO
    ! Scan list from left end until element >= reference is found
    DO
      i = i + 1
      IF (list(i) >= reference) EXIT
    END DO
    ! Scan list from right end until element <= reference is found
    DO
      j = j - 1
      IF (list(j) <= reference) EXIT
    END DO


    IF (i < j) THEN
      ! Swap two out-of-order elements
      temp = list(i); list(i) = list(j); list(j) = temp
      itemp = order(i); order(i) = order(j); order(j) = itemp
    ELSE IF (i == j) THEN
      i = i + 1
      EXIT
    ELSE
      EXIT
    END IF
  END DO

  IF (left_end < j) CALL quick_sort_1(left_end, j)
  IF (i < right_end) CALL quick_sort_1(i, right_end)
END IF

END SUBROUTINE quick_sort_1

SUBROUTINE interchange_sort(left_end, right_end) 

implicit none
INTEGER, INTENT(IN) :: left_end, right_end

!     Local variables
INTEGER             :: i, j, itemp, temp

DO i = left_end, right_end - 1
  DO j = i+1, right_end
    IF (list(i) > list(j)) THEN
      temp = list(i); list(i) = list(j); list(j) = temp
      itemp = order(i); order(i) = order(j); order(j) = itemp
    END IF
  END DO
END DO

END SUBROUTINE interchange_sort

END SUBROUTINE Iquicksort


RECURSIVE SUBROUTINE Dquicksort(n,list, order) 

implicit none

REAL*8, DIMENSION (n), INTENT(INOUT)  :: list
INTEGER, DIMENSION (n), INTENT(INOUT)  :: order

! Local variable
INTEGER :: i,n

CALL quick_sort_1(1, SIZE(list))

CONTAINS

RECURSIVE SUBROUTINE quick_sort_1(left_end, right_end) 

implicit none
INTEGER, INTENT(IN) :: left_end, right_end

!     Local variables
INTEGER             :: i, j, itemp
REAL*8                :: reference, temp
INTEGER, PARAMETER  :: max_simple_sort_size = 6

IF (right_end < left_end + max_simple_sort_size) THEN
  ! Use interchange sort for small lists
  CALL interchange_sort(left_end, right_end)

ELSE
  ! Use partition ("quick") sort
  reference = list((left_end + right_end)/2)
  i = left_end - 1; j = right_end + 1

  DO
    ! Scan list from left end until element >= reference is found
    DO
      i = i + 1
      IF (list(i) >= reference) EXIT
    END DO
    ! Scan list from right end until element <= reference is found
    DO
      j = j - 1
      IF (list(j) <= reference) EXIT
    END DO


    IF (i < j) THEN
      ! Swap two out-of-order elements
      temp = list(i); list(i) = list(j); list(j) = temp
      itemp = order(i); order(i) = order(j); order(j) = itemp
    ELSE IF (i == j) THEN
      i = i + 1
      EXIT
    ELSE
      EXIT
    END IF
  END DO

  IF (left_end < j) CALL quick_sort_1(left_end, j)
  IF (i < right_end) CALL quick_sort_1(i, right_end)
END IF

END SUBROUTINE quick_sort_1



SUBROUTINE interchange_sort(left_end, right_end) 

implicit none
INTEGER, INTENT(IN) :: left_end, right_end

!     Local variables
INTEGER             :: i, j, itemp
REAL                :: temp

DO i = left_end, right_end - 1
  DO j = i+1, right_end
    IF (list(i) > list(j)) THEN
      temp = list(i); list(i) = list(j); list(j) = temp
      itemp = order(i); order(i) = order(j); order(j) = itemp
    END IF
  END DO
END DO

END SUBROUTINE interchange_sort

END SUBROUTINE Dquicksort