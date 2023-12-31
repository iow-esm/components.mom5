! -*-f90-*-
! $Id: drifters_set_field.h,v 13.0 2006/03/28 21:38:37 fms Exp $

!<CONTACT EMAIL="torsten.seifert@io-warnemuende.de"> Torsten Seifert
!</CONTACT>
!<CONTACT EMAIL="klaus-ketelsen@t-online.de"> Klaus Ketelsen 
!</CONTACT>
!<DESCRIPTION>
!  IOW version 3.0 from 2014/11/28
!  code changes:  _DEBUG option
!  clip depth range to minval(z) ... maxval(z) (experimental)
!  timing & simple+fast call cld_ntrp_get_cell_values_2D by Klaus Ketelsen  
!  avoid _FLATTEN for 2D drifters
!</DESCRIPTION>
! iow ! #define _DEBUG 

subroutine drifters_set_field_XXX(self, index_field, x, y, &
#if _DIMS >= 3
     & z, &
#endif
     & data, ermesg)
  use cloud_interpolator_mod
  type(drifters_type) :: self
  ! field index must be consistent with field_names from input file
  integer, intent(in) :: index_field
  real, intent(in)    :: x(:)
  real, intent(in)    :: y(:)
#if _DIMS == 2
  real, intent(in)    :: data(:,:)
#endif
#if _DIMS == 3
  real, intent(in)    :: z(:)
  real, intent(in)    :: data(:,:,:)
#endif
  character(len=*), intent(out) :: ermesg

  integer i, j, ip, ier, ij(self%core%nd), nsizes(self%core%nd), nf
  real fvals(2**self%core%nd), ts(self%core%nd)

  call mpp_clock_begin(id_fields)

  ermesg = ''
  ! only interpolate field if RK step is complete
  if(self%rk4_step > 1) return

  ! interpolate onto new positions
  nsizes(1) = size(x)
  nsizes(2) = size(y)
#if _DIMS >= 3
  nsizes(3) = size(z)
#endif

  if(nsizes(1) /= size(data, 1) .or. nsizes(2) /= size(data, 2)) then
     ermesg = 'drifters_set_field_XXX: ERROR size mismatch between data and x or y'
     return
  end if
#if _DIMS >=3 
  if(nsizes(3) /= size(data, 3)) then
     ermesg = 'drifters_set_field_XXX: ERROR size mismatch between data and z'
     return
  endif
#endif

  if(size(self%fields, 2) < self%core%np) then
     ! resize
     deallocate(self%fields, stat=ier)
     nf = size(self%input%field_names)
     allocate(self%fields(nf, self%core%npdim))
     self%fields = -huge(1.)
  endif

  do ip = 1, self%core%np
     call cld_ntrp_locate_cell(x, self%core%positions(1,ip), i, ier)
     ij(1) = i
#ifdef _DEBUG
     if(i<1) then
        print *,'drifters_set_field_XXX::OUT OF BOUND ERROR. xmin, x, xmax=', minval(x), self%core%positions(1,ip), maxval(x)
     endif
#endif
     ts(1) = (self%core%positions(1,ip) - x(i))/(x(i+1) - x(i))

     call cld_ntrp_locate_cell(y, self%core%positions(2,ip), j, ier)
     ij(2) = j
#ifdef _DEBUG
     if(j<1) then
        print *,'drifters_set_field_XXX::OUT OF BOUND ERROR. ymin, y, ymax=', minval(y), self%core%positions(2,ip), maxval(y)
     endif
#endif
     ts(2) = (self%core%positions(2,ip) - y(j))/(y(j+1) - y(j))

#if _DIMS >= 3
!!     call cld_ntrp_locate_cell_z(z, self%core%positions(3,ip), j, ier)
!!     ij(3) = j
!! #ifdef _DEBUG
!!      if(j<1) then
!!         print *,'drifters_set_field_XXX::OUT OF BOUND ERROR. zmin, z, zmax=', minval(z), self%core%positions(3,ip), maxval(z)
!!      endif
!! #endif
!!      ts(3) = (self%core%positions(3,ip) - z(j))/(z(j+1) - z(j))
! iow clip depth range to minval(z) ... maxval(z)
     if (self%core%positions(3,ip).le.minval(z)) then
       ij(3) = 1
       ts(3) = 0.0
     elseif (self%core%positions(3,ip).ge.maxval(z)) then
       ij(3) = nsizes(3)
       ts(3) = 1.0
     else
       call cld_ntrp_locate_cell(z, self%core%positions(3,ip), j, ier)
       ij(3) = j
       ts(3) = (self%core%positions(3,ip) - z(j))/(z(j+1) - z(j))
     endif
#endif

#if _DIMS == 2
     call cld_ntrp_get_cell_values_2D(nsizes, data, ij, fvals, ier)
#else
     call cld_ntrp_get_cell_values(nsizes, _FLATTEN(data), ij, fvals, ier)
#endif
     call cld_ntrp_linear_cell_interp(fvals, ts, self%fields(index_field, ip), ier)
  enddo

  call mpp_clock_end(id_fields)

end subroutine drifters_set_field_XXX

