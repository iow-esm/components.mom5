! -*-f90-*-
! $Id: drifters_compute_k.h,v 13.0 2006/03/28 21:38:20 fms Exp $

!<CONTACT EMAIL="torsten.seifert@io-warnemuende.de"> Torsten Seifert 
!</CONTACT>
!<CONTACT EMAIL="klaus-ketelsen@t-online.de"> Klaus Ketelsen 
!</CONTACT>
!<DESCRIPTION>
!  IOW version 3.0 from 2014/11/28
!  code changes: _DEBUG option
!  locate tracer cells of drifter positions
!  clip depth range to minval(zv) ... maxval(zv) (experimental)
!  avoid _FLATTEN for 2D drifters by Klaus Ketelsen
!</DESCRIPTION>
! iow !#define _DEBUG

subroutine drifters_compute_k_XXX(self, positions, u, v, &
#if _DIMS >= 3
     & w, &
#endif
     & k, ermesg)

#ifdef _DEBUG
  use mpp_mod,only: stdlog
#endif

  use cloud_interpolator_mod
  type(drifters_type) :: self
  real, intent(in)    :: positions(:,:)
#if _DIMS == 2
  real, intent(in)    :: u(:,:)
  real, intent(in)    :: v(:,:)
#endif
#if _DIMS == 3
  real, intent(in)    :: u(:,:,:)
  real, intent(in)    :: v(:,:,:)
  real, intent(in)    :: w(:,:,:)
#endif
  real, intent(out)   :: k(:,:)
  character(len=*), intent(out) :: ermesg

  integer, parameter :: nd = _DIMS ! number of dims
  integer i, ip, np, ij(nd), ier, nsizes_u(nd), nsizes_v(nd)
#if _DIMS >= 3
  integer nsizes_w(nd)
#endif
  real fvals(2**nd), ts(nd)
  real pos(nd, self%core%np)

  call mpp_clock_begin(id_compk)

  ermesg = ''

  nsizes_u(1) = size(u, 1)
  nsizes_u(2) = size(u, 2)

  nsizes_v(1) = size(v, 1)
  nsizes_v(2) = size(v, 2)

#if _DIMS >= 3
  nsizes_u(3) = size(u, 3)
  nsizes_v(3) = size(v, 3)
  nsizes_w(1) = size(w, 1)
  nsizes_w(2) = size(w, 2)
  nsizes_w(3) = size(w, 3)
#endif

  np = self%core%np

  ! correct for periodicity
  if(self%comm%xperiodic) then
     do ip = 1, np
        pos(1,ip) = self%comm%xgmin + modulo(positions(1,ip)-self%comm%xgmin, self%comm%xgmax-self%comm%xgmin)
     enddo
  else
     pos(1,:) = positions(1,1:np)
  endif
  if(self%comm%yperiodic) then
     do ip = 1, np
        pos(2,ip) = self%comm%ygmin + modulo(positions(2,ip)-self%comm%ygmin, self%comm%ygmax-self%comm%ygmin)
     enddo
  else
     pos(2,:) = positions(2,1:np)
  endif

#if _DIMS >= 3
  pos(3,:) = positions(3,1:self%core%np)
#endif

  do ip = 1, np

     ! iterate over particles

     k(:, ip) = huge(1.)

     ! u-component...
     call cld_ntrp_locate_cell(self%xu, pos(1,ip), i, ier)

     if(i==-1) self%remove(ip) = .TRUE.
#ifdef _DEBUG
     if(i<1) then
        write(stdlog(),*) '***PE: ', _MPP_PE,' i=', i,' itt=',self%core%it,' id=',self%core%ids(ip), 'pos_xu=', pos(1,ip), ' axis min/max=', minval(self%xu), maxval(self%xu)
     endif
#endif
     i = max(1, i)
     ts(1) = (pos(1,ip) - self%xu(i))/(self%xu(i+1)-self%xu(i))
     ij(1) = i

     call cld_ntrp_locate_cell(self%yu, pos(2,ip), i, ier)
     if(i==-1) self%remove(ip) = .TRUE.
#ifdef _DEBUG
     if(i<1) then
        write(stdlog(),*) '***PE: ', _MPP_PE,' i=', i,' itt=',self%core%it,' id=',self%core%ids(ip), 'pos_yu=', pos(2,ip), ' axis min/max=', minval(self%yu), maxval(self%yu)
     endif
#endif
     i = max(1, i)
     ts(2) = (pos(2,ip) - self%yu(i))/(self%yu(i+1)-self%yu(i))
     ij(2) = i

#if _DIMS >= 3
!!      call cld_ntrp_locate_cell(self%zu, pos(3,ip), i, ier)
!!      if(i==-1) self%remove(ip) = .TRUE.
!! #ifdef _DEBUG
!!      if(i<1) then
!!         write(stdlog(),*) '***PE: ', _MPP_PE,' i=', i,' itt=',self%core%it,' id=',self%core%ids(ip), 'pos_zu=', pos(3,ip), ' axis min/max=', minval(self%zu), maxval(self%zu)
!!      endif
!! #endif
!!      i = max(1, i)
!!      ts(3) = (pos(3,ip) - self%zu(i))/(self%zu(i+1)-self%zu(i))
!!      ij(3) = i
! iow ! clip depth range to minval(zu) ... maxval(zu)
     if (pos(3,ip).le.minval(self%zu)) then
       ij(3) = 1
       ts(3) = 0.0
     elseif (pos(3,ip).ge.maxval(self%zu)) then
       ij(3) = size(self%zu)
       ts(3) = 1.0
     else
       call cld_ntrp_locate_cell(self%zu, pos(3,ip), i, ier)
       ij(3) = i
       ts(3) = (pos(3,ip) - self%zu(i))/(self%zu(i+1)-self%zu(i))
     endif
#endif

! iow ! this works for local indices
! cld_ntrp_locate_cell finds t-cell ij corresponding to self%xu(i) < pos(1) <= self%xu(i+1) etc. 
! therefore interpolation uses fvals at (i,j), (i+1,j), (i,j+1), (i+1,j+1)  
#if _DIMS == 2
     call cld_ntrp_get_cell_values_2D (nsizes_u, u, ij, fvals, ier)
#else
!kk  call to cld_ntrp_get_cell_values without reshape (_FLATTEN) has to be programmed
     call cld_ntrp_get_cell_values(nsizes_u, _FLATTEN(u), ij, fvals, ier)
#endif
     call cld_ntrp_linear_cell_interp(fvals, ts, k(1, ip), ier)
     k(1, ip) = self%dt * k(1, ip)

!kk     if(ip == 3) write(stdlog(),*) 'cell_interp ', ip, k(1, ip),sum(fvals)
!kk     if(ip == 3) write(stdlog(),'(a,8i7)') 'Aftercld_ntrp_linear_cell_interp ', ip, nsizes_u,ij,shape(u),_DIMS

#ifdef _DEBUG
     write(stdlog(),'(a,i5,1pe18.10,i5,1p2e18.10)') ' drifters_compute_xu',ip,pos(1,ip),ij(1),ts(1),k(1,ip)
     write(stdlog(),'(a,i5,1pe18.10,i5,1p2e18.10)') ' drifters_compute_yu',ip,pos(2,ip),ij(2),ts(2)
#if _DIMS >= 3
     write(stdlog(),'(a,i5,1pe18.10,i5,1p2e18.10)') ' drifters_compute_zu',ip,pos(3,ip),ij(3),ts(3)
#endif
     write(stdlog(),'(a,i5,1p4e18.10)') ' fvals_u',ip,fvals
#endif 

     ! v-component...
     call cld_ntrp_locate_cell(self%xv, pos(1,ip), i, ier)
     if(i==-1) self%remove(ip) = .TRUE.
#ifdef _DEBUG
     if(i<1) then
        write(stdlog(),*) '***PE: ', _MPP_PE,' i=', i,' itt=',self%core%it,' id=',self%core%ids(ip), 'pos_xv=', pos(1,ip), ' axis min/max=', minval(self%xv), maxval(self%xv)
     endif
#endif
     i = max(1, i)
     ts(1) = (pos(1,ip) - self%xv(i))/(self%xv(i+1)-self%xv(i))
     ij(1) = i

     call cld_ntrp_locate_cell(self%yv, pos(2,ip), i, ier)
     if(i==-1) self%remove(ip) = .TRUE.
#ifdef _DEBUG
     if(i<1) then
        write(stdlog(),*) '***PE: ', _MPP_PE,' i=', i,' itt=',self%core%it,' id=',self%core%ids(ip), 'pos_yv=', pos(2,ip), ' axis min/max=', minval(self%yv), maxval(self%yv)
     endif
#endif
     i = max(1, i)
     ts(2) = (pos(2,ip) - self%yv(i))/(self%yv(i+1)-self%yv(i))
     ij(2) = i

#if _DIMS >= 3
!!      call cld_ntrp_locate_cell(self%zv, pos(3,ip), i, ier)
!!      if(i==-1) self%remove(ip) = .TRUE.
!! #ifdef _DEBUG
!!      if(i<1) then
!!         write(stdlog(),*) '***PE: ', _MPP_PE,' i=', i,' itt=',self%core%it,' id=',self%core%ids(ip), 'pos_zv=', pos(3,ip), ' axis min/max=', minval(self%zv), maxval(self%zv)
!!      endif
!! #endif
!!      i = max(1, i)
!!      ts(3) = (pos(3,ip) - self%zv(i))/(self%zv(i+1)-self%zv(i))
!!      ij(3) = i
! iow ! clip depth range to minval(zv) ... maxval(zv)
     if (pos(3,ip).le.minval(self%zv)) then
       ij(3) = 1
       ts(3) = 0.0
     elseif (pos(3,ip).ge.maxval(self%zv)) then
       ij(3) = size(self%zv)
       ts(3) = 1.0
     else
       call cld_ntrp_locate_cell(self%zv, pos(3,ip), i, ier)
       ij(3) = i
       ts(3) = (pos(3,ip) - self%zv(i))/(self%zv(i+1)-self%zv(i))
     endif
#endif

#if _DIMS == 2
     call cld_ntrp_get_cell_values_2D (nsizes_v, v, ij, fvals, ier)
#else
     call cld_ntrp_get_cell_values(nsizes_v, _FLATTEN(v), ij, fvals, ier)
#endif
     call cld_ntrp_linear_cell_interp(fvals, ts, k(2, ip), ier)
     k(2, ip) = self%dt * k(2, ip)

#ifdef _DEBUG
     write(stdlog(),'(a,i5,1pe18.10,i5,1p2e18.10)') ' drifters_compute_xv',ip,pos(1,ip),ij(1),ts(1)
     write(stdlog(),'(a,i5,1pe18.10,i5,1p2e18.10)') ' drifters_compute_yv',ip,pos(2,ip),ij(2),ts(2),k(2,ip)
#if _DIMS >= 3
     write(stdlog(),'(a,i5,1pe18.10,i5,1p2e18.10)') ' drifters_compute_zv',ip,pos(3,ip),ij(3),ts(3)
#endif
     write(stdlog(),'(a,i5,1p4e18.10)') ' fvals_v',ip,fvals
#endif

#if _DIMS >= 3
     ! w-component...
     call cld_ntrp_locate_cell(self%xw, pos(1,ip), i, ier)
     if(i==-1) self%remove(ip) = .TRUE.
#ifdef _DEBUG
     if(i<1) then
        write(stdlog(),*) '***PE: ', _MPP_PE,' i=', i,' itt=',self%core%it,' id=',self%core%ids(ip), 'pos_xw=', pos(1,ip), ' axis min/max=', minval(self%xw), maxval(self%xw)
     endif
#endif
     i = max(1, i)
     ts(1) = (pos(1,ip) - self%xw(i))/(self%xw(i+1)-self%xw(i))
     ij(1) = i

     call cld_ntrp_locate_cell(self%yw, pos(2,ip), i, ier)
     if(i==-1) self%remove(ip) = .TRUE.
#ifdef _DEBUG
     if(i<1) then
        write(stdlog(),*) '***PE: ', _MPP_PE,' i=', i,' itt=',self%core%it,' id=',self%core%ids(ip), 'pos_yw=', pos(2,ip), ' axis min/max=', minval(self%yw), maxval(self%yw)
     endif
#endif
     i = max(1, i)
     ts(2) = (pos(2,ip) - self%yw(i))/(self%yw(i+1)-self%yw(i))
     ij(2) = i
!!      call cld_ntrp_locate_cell(self%zw, pos(3,ip), i, ier)
!!      if(i==-1) self%remove(ip) = .TRUE.
!! #ifdef _DEBUG
!!      if(i<1) then
!!         write(stdlog(),*) '***PE: ', _MPP_PE,' i=', i,' itt=',self%core%it,' id=',self%core%ids(ip), 'pos_zw=', pos(3,ip), ' axis min/max=', minval(self%zw), maxval(self%zw)
!!      endif
!! #endif
!!      i = max(1, i)
!!      ts(3) = (pos(3,ip) - self%zw(i))/(self%zw(i+1)-self%zw(i))
!!      ij(3) = i
! iow ! clip depth range to minval(zv) ... maxval(zv)
     if (pos(3,ip).le.minval(self%zw)) then
       ij(3) = 1
       ts(3) = 0.0
     elseif (pos(3,ip).ge.maxval(self%zw)) then
       ij(3) = size(self%zw)
       ts(3) = 1.0
     else
       call cld_ntrp_locate_cell(self%zw, pos(3,ip), i, ier)
       ij(3) = i
       ts(3) = (pos(3,ip) - self%zw(i))/(self%zw(i+1)-self%zw(i))
     endif

     call cld_ntrp_get_cell_values(nsizes_w, _FLATTEN(w), ij, fvals, ier)
     call cld_ntrp_linear_cell_interp(fvals, ts, k(3, ip), ier)
     k(3, ip) = self%dt * k(3, ip)

#ifdef _DEBUG
     write(stdlog(),'(a,1pe18.10,i5,1p2e18.10)') ' drifters_compute_xw',pos(1,ip),ij(1),ts(1)
     write(stdlog(),'(a,1pe18.10,i5,1p2e18.10)') ' drifters_compute_yw',pos(2,ip),ij(2),ts(2)
     write(stdlog(),'(a,1pe18.10,i5,1p2e18.10)') ' drifters_compute_zw',pos(3,ip),ij(3),ts(3),k(3,ip)
     write(stdlog(),'(a,i5,1p4e18.10)') ' fvals_w',ip,fvals
#endif

#endif

  enddo

  call mpp_clock_end(id_compk)

end subroutine drifters_compute_k_XXX

