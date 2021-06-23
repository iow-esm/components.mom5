! -*-f90-*-
! $Id: drifters_push.h,v 13.0 2006/03/28 21:38:35 fms Exp $

!<CONTACT EMAIL="torsten.seifert@io-warnemuende.de"> Torsten Seifert 
!</CONTACT>
!<DESCRIPTION>
!  IOW version 3.0 from 2014/11/28
!  code changes: _DEBUG option
!  initial buffer size max_add_remove=1000 (will be increased if necessary)
!  reset Runge-Kutta contributions to zero
!  clip depth range to z >=0 for 3D (experimental)
!</DESCRIPTION>
! iow ! #define _DEBUG
subroutine drifters_push_XXX(self, u, v, &
#if _DIMS >= 3
  & w, &
#endif
  & ermesg)
  use mpp_mod,only: stdlog !iow!
  type(drifters_type) :: self
#if _DIMS == 2
  real, intent(in)    :: u(:,:)
  real, intent(in)    :: v(:,:)
#endif
#if _DIMS == 3
  real, intent(in)    :: u(:,:,:)
  real, intent(in)    :: v(:,:,:)
  real, intent(in)    :: w(:,:,:)
#endif
  character(len=*), intent(out) :: ermesg

  integer i, np, nf, max_add_remove
  
  ermesg = ''
  np = self%core%np
  nf = size(self%input%field_names)

  if(self%core%nd /= _DIMS) then
     ermesg = 'drifters_push: wrong number of dimensions'
     return
  endif

  select case(self%rk4_step)

  case (0)

     ! only invoked at the first step

     call drifters_reset_rk4(self, ermesg=ermesg)

     call drifters_compute_k(self, self%core%positions, &
          & u=u, v=v, &
#if _DIMS >= 3
          & w=w, &
#endif
          & k=self%rk4_k1, ermesg=ermesg)

     self%temp_pos(:,1:np) = self%core%positions(:,1:np) + 0.5*self%rk4_k1(:,1:np)

     self%rk4_step = 1
     self%rk4_completed = .FALSE.

  case (1)

     call drifters_compute_k(self, self%temp_pos, &
          & u=u, v=v, &
#if _DIMS >= 3
          & w=w, &
#endif
          & k=self%rk4_k2, ermesg=ermesg)

     self%temp_pos(:,1:np) = self%core%positions(:,1:np) + 0.5*self%rk4_k2(:,1:np)

     call drifters_compute_k(self, self%temp_pos, &
          & u=u, v=v, &
#if _DIMS >= 3
          & w=w, &
#endif
          & k=self%rk4_k3, ermesg=ermesg)

     self%temp_pos(:,1:np) = self%core%positions(:,1:np) +     self%rk4_k3(:,1:np)

     self%rk4_step = 2
     self%rk4_completed = .FALSE.

  case (2)

     call drifters_compute_k(self, self%temp_pos, &
          & u=u, v=v, &
#if _DIMS >= 3
          & w=w, &
#endif
          & k=self%rk4_k4, ermesg=ermesg)

     ! This completes the RK4 steps
     do i = 1, np
        self%temp_pos(:,i) = self%core%positions(:,i) + &
             & (self%rk4_k1(:,i) + 2*self%rk4_k2(:,i) + 2*self%rk4_k3(:,i) + self%rk4_k4(:,i))/6.

#if _DIMS >= 3
#ifdef _DEBUG
        if (self%temp_pos(3,i).lt.0.0) write(*,*) 'np,z<0',np,self%temp_pos(3,i)
#endif
        self%temp_pos(3,i) = max(self%temp_pos(3,i),0e0) ! iow ! keep z>=0
#endif
     enddo

     ! correct for periodic domain, if necessary
     call drifters_modulo(self, self%temp_pos, ermesg=ermesg)

     ! estimate of max number of particles to add/remove
     ! this may need to be adjusted over time... 
     ! [_MPP_NPES is a macro for mpp_npes()]
     if(self%comm%pe_end < 0) self%comm%pe_end = _MPP_NPES - 1
!!      max_add_remove = &
!!           max(  10, &
!! 	  &     int(0.2*self%core%npdim/ &
!!           &        (self%comm%pe_end-self%comm%pe_beg+1)), &
!! 	  &     int(np * 2*(self%nx+self%ny)/real(self%nx*self%ny)) &
!! 	  &  ) 
!! write(*,'(a,10i6)') 'pe,max_add_remove,npdim,pe_end,pe_beg,np,nx,ny,int1,int2', &
!! self%comm%pe,max_add_remove,self%core%npdim,self%comm%pe_end,self%comm%pe_beg,np,self%nx,self%ny, &
!! int(0.2*self%core%npdim/(self%comm%pe_end-self%comm%pe_beg+1)), &
!! int(np * 2*(self%nx+self%ny)/real(self%nx*self%ny))

     ! copy/move drifters across domain boundaries, if necessary
! iow ! initial buffer size for data_send and indices_to_remove
     max_add_remove = 1000
     call drifters_comm_update(self%comm, self%core, self%temp_pos(:,1:np), &
          & remove=self%remove(1:np), max_add_remove=max_add_remove)

     np = self%core%np

     ! update time
     self%time = self%time + self%dt

     if(self%core%npdim < self%core%np) self%core%npdim = int(1.2 * self%core%np)

     ! resize local arrays if necessary
     call drifters_reset_rk4(self, ermesg=ermesg)

! iow ! reset rk4_kn=0
     self%rk4_k1(:,:) = 0e0
     self%rk4_k2(:,:) = 0e0
     self%rk4_k3(:,:) = 0e0
     self%rk4_k4(:,:) = 0e0

     self%remove = .FALSE. ! by definition no drifters outside compute domain

     ! prepare for next step...
     call drifters_compute_k(self, self%core%positions, &
          & u=u, v=v, &
#if _DIMS >= 3
          & w=w, &
#endif
          & k=self%rk4_k1, ermesg=ermesg)

     self%temp_pos(:, 1:np) = self%core%positions(:, 1:np) + 0.5*self%rk4_k1(:, 1:np)

     self%rk4_step = 1
     self%rk4_completed = .TRUE.

  case default
     ermesg = 'drifters_push: invalid rk4_step'

  end select

end subroutine drifters_push_XXX
