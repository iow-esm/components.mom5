#IFDEF COUP_OAS
SUBROUTINE oas_send_field(Ice,Time_start,Timet)

!!---------------------------------------------------------------------
!!              ***  ROUTINE send_fld  ***
!!
!! ** Purpose : Prepare and send coupling fields to OASIS
!!      
!!----------------------------------------------------------------------

USE oas_vardef

USE mpp_mod,         only: mpp_pe, mpp_root_pe
USE ice_grid_mod
USE ice_type_mod,    only: ice_data_type    
USE time_manager_mod, only : get_time,get_time,operator(-),time_type
USE netcdf

IMPLICIT NONE
!
! local parameters, variables and arrays
!
INTEGER :: &
  isec,            &
  kinfo,           &
  il_var_id(16),   &  
  dimids(2),       &
  k,j,i,jn,sc,dy,  &
  dt,tcks,sec,day
       

REAL (KIND=8)      :: &
  ztmp1  (iec,jec),    &
  ztmp2  (iec,jec)

TYPE (ice_data_type),      INTENT(IN) :: Ice
type(time_type),       intent(in) :: Time_start,Timet ! for coupling sandra

! ! call get_time(Ice%Time-Ice%Time_init,sc,dy)
! ! isec = sc+(dy*86400) 

! ! call get_time(Ice%Time_step_slow,scs,dys)
! !  call get_time(Ice%Time-Time_start, sc, dy)
!   call get_time(Timet-Time_start, sc, dy)
!   isec = (864e2*dy+sc)
! !   write(*,*) "Time",Ice%Time,Time_start,Ice%Time_init,isec


! !----------------------------------------------------------------------------
! ! handling of sent fields
! !----------------------------------------------------------------------------
!   jn = 1.0 
!   ztmp1(:,:)=0.0
!   ztmp2(:,:)=0.0

! !write(*,*)" oasis send laction: " , ssnd(jn)%laction, isec 
!  !! jn = jn + 1
! !  IF( ssnd(jn)%laction ) THEN                                         ! sea surface temperature        
! !    write(*,*) "max t_surf mom: ",MAXVAL(Ice%t_surf(isc:iec,jsc:jec,1)), isec
! !    write(*,*) "min t_surf mom: ",MINVAL(Ice%t_surf(isc:iec,jsc:jec,1)), isec

! !    ztmp1(:,:)=0.
! !    do k=1,1 !,km
! !       do j = jsc, jec
! !          do i = isc, iec
! !             ztmp1(i,j) = ztmp1(i,j) + Ice%t_surf(i,j,k)
! !          enddo
! !       enddo
! !    enddo

!    ztmp1(isc:iec,jsc:jec) =Ice%t_surf(isc:iec,jsc:jec,1)             

!    CALL oas_send (jn, isec,ztmp1, kinfo )
! !  ENDIF

!   jn = jn + 1
!  ! IF( ssnd(jn)%laction ) THEN                                         ! sea ice area fraction
!     ztmp2(:,:)=0.
!     do k=2,km
!        do j = jsc, jec
!           do i = isc, iec
!              ztmp2(i,j) = ztmp2(i,j) + Ice%part_size(i,j,k)
!           enddo
!        enddo
!     enddo
!     CALL oas_send (jn, isec, ztmp2, kinfo )
!  ! ENDIF

END SUBROUTINE oas_send_field
#ENDIF
