#IFDEF OASIS_IOW_ESM
SUBROUTINE oas_exchange_fields(Ice,Ice_boundary,Ice_Ocean_Boundary,Time_start,Timet)

!!---------------------------------------------------------------------
!!              ***  ROUTINE send_fld  ***
!!
!! ** Purpose : Prepare and send coupling fields to OASIS
!!      
!!----------------------------------------------------------------------

USE oas_vardef

USE mpp_mod,         only: mpp_pe, mpp_root_pe
USE ice_grid_mod
USE ice_model_mod,    only: ice_data_type,atmos_ice_boundary_type
USE ocean_model_mod, only: ice_ocean_boundary_type 
USE time_manager_mod, only : get_time,get_time,operator(-),time_type

IMPLICIT NONE
!
! local parameters, variables and arrays
!
INTEGER :: &
  isec,            &
  kinfo,           &
  jn,              &
  sc,dy,dummy
       

REAL (KIND=8)      :: &
  ztmp1  (iec,jec)

TYPE (ice_data_type),      INTENT(IN) :: Ice
TYPE (atmos_ice_boundary_type), INTENT(INOUT) :: Ice_boundary
type(ice_ocean_boundary_type), intent(inout) :: Ice_Ocean_Boundary
type(time_type),           intent(in) :: Time_start,Timet ! for coupling sandra
integer, parameter :: dt_cpld=600

  call get_time(Timet-Time_start, sc, dy)
  isec = (864e2*dy+sc)

  write(*,*) 'inside oas_exchange_fields at isec=',isec
  if (mod(isec,dt_cpld)==0) then
    write(*,*) 'calling MPI_BARRIER'
    call MPI_BARRIER(MPI_COMM_WORLD,dummy)
    write(*,*) 'passed MPI_BARRIER'
  endif
!----------------------------------------------------------------------------
! STEP 1: Send area fraction of each ice class and surface temperature of each ice class 
!----------------------------------------------------------------------------
  write(*,*) 'Sending MSFARE01..MSFARE06 at isec=',isec
  DO jn=1,6  
    ztmp1(isc:iec,jsc:jec) =Ice%part_size(isc:iec,jsc:jec,jn)             
    CALL oas_send (jn, isec,ztmp1, kinfo )   ! send MSFARE01..MSFARE06
  ENDDO
  write(*,*) 'Sending MSTSUR01..MSTSUR06 at isec=',isec
  DO jn=1,6
    ztmp1(isc:iec,jsc:jec) =Ice%t_surf(isc:iec,jsc:jec,jn)             
    CALL oas_send (jn+6, isec,ztmp1, kinfo )   ! send MSTSUR01..MSTSUR06
  ENDDO
  write(*,*) 'oas_exchange_fields finished.'
END SUBROUTINE oas_exchange_fields
#ENDIF
