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

USE netcdf

IMPLICIT NONE
!
! local parameters, variables and arrays
!
INTEGER         :: &
  maskt(iec-isc+1,jec-jsc+1),  & ! mask array
  maskc(iec-isc+1,jec-jsc+1) ! mask array

INTEGER :: &
  isec,            &
  kinfo,           &
  jn,              &
  ncfileid, ncvarid,  & ! NetCDF IDs
  sc,dy,k,dummy,   &
  nrcvinfo(nfld_rcv_tot)  ! OASIS info argument
       

REAL (KIND=8)      :: &
  ztmp1 (iec,jec)

TYPE (ice_data_type),      INTENT(IN) :: Ice
TYPE (atmos_ice_boundary_type), INTENT(INOUT) :: Ice_boundary
type(ice_ocean_boundary_type), intent(inout) :: Ice_Ocean_Boundary
type(time_type),           intent(in) :: Time_start,Timet ! for coupling sandra
integer, parameter :: dt_cpld=600

  ztmp1 (:,:) = 0
  nrcvinfo (:) = OASIS_idle

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
  write(*,*) 'Sending MSALBE01..MSALBE06 at isec=',isec
  DO jn=1,6
    ztmp1(isc:iec,jsc:jec) =Ice%albedo_vis_dir(isc:iec,jsc:jec,jn)             
    CALL oas_send (jn+12, isec,ztmp1, kinfo )   ! send MSALBE01..MSALBE06
  ENDDO
  

!----------------------------------------------------------------------------
! STEP 2: Receive all fluxes from flux calculator 
!----------------------------------------------------------------------------
  
  ! 2.1 Receive fields via coupler from flux calculator

  DO jn = 1, nfld_rcv_tot
  ! IF( srcv(jn)%laction ) THEN
      CALL oas_recieve( jn, isec, ztmp1(:,:), nrcvinfo(jn) )
      !write(*,*) " nrcvinfo(jn) ", nrcvinfo(jn), jn, OASIS_Rcv, isec
      IF( nrcvinfo(jn) == OASIS_Rcv ) frcv(:,:,jn)=ztmp1(:,:)
  ! ENDIF
  ENDDO
  !write(*,*)" oasis receive laction: " , srcv(jn)%laction, isec
    
  !IF (ltime) CALL get_timings (i_cpl_get, ntstep, dt, izerror)
    
  
  ! 2.2 Obtain coupling mask for fileds

  istatus=nf90_open('masks.nc', NF90_NOWRITE, ncfileid)
  istatus=nf90_inq_varid(ncfileid, 'tmom.msk' , ncvarid)
  istatus=nf90_get_var(ncfileid, ncvarid, maskt, &
          (/ isc, jsc /), (/ iec-isc+1,jec-jsc+1 /))
  istatus=nf90_inq_varid(ncfileid, 'cmom.msk' , ncvarid)
  istatus=nf90_get_var(ncfileid, ncvarid, maskc, &
          (/ isc, jsc /), (/ iec-isc+1,jec-jsc+1 /))

  istatus=nf90_close(ncfileid)
   
  
  ! 2.3 Store the received fluxes

  jn =  1
  !IF( nrcvinfo(jn) == OASIS_Rcv  ) THEN                                                              ! rain 1)
    do k = 1, size(Ice_boundary%lprec,3)
      WHERE (maskt == 0) Ice_boundary%lprec(isc:iec,jsc:jec,k) = -frcv(isc:iec,jsc:jec,jn)
    enddo
  !ENDIF
  
  jn = 2
  !IF( nrcvinfo(jn) == OASIS_Rcv  ) THEN                                                              ! evaporation 2..7)
    do k = 1, size(Ice_boundary%q_flux,3)
      WHERE (maskt == 0) Ice_boundary%q_flux(isc:iec,jsc:jec,k) =  frcv(isc:iec,jsc:jec,jn)
      jn = jn + 1 ! count up for different surface types
    enddo
  !ENDIF 
  
  jn = 8
  !IF( nrcvinfo(jn) == OASIS_Rcv  ) THEN                                                              ! snow 8)
    do k = 1, size(Ice_boundary%fprec,3)
      WHERE (maskt == 0) Ice_boundary%fprec(isc:iec,jsc:jec,k) = -frcv(isc:iec,jsc:jec,jn)
    enddo
  !ENDIF

  jn = 9
  !IF( nrcvinfo(jn) == OASIS_Rcv  ) THEN                                                              ! sea level pressure 9)
    do k = 1, size(Ice_boundary%p,3)
      WHERE (maskt == 0) Ice_boundary%p(isc:iec,jsc:jec,k) =  frcv(isc:iec,jsc:jec,jn)
    enddo
  !ENDIF

  jn = 10
  !IF( nrcvinfo(jn) == OASIS_Rcv  ) THEN                                                              ! u velocity 10)
    WHERE (maskc == 0) Ice_Ocean_Boundary%u_wind(isc:iec,jsc:jec) =  frcv(isc:iec,jsc:jec,jn)
  !ENDIF

  jn = 11
  !IF( nrcvinfo(jn) == OASIS_Rcv  ) THEN                                                              ! v velocity 11)
    WHERE (maskc == 0) Ice_Ocean_Boundary%v_wind(isc:iec,jsc:jec) = frcv(isc:iec,jsc:jec,jn)
  !ENDIF

  jn = 12
  !IF( nrcvinfo(jn) == OASIS_Rcv  ) THEN                                                              ! u wind stress 12..17)
    do k = 1, size(Ice_boundary%u_flux,3)
      WHERE (maskc == 0) Ice_boundary%u_flux(isc:iec,jsc:jec,k) =  frcv(isc:iec,jsc:jec,jn)
      jn = jn + 1
    enddo
  !ENDIF

  jn = 18
  !IF( nrcvinfo(jn) == OASIS_Rcv  ) THEN                                                              ! v wind stress 18..23)
    do k = 1, size(Ice_boundary%v_flux,3)
      WHERE (maskc == 0) Ice_boundary%v_flux(isc:iec,jsc:jec,k) =  frcv(isc:iec,jsc:jec,jn)
      jn = jn + 1
    enddo
  !ENDIF

  jn = 24
  !IF( nrcvinfo(jn) == OASIS_Rcv  ) THEN                                                              ! longwave radiation upward 24..29)
    do k = 1, size(Ice_boundary%lw_flux,3)
      WHERE (maskt == 0) Ice_boundary%lw_flux(isc:iec,jsc:jec,k) =  frcv(isc:iec,jsc:jec,jn)
      jn = jn + 1
    enddo
  !ENDIF

  jn = 30
  !IF( nrcvinfo(jn) == OASIS_Rcv  ) THEN                                                              ! longwave radiation downward 30)
    do k = 1, size(Ice_boundary%lw_flux,3)
      WHERE (maskt == 0) Ice_boundary%lw_flux(isc:iec,jsc:jec,k) = -Ice_boundary%lw_flux(isc:iec,jsc:jec,k) - frcv(isc:iec,jsc:jec,jn)
    enddo
  !ENDIF

  jn = 31
  !IF( nrcvinfo(jn) == OASIS_Rcv ) THEN                                                               ! shortwave radiation direct 31..36
    do k = 1, size(Ice_boundary%sw_flux_vis_dir,3)
      WHERE (maskt == 0) Ice_boundary%sw_flux_vis_dir(isc:iec,jsc:jec,k) =  -frcv(isc:iec,jsc:jec,jn)
      jn = jn + 1
    enddo
  !ENDIF

  jn = 37
  !IF( nrcvinfo(jn) == OASIS_Rcv ) THEN                                                              ! shortware radiation diffusive 37)
    do k = 1, size(Ice_boundary%sw_flux_vis_dif,3)
      WHERE (maskt == 0) Ice_boundary%sw_flux_vis_dif(isc:iec,jsc:jec,k) =  -frcv(isc:iec,jsc:jec,jn)
    enddo
  !ENDIF

  jn = 38
  !IF( nrcvinfo(jn) == OASIS_Rcv ) THEN                                                              ! latent heat flux 38..43)
    do k = 1, size(Ice_boundary%lh_flux,3)
      WHERE (maskt == 0) Ice_boundary%lh_flux(isc:iec,jsc:jec,k) =  frcv(isc:iec,jsc:jec,jn)
      jn = jn + 1
    enddo 
  !ENDIF

  jn = 44
  !IF( nrcvinfo(jn) == OASIS_Rcv ) THEN                                                              ! sensible heat flux 44..49)
    do k = 1, size(Ice_boundary%t_flux,3)
      WHERE (maskt == 0) Ice_boundary%t_flux(isc:iec,jsc:jec,k) =  frcv(isc:iec,jsc:jec,jn)
      jn = jn + 1
    enddo 
  !ENDIF
  
  write(*,*) 'oas_exchange_fields finished.'

END SUBROUTINE oas_exchange_fields
#ENDIF
