#IFDEF COUP_OAS
SUBROUTINE oas_receive_field(Ice,Ice_boundary,Ice_Ocean_Boundary,Time_start,Timet)

USE oas_vardef

USE mpp_mod,         only: mpp_pe, mpp_root_pe
USE time_manager_mod, only : get_time, get_time,operator(-),time_type
USE ice_model_mod,    only: ice_data_type,atmos_ice_boundary_type
USE ocean_model_mod, only: ice_ocean_boundary_type
USE netcdf
USE ice_grid_mod

IMPLICIT NONE

INTEGER         :: &
  maskt(iec-isc+1,jec-jsc+1),  & ! mask array
  maskc(iec-isc+1,jec-jsc+1) ! mask array

type(time_type),       intent(in) :: Time_start,Timet ! for coupling sandra
TYPE (ice_data_type),      INTENT(IN) :: Ice
TYPE (atmos_ice_boundary_type), INTENT(INOUT) :: Ice_boundary
type(ice_ocean_boundary_type), intent(inout) :: Ice_Ocean_Boundary

INTEGER :: &
  isec,sec,day, tcks,       &
  ncfileid, ncvarid,  & ! NetCDF IDs
  sc,dy,jn,k,      &
  nrcvinfo(nfld_rcv_tot)  ! OASIS info argument

REAL(kind=8)      :: &
  ztmp1  (iec,jec)          

  nrcvinfo (:) = OASIS_idle
  ztmp1 (:,:) = 0.0


! call get_time(Ice%Time-Ice%Time_init,sc,dy)
! isec = sc+(dy*86400) 


! call get_time(Ice%Time_step_slow,scs,dys)
!  call get_time(Ice%Time-Time_start, sc, dy)
 call get_time(Timet-Time_start, sc, dy)
  isec = (864e2*dy+sc)!-(864e2*dys+scs)


!  isec     =  (nsteps) * dt_clpd
  
  !IF (ltime) CALL get_timings (i_cpl_add_comp, ntstep, dt, izerror)

!------------------------------------------------------------------------------
! Receive all coupling fields  
!-------------------------------------------------------------------------------


#ifdef OASIS_IOW_ESM
! Do not read here, use instead oas_exchange_fields module
#else
  DO jn = 1, nfld_rcv_tot
    !    IF( srcv(jn)%laction ) THEN
          CALL oas_recieve( jn, isec, ztmp1(:,:), nrcvinfo(jn) )
          ! write(*,*) " nrcvinfo(jn) ", nrcvinfo(jn), jn, OASIS_Rcv, isec
          IF( nrcvinfo(jn) == OASIS_Rcv ) frcv(:,:,jn)=ztmp1(:,:)
    !    ENDIF
      ENDDO
!write(*,*)" oasis receive laction: " , srcv(jn)%laction, isec

!IF (ltime) CALL get_timings (i_cpl_get, ntstep, dt, izerror)

  istatus=nf90_open('masks.nc', NF90_NOWRITE, ncfileid)
  istatus=nf90_inq_varid(ncfileid, 'tmom.msk' , ncvarid)
  istatus=nf90_get_var(ncfileid, ncvarid, maskt, &
         (/ isc, jsc /), (/ iec-isc+1,jec-jsc+1 /))
  istatus=nf90_inq_varid(ncfileid, 'cmom.msk' , ncvarid)
  istatus=nf90_get_var(ncfileid, ncvarid, maskc, &
         (/ isc, jsc /), (/ iec-isc+1,jec-jsc+1 /))

  istatus=nf90_close(ncfileid)

  jn =  1
!  IF( nrcvinfo(jn) == OASIS_Rcv  ) THEN                                                              ! rain 1)
     do k = 1, size(Ice_boundary%lprec,3)
         WHERE (maskt == 0) Ice_boundary%lprec(isc:iec,jsc:jec,k) = frcv(isc:iec,jsc:jec,jn)
     enddo
! ENDIF
  jn = jn + 1 
!  IF( nrcvinfo(jn) == OASIS_Rcv  ) THEN                                                              ! evaporation 2)
     do k = 1, size(Ice_boundary%q_flux,3)
        WHERE (maskt == 0) Ice_boundary%q_flux(isc:iec,jsc:jec,k) =  -frcv(isc:iec,jsc:jec,jn)
     enddo
!  ENDIF 
  jn = jn + 1
!  IF( nrcvinfo(jn) == OASIS_Rcv  ) THEN                                                              ! snow 3)
     do k = 1, size(Ice_boundary%fprec,3)
        WHERE (maskt == 0) Ice_boundary%fprec(isc:iec,jsc:jec,k) = frcv(isc:iec,jsc:jec,jn)
     enddo
!  ENDIF
  jn = jn + 1
!  IF( nrcvinfo(jn) == OASIS_Rcv  ) THEN                                                              ! sea level pressure 4)
     do k = 1, size(Ice_boundary%p,3)
        WHERE (maskt == 0) Ice_boundary%p(isc:iec,jsc:jec,k) =  frcv(isc:iec,jsc:jec,jn)
     enddo
!  ENDIF
  jn = jn + 1
!  IF( nrcvinfo(jn) == OASIS_Rcv  ) THEN                                                              ! u velocity 5)
        WHERE (maskc == 0) Ice_Ocean_Boundary%u_wind(isc:iec,jsc:jec) =  frcv(isc:iec,jsc:jec,jn)
!  ENDIF
  jn = jn + 1
!  IF( nrcvinfo(jn) == OASIS_Rcv  ) THEN                                                              ! v velocity 6)
        WHERE (maskc == 0) Ice_Ocean_Boundary%v_wind(isc:iec,jsc:jec) = frcv(isc:iec,jsc:jec,jn)
!  ENDIF
  jn = jn + 1
!  IF( nrcvinfo(jn) == OASIS_Rcv  ) THEN                                                              ! u wind stress 7)
     do k = 1, size(Ice_boundary%u_flux,3)
        WHERE (maskc == 0) Ice_boundary%u_flux(isc:iec,jsc:jec,k) =  -frcv(isc:iec,jsc:jec,jn)
     enddo
!  ENDIF
  jn = jn + 1
!  IF( nrcvinfo(jn) == OASIS_Rcv  ) THEN                                                              ! v wind stress 8)
     do k = 1, size(Ice_boundary%v_flux,3)
        WHERE (maskc == 0) Ice_boundary%v_flux(isc:iec,jsc:jec,k) =  -frcv(isc:iec,jsc:jec,jn)
     enddo
!  ENDIF
  jn = jn + 1
!  IF( nrcvinfo(jn) == OASIS_Rcv  ) THEN                                                              ! longwave radiation upward 9)
     do k = 1, size(Ice_boundary%lw_flux,3)
          WHERE (maskt == 0) Ice_boundary%lw_flux(isc:iec,jsc:jec,k) =  frcv(isc:iec,jsc:jec,jn)
      enddo
!  ENDIF
  jn = jn + 1
!  IF( nrcvinfo(jn) == OASIS_Rcv  ) THEN                                                              ! longwave radiation downward 10)
     do k = 1, size(Ice_boundary%lw_flux,3)
         WHERE (maskt == 0) Ice_boundary%lw_flux(isc:iec,jsc:jec,k) = -Ice_boundary%lw_flux(isc:iec,jsc:jec,k) + frcv(isc:iec,jsc:jec,jn)
     enddo
! ENDIF
  jn = jn + 1
!  IF( nrcvinfo(jn) == OASIS_Rcv ) THEN                                                              ! shortwave radiation direct 11)
     do k = 1, size(Ice_boundary%sw_flux_vis_dir,3)
         WHERE (maskt == 0) Ice_boundary%sw_flux_vis_dir(isc:iec,jsc:jec,k) =  frcv(isc:iec,jsc:jec,jn)
     enddo
!  ENDIF
    jn = jn + 1
!  IF( nrcvinfo(jn) == OASIS_Rcv ) THEN                                                              ! shortware radiation diffusive 12)
     do k = 1, size(Ice_boundary%sw_flux_vis_dif,3)
         WHERE (maskt == 0) Ice_boundary%sw_flux_vis_dif(isc:iec,jsc:jec,k) =  frcv(isc:iec,jsc:jec,jn)
     enddo
!  ENDIF
  jn = jn + 1
!  IF( nrcvinfo(jn) == OASIS_Rcv ) THEN                                                              ! latent heat flux 13)
     do k = 1, size(Ice_boundary%lh_flux,3)
        WHERE (maskt == 0) Ice_boundary%lh_flux(isc:iec,jsc:jec,k) =  -frcv(isc:iec,jsc:jec,jn)
     enddo 
!  ENDIF
  jn = jn + 1
!  IF( nrcvinfo(jn) == OASIS_Rcv ) THEN                                                              ! sensible heat flux 14)
     do k = 1, size(Ice_boundary%t_flux,3)
         WHERE (maskt == 0) Ice_boundary%t_flux(isc:iec,jsc:jec,k) =  -frcv(isc:iec,jsc:jec,jn)
     enddo 
!  ENDIF
#endif

END SUBROUTINE oas_receive_field
#ENDIF
