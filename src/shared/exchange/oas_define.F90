#ifdef COUP_OAS
SUBROUTINE oas_cos_define
!
! grid definitions for OASIS
!

USE oas_vardef
USE mod_oasis_namcouple    ! OASIS3-MCT namcouple variables: e.g. coupling time step
USE ice_grid_mod 
USE mpp_domains_mod, only: mpp_get_compute_domain, mpp_get_compute_domains, &
                           Domain2d, mpp_get_global_domain 
USE mpp_mod,         only: mpp_pe, mpp_root_pe

IMPLICIT NONE

INTEGER               :: write_aux_files=1

CHARACTER(len=4)       :: &
  cl_grd_t = 'momt',        & ! name of t grid of mom
  cl_grd_c = 'momc'           ! name of c grid of mom

!INTEGER                           :: isc, iec, jsc, jec ! compute domain
!INTEGER                           :: isg, ieg, jsg, jeg ! global domain

!CHARACTER(len=30)       :: &
!  data_gridname='grids.nc', & ! file with the grids
!  data_maskname='masks.nc'    ! file with the masks


  call mpp_get_compute_domain( Domain, isc, iec, jsc, jec )
  call mpp_get_global_domain( Domain, isg, ieg, jsg, jeg )

! -----------------------------------------------------------------
! ... Write auxiliary files by the master processor 
! -----------------------------------------------------------------

IF (mpp_pe().EQ.mpp_root_pe()) THEN

  CALL oasis_start_grids_writing(write_aux_files)
      
  CALL oasis_write_grid (momt, iec, jec, geo_lon, geo_lat, part_id)
  CALL oasis_write_grid (momc, iec, jec, geo_lonv_ib, geo_latv_ib, part_id)

  CALL oasis_write_mask(momt, iec, jec, wett, part_id)
  CALL oasis_write_mask(momc, iec, jec, wetv, part_id)

  CALL oasis_write_area(momt, iec, jec ,cell_area, part_id)
  CALL oasis_write_area(momc, iec, jec ,cell_area, part_id)

  CALL oasis_terminate_grids_writing()


ENDIF


 ! oasis pertitioning

  il_paral(1) = 2  ! box partitioning
  ! Global extent in x
  il_paral(5) = ieg 
  ! Upper left corner global offset
  il_paral(2) = (iec-1) + ((jec-1) * il_paral(5))
  ! Local extent in x
  il_paral(3) = iec-isc+1 
  ! Local extent in y
  il_paral(4) = jec-jsc+1


CALL oasis_def_partition(part_id, il_paral, ierror)

  IF (ierror /= OASIS_Success) THEN
    CALL oasis_abort(comp_id, 'oas_define', 'Failure in oasis_def_partition' )
  ENDIF

! -----------------------------------------------------------------
! ... Initialize some variables
! ----------------------------------------------------------------

  nfld_snd_oce = 0
  nfld_rcv_oce = 0
  nfld_snd_tot = 0
  nfld_rcv_tot = 0

! -----------------------------------------------------------------
! Define list of RECIEVE variables per coupling
! -----------------------------------------------------------------
    nfld=0

    nfld_rcv_oce = 12
    ALLOCATE ( nlev_rcv_oce (nfld_rcv_oce), STAT=ierror )
    ALLOCATE ( nam_rcv_oce  (nfld_rcv_oce), STAT=ierror )
    nlev_rcv_oce = (/1,1,1, &
                     1,1,1, &
                     1,1, &
                     1,1,1,1/)
    nam_rcv_oce(1) =  'OCEREPRE'      ! precipitation
    nam_rcv_oce(2) =  'OCEREEVA'      ! evaporation
    nam_rcv_oce(3) =  'OCERESNO'      ! snow
    nam_rcv_oce(4) =  'OCERESLP'      ! sea level pressure
    nam_rcv_oce(5) =  'OCERETBO'      ! 2m temperature (temperature at bottom)
    nam_rcv_oce(6) =  'OCERESHU'      ! specific humidity
    nam_rcv_oce(7) =  'OCEREUVL'      ! u velocity 10m
    nam_rcv_oce(8) =  'OCEREVVL'      ! v velocity 10m
    nam_rcv_oce(9) =  'OCERELWR'      ! net longwave radiation
    nam_rcv_oce(10) =  'OCERESWR'     ! net shortwave radiation
    nam_rcv_oce(11) =  'OCERELHF'     ! latent heaf flux
    nam_rcv_oce(12) =  'OCERESHF'     ! sensible heat flux

    DO ji = 1, nfld_rcv_oce
      nfld = nfld + nlev_rcv_oce(ji)
    ENDDO

    ! total number of fields to be received:
    nfld_rcv_tot = nfld


    nfld = 0

    nfld_snd_oce = 2
    ALLOCATE ( nlev_snd_oce (nfld_snd_oce), STAT=ierrstat )
    ALLOCATE ( nam_snd_oce  (nfld_snd_oce), STAT=ierrstat )
    nlev_snd_oce = (/1,1/)
    nam_snd_oce(1)  = 'OCESESST'
    nam_snd_oce(2) =  'OCESEICE'
    DO ji = 1, nfld_snd_oce
      nfld = nfld + nlev_snd_oce(ji)
    ENDDO

  ! total number of fields to be sent:
  nfld_snd_tot = nfld

! -----------------------------------------------------------------
! ... write variables names for all active couplings to one structure
! ----------------------------------------------------------------

  ! Allocate memory for data exchange:
  ALLOCATE( ssnd(nfld_snd_tot), stat = ierror )
  ALLOCATE( srcv(nfld_rcv_tot), stat = ierror )
  ALLOCATE( exfld(isc:iec, jsc:jec), stat = ierror )

  IF ( ierror > 0 ) THEN
    CALL oasis_abort( comp_id, 'oas_define', 'Failure in allocating exfld, ssnd or srcv' )
    RETURN
  ENDIF

  ! fill ssnd with names of fields to be sent
  nfld = 0

    DO ji = 1, nfld_snd_oce
      nfld = nfld + 1
      ssnd(nfld)%vari_name = TRIM(nam_snd_oce(ji))
    ENDDO

  ! fill srcv with names of fields to be received:
  nfld = 0

    DO ji = 1, nfld_rcv_oce
      nfld = nfld + 1
      srcv(nfld)%vari_name = TRIM(nam_rcv_oce(ji))
    ENDDO

  var_nodims(1) = 2   ! Dimension number of exchanged arrays
  var_nodims(2) = 1   ! number of bundles (always 1 for OASIS3)

  ishape(1,1) = 1
  ishape(2,1) = iec-isc+1
  ishape(1,2) = 1
  ishape(2,2) = jec-jsc+1

    ! Announce variables to be sent:
  DO ji = 1, nfld_snd_tot
    CALL oasis_def_var( ssnd(ji)%vari_id, ssnd(ji)%vari_name, part_id, &
      var_nodims, OASIS_Out, ishape, OASIS_REAL, ierror )
    IF ( ierror /= OASIS_Success ) CALL oasis_abort( ssnd(ji)%vari_id, &
      'oas_define', 'Failure in oasis_def_var for '//TRIM(ssnd(ji)%vari_name) )
  ENDDO

 ji = 1, nfld_rcv_tot
    CALL oasis_def_var( srcv(ji)%vari_id, srcv(ji)%vari_name, part_id, &
      var_nodims, OASIS_In, ishape, OASIS_REAL, ierror )
    IF ( ierror /= OASIS_Success ) CALL oasis_abort( srcv(ji)%vari_id, &
      'oas_define', 'Failure in oasis_def_var for '//TRIM(srcv(ji)%vari_name) )
  ENDDO

  ! Allocate array to store received fields between two coupling steps
  ALLOCATE( frcv(iec, jec, nfld_rcv_tot), stat = ierror )
  IF ( ierror > 0 ) THEN
    CALL oasis_abort( comp_id, 'oas_define', 'Failure in allocating frcv' )
    RETURN
  ENDIF

  CALL oasis_enddef(ierror)

  IF ( ierror /= OASIS_Success ) CALL oasis_abort ( comp_id, &
    'oas_define', 'Failure in oasis_enddef')


! -----------------------------------------------------------------
! ... Deallocate arrays that are not needed anymore
! -----------------------------------------------------------------

  IF (ALLOCATED (nlev_snd_lsm)) THEN
    DEALLOCATE ( nlev_snd_lsm, STAT=ierrstat )
    DEALLOCATE ( nam_snd_lsm , STAT=ierrstat )
  ENDIF
  IF (ALLOCATED (nlev_rcv_lsm)) THEN
    DEALLOCATE ( nlev_rcv_lsm, STAT=ierrstat )
    DEALLOCATE ( nam_rcv_lsm , STAT=ierrstat )
  ENDIF

  IF (ALLOCATED (nam_snd_2wn)) DEALLOCATE ( nam_snd_2wn, STAT=ierrstat )
  IF (ALLOCATED (nam_rcv_2wn)) DEALLOCATE ( nam_rcv_2wn, STAT=ierrstat )

  IF (ALLOCATED (nlev_snd_oce)) THEN
    DEALLOCATE ( nlev_snd_oce, STAT=ierrstat )
    DEALLOCATE ( nam_snd_oce , STAT=ierrstat )
  ENDIF
  IF (ALLOCATED (nlev_rcv_oce)) THEN
    DEALLOCATE ( nlev_rcv_oce, STAT=ierrstat )
    DEALLOCATE ( nam_rcv_oce , STAT=ierrstat )
  ENDIF

ENDIF ! lpe_cpl 

CALL MPI_Barrier(local_comm, ierror)

END SUBROUTINE oas_define
#endif
