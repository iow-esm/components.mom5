#IFDEF COUP_OAS
SUBROUTINE oas_define(local_comm)
  !
  ! grid definitions for OASIS
  !

  USE oas_vardef
  USE mod_oasis_namcouple    ! OASIS3-MCT namcouple variables: e.g. coupling time step
  USE ice_grid_mod 
  USE mpp_domains_mod, only: mpp_get_compute_domain, domain2D, &
                            mpp_get_global_domain, mpp_global_field !, &
  !                           XUPDATE,YUPDATE
  USE mpp_parameter_mod,  only : GLOBAL_ROOT_ONLY, XUPDATE, YUPDATE 
  USE mpp_mod,         only: mpp_pe, mpp_root_pe
  USE ice_type_mod,    only: ice_data_type
  USE netcdf
  USE constants_mod,    only: radius, pi
  USE fms_mod,          only: read_data

  IMPLICIT NONE

  INTEGER               :: write_aux_files=1, ji,k,m
  INTEGER, INTENT(IN)   :: local_comm

  !TYPE (ice_data_type),      INTENT(IN) :: Ice

  CHARACTER(len=128)                  :: grid_file

  CHARACTER(len=4)       :: &
    cl_grd_t = 'tmom' ,        & ! name of t grid of mom
    cl_grd_c = 'cmom'           ! name of c grid of mom

  INTEGER                           :: isg, ieg, jsg, jeg ! global domain

  INTEGER, ALLOCATABLE ::   tmsk(:,:) , vmsk(:,:)
  REAL(KIND=8), ALLOCATABLE   :: tmp1(:,:),tmp2(:,:), x_T(:,:), &
                                y_T(:,:),x_C(:,:),y_C(:,:),  &
                                area_t(:,:), area_c(:,:), wet(:,:), &
                                x_vert(:,:,:), y_vert(:,:,:), &
                                x_vert_tot(:,:,:), y_vert_tot(:,:,:) 


  call mpp_get_compute_domain(Domain, isc, iec, jsc, jec )
  call mpp_get_global_domain( Domain, isg, ieg, jsg, jeg )

  ALLOCATE(tmsk(isg:ieg,jsg:jeg), stat = ierror )
  ALLOCATE(vmsk(isg:ieg,jsg:jeg), stat = ierror )
  ALLOCATE(wet(isc:iec,jsc:jec), stat = ierror )
  ALLOCATE(area_t(isc:iec,jsc:jec), stat = ierror )
  ALLOCATE(area_c(isc:iec,jsc:jec), stat = ierror )
  ALLOCATE(tmp1(isg:ieg,jsg:jeg), stat = ierror )
  ALLOCATE(tmp2(isg:ieg,jsg:jeg), stat = ierror )
  ALLOCATE(x_T(isc:iec,jsc:jec), stat = ierror )
  ALLOCATE(y_T(isc:iec,jsc:jec), stat = ierror )
  ALLOCATE(x_C(isc:iec,jsc:jec), stat = ierror )
  ALLOCATE(y_C(isc:iec,jsc:jec), stat = ierror )
  ALLOCATE(x_vert_tot(isg:ieg,jsg:jeg,4), stat = ierror )
  ALLOCATE(y_vert_tot(isg:ieg,jsg:jeg,4), stat = ierror )
  ALLOCATE(x_vert(isc:iec,jsc:jec,4), stat = ierror )
  ALLOCATE(y_vert(isc:iec,jsc:jec,4), stat = ierror )

  grid_file = 'INPUT/grid_spec.nc'

  tmsk(:,:)=0.
  vmsk(:,:)=0.
  tmp1(:,:)=0.
  tmp2(:,:)=0.

  CALL read_data(grid_file, 'wet', wet,  Domain)

  CALL mpp_global_field(Domain, wet, tmp1, flags=XUPDATE+YUPDATE+GLOBAL_ROOT_ONLY )  
  CALL mpp_global_field(Domain, wet, tmp2, flags=XUPDATE+YUPDATE+GLOBAL_ROOT_ONLY )  

  DO k = isg,ieg
    DO m = jsg,jeg
      IF (tmp1(k,m) <= 0.5) THEN
        tmsk(k,m)=1
      END IF
  !    IF (tmp2(k,m) <= 0.5) THEN
  !      vmsk(k,m)=1
  !    END IF
    END DO
  END DO


  ! -----------------------------------------------------------------
  ! ... paritioning 
  ! -----------------------------------------------------------------


  il_paral(1) = 2  ! box partitioning
  ! Global extent in x
  il_paral(5) = ieg-isg+1 
  ! Upper left corner global offset
  il_paral(2) = (isc-1) + ((jsc-1) * il_paral(5))
  ! Local extent in x
  il_paral(3) = iec-isc+1 
  ! Local extent in y
  il_paral(4) = jec-jsc+1

!  write(*,*) "il_paral",  il_paral(5),il_paral(2), il_paral(3), il_paral(4), (jeg-isg+1)
  CALL oasis_def_partition(part_id, il_paral, ierror, (ieg-isg+1)*(jeg-jsg+1))

  IF (ierror /= OASIS_Success) THEN
    CALL oasis_abort(comp_id, 'oas_define', 'Failure in oasis_def_partition' )
  ENDIF

 !write(*,*) "finished oasis patitioning"


!-----------------------------------------------------------
! ... Write auxiliary files by the master processor 
! -----------------------------------------------------------------

  CALL oasis_start_grids_writing(write_aux_files)

  tmp1(:,:)=0.
  tmp2(:,:)=0.
  
!  CALL mpp_global_field(Domain, x_T, tmp1, flags=XUPDATE+YUPDATE+GLOBAL_ROOT_ONLY )
!  CALL mpp_global_field(Domain, x_T, tmp2,flags=XUPDATE+YUPDATE+GLOBAL_ROOT_ONLY  )
  IF (mpp_pe().EQ.mpp_root_pe()) THEN
    CALL read_data(grid_file, 'x_T', tmp1)
    CALL read_data(grid_file, 'y_T', tmp2)
    CALL oasis_write_grid(cl_grd_t, ieg-isg+1, jeg-jsg+1, tmp1, tmp2)
  ENDIF

  tmp1(:,:)=0.
  tmp2(:,:)=0.
!  CALL mpp_global_field(Domain, x_C, tmp1, flags=XUPDATE+YUPDATE+GLOBAL_ROOT_ONLY )
!  CALL mpp_global_field(Domain, y_C, tmp2,flags=XUPDATE+YUPDATE+GLOBAL_ROOT_ONLY  )
  IF (mpp_pe().EQ.mpp_root_pe()) THEN
    CALL read_data(grid_file, 'x_C', tmp1)
    CALL read_data(grid_file, 'y_C', tmp2)
    CALL oasis_write_grid(cl_grd_c, ieg-isg+1, jeg-jsg+1, tmp1, tmp2)
  ENDIF

! -----------------------------------------------------------------
!  get the corners grid points
! -----------------------------------------------------------------


  x_vert(:,:,:)=0
  y_vert(:,:,:)=0
  x_vert_tot(:,:,:)=0
  y_vert_tot(:,:,:)=0


!  CALL mpp_global_field(Domain, x_vert, x_vert_tot,flags=XUPDATE+YUPDATE+GLOBAL_ROOT_ONLY )
!  CALL mpp_global_field(Domain, y_vert, y_vert_tot,flags=XUPDATE+YUPDATE+GLOBAL_ROOT_ONLY  )

  IF (mpp_pe().EQ.mpp_root_pe()) THEN
    CALL read_data(grid_file, 'y_vert_T', y_vert_tot)
    CALL read_data(grid_file, 'x_vert_T', x_vert_tot)
    CALL oasis_write_corner(cl_grd_t, ieg-isg+1, jeg-jsg+1, 4,  x_vert_tot, y_vert_tot)
  ENDIF

  x_vert(:,:,:)=0
  y_vert(:,:,:)=0
  x_vert_tot(:,:,:)=0
  y_vert_tot(:,:,:)=0


!  CALL mpp_global_field(Domain, x_vert, x_vert_tot,flags=XUPDATE+YUPDATE+GLOBAL_ROOT_ONLY )
!  CALL mpp_global_field(Domain, y_vert, y_vert_tot,flags=XUPDATE+YUPDATE+GLOBAL_ROOT_ONLY  )

  IF (mpp_pe().EQ.mpp_root_pe()) THEN
    CALL read_data(grid_file, 'y_vert_C', y_vert_tot)
    CALL read_data(grid_file, 'x_vert_C', x_vert_tot)
    CALL oasis_write_corner(cl_grd_c, ieg-isg+1, jeg-jsg+1, 4,  x_vert_tot, y_vert_tot)
  ENDIF

! write masks
   
  IF (mpp_pe().EQ.mpp_root_pe()) THEN
    CALL oasis_write_mask(cl_grd_c, ieg-isg+1, jeg-jsg+1, tmsk)
    CALL oasis_write_mask(cl_grd_t, ieg-isg+1, jeg-jsg+1, tmsk)
  ENDIF

! write areas


  tmp1(:,:)=0.
  tmp2(:,:)=0.
 !  CALL mpp_global_field(Domain, area_t, tmp1, flags=XUPDATE+YUPDATE+GLOBAL_ROOT_ONLY )
 !  CALL mpp_global_field(Domain, area_c, tmp2, flags=XUPDATE+YUPDATE+GLOBAL_ROOT_ONLY )
  
  IF (mpp_pe().EQ.mpp_root_pe()) THEN
   CALL read_data(grid_file, 'AREA_T', tmp1)
   CALL read_data(grid_file, 'AREA_C', tmp2)
   CALL oasis_write_area(cl_grd_t, ieg-isg+1, jeg-jsg+1 ,tmp1*10000.) ! *1000 to get to cm^2
   CALL oasis_write_area(cl_grd_c, ieg-isg+1, jeg-jsg+1 ,tmp2*10000.)
  ENDIF

    CALL oasis_terminate_grids_writing()

!write(*,*) "finished oasis file writings"

! -----------------------------------------------------------------
! ... write variables names for all active couplings to one structure
! ----------------------------------------------------------------

  ! Store info which fields to send
  ALLOCATE( ssnd(nfld_snd_tot), stat = ierror )
#ifdef OASIS_IOW_ESM
  ssnd( 1)%var_name = MOM5_instance_letter//'SFARE01' ! cell area fraction of this ice class (1=open sea, 2..6=ice)
  ssnd( 2)%var_name = MOM5_instance_letter//'SFARE02'
  ssnd( 3)%var_name = MOM5_instance_letter//'SFARE03'
  ssnd( 4)%var_name = MOM5_instance_letter//'SFARE04'
  ssnd( 5)%var_name = MOM5_instance_letter//'SFARE05'
  ssnd( 6)%var_name = MOM5_instance_letter//'SFARE06'
  ssnd( 7)%var_name = MOM5_instance_letter//'STSUR01' ! surface temperature of this ice class (1=open sea, 2..6=ice)
  ssnd( 8)%var_name = MOM5_instance_letter//'STSUR02'
  ssnd( 9)%var_name = MOM5_instance_letter//'STSUR03'
  ssnd(10)%var_name = MOM5_instance_letter//'STSUR04'
  ssnd(11)%var_name = MOM5_instance_letter//'STSUR05'
  ssnd(12)%var_name = MOM5_instance_letter//'STSUR06'
  ssnd(13)%var_name = MOM5_instance_letter//'SALBE01' ! shortwave albedo of this ice class (1=open sea, 2..6=ice)
  ssnd(14)%var_name = MOM5_instance_letter//'SALBE02'
  ssnd(15)%var_name = MOM5_instance_letter//'SALBE03'
  ssnd(16)%var_name = MOM5_instance_letter//'SALBE04'
  ssnd(17)%var_name = MOM5_instance_letter//'SALBE05'
  ssnd(18)%var_name = MOM5_instance_letter//'SALBE06'
#endif

  ! Store info which fields to receive
  ALLOCATE( srcv(nfld_rcv_tot), stat = ierror )
#ifdef OASIS_IOW_ESM
  srcv( 1)%var_name = 'MRMRAI01'      ! precipitation
  srcv( 2)%var_name = 'MRMEVA01'      ! evaporation
  srcv( 3)%var_name = 'MRMEVA02'      ! evaporation
  srcv( 4)%var_name = 'MRMEVA03'      ! evaporation
  srcv( 5)%var_name = 'MRMEVA04'      ! evaporation
  srcv( 6)%var_name = 'MRMEVA05'      ! evaporation
  srcv( 7)%var_name = 'MRMEVA06'      ! evaporation
  srcv( 8)%var_name = 'MRMSNO01'      ! snow
  srcv( 9)%var_name = 'MRPSUR01'      ! sea level pressure
  srcv(10)%var_name = 'MRU10M01'      ! u velocity 10m
  srcv(11)%var_name = 'MRV10M01'      ! v velocity 10m
  srcv(12)%var_name = 'MRUMOM01'      ! u wind stress
  srcv(13)%var_name = 'MRUMOM02'      ! u wind stress
  srcv(14)%var_name = 'MRUMOM03'      ! u wind stress
  srcv(15)%var_name = 'MRUMOM04'      ! u wind stress
  srcv(16)%var_name = 'MRUMOM05'      ! u wind stress
  srcv(17)%var_name = 'MRUMOM06'      ! u wind stress
  srcv(18)%var_name = 'MRVMOM01'      ! v wind stress
  srcv(19)%var_name = 'MRVMOM02'      ! v wind stress
  srcv(20)%var_name = 'MRVMOM03'      ! v wind stress
  srcv(21)%var_name = 'MRVMOM04'      ! v wind stress
  srcv(22)%var_name = 'MRVMOM05'      ! v wind stress
  srcv(23)%var_name = 'MRVMOM06'      ! v wind stress
  srcv(24)%var_name = 'MRRBBR01'      ! longwave radiation upward
  srcv(25)%var_name = 'MRRBBR02'      ! longwave radiation upward
  srcv(26)%var_name = 'MRRBBR03'      ! longwave radiation upward
  srcv(27)%var_name = 'MRRBBR04'      ! longwave radiation upward
  srcv(28)%var_name = 'MRRBBR05'      ! longwave radiation upward
  srcv(29)%var_name = 'MRRBBR06'      ! longwave radiation upward
  srcv(30)%var_name = 'MRRLWD01'      ! longwave radiation downward
  srcv(31)%var_name = 'MRRSDD01'      ! shortwave radiation downward direct
  srcv(32)%var_name = 'MRRSIN01'      ! shortwave radiation downward diffusive
  srcv(33)%var_name = 'MRHLAT01'      ! latent heat flux
  srcv(34)%var_name = 'MRHLAT02'      ! latent heat flux
  srcv(35)%var_name = 'MRHLAT03'      ! latent heat flux
  srcv(36)%var_name = 'MRHLAT04'      ! latent heat flux
  srcv(37)%var_name = 'MRHLAT05'      ! latent heat flux
  srcv(38)%var_name = 'MRHLAT06'      ! latent heat flux
  srcv(39)%var_name = 'MRHSEN01'      ! sensible heat flux
  srcv(40)%var_name = 'MRHSEN02'      ! sensible heat flux
  srcv(41)%var_name = 'MRHSEN03'      ! sensible heat flux
  srcv(42)%var_name = 'MRHSEN04'      ! sensible heat flux
  srcv(43)%var_name = 'MRHSEN05'      ! sensible heat flux
  srcv(44)%var_name = 'MRHSEN06'      ! sensible heat flux

#else
  srcv( 1)%var_name = 'OCEREPRE'      ! precipitation
  srcv( 2)%var_name = 'OCEREEVA'      ! evaporation
  srcv( 3)%var_name = 'OCERESNO'      ! snow
  srcv( 4)%var_name = 'OCERESLP'      ! sea level pressure
  srcv( 5)%var_name = 'OCEREUVL'      ! u velocity 10m
  srcv( 6)%var_name = 'OCEREVVL'      ! v velocity 10m
  srcv( 7)%var_name = 'OCEREUWS'      ! u wind stress
  srcv( 8)%var_name = 'OCEREVWS'      ! v wind stress
  srcv( 9)%var_name = 'OCERELWU'      ! longwave radiation upward
  srcv(10)%var_name = 'OCERELWD'      ! longwave radiation downward
  srcv(11)%var_name = 'OCERESWD'      ! shortwave radiation downward direct
  srcv(12)%var_name = 'OCERESWF'      ! shortwave radiation downward diffusive
  srcv(13)%var_name = 'OCERELHF'      ! latent heat flux
  srcv(14)%var_name = 'OCERESHF'      ! sensible heat flux
#endif

  ! Allocate field for temporary storage of coupler data
  ALLOCATE( exfld1(isc:iec, jsc:jec), stat = ierror )

  IF ( ierror > 0 ) THEN
    CALL oasis_abort( comp_id, 'oas_define', 'Failure in allocating exfld, ssnd or srcv' )
    RETURN
  ENDIF

  ! This is still preliminary: set laction for all (this way no subsets of
  ! variables can be selected; this has to be changed yet)
  ssnd(:)%laction = .TRUE.
  srcv(:)%laction = .TRUE.

  var_nodims(1) = 2   ! Dimension number of exchanged arrays
  var_nodims(2) = 1   ! number of bundles (always 1 for OASIS3)

  ishape(1,1) = 1
  ishape(2,1) = iec-isc+1
  ishape(1,2) = 1
  ishape(2,2) = jec-jsc+1

    ! Announce variables to be sent:
  DO ji = 1, nfld_snd_tot
!    write(*,*) "var_name", ssnd(ji)%var_name
    CALL oasis_def_var( ssnd(ji)%vari_id, ssnd(ji)%var_name, part_id, &
      var_nodims, OASIS_Out, ishape, OASIS_REAL, ierror )
    IF ( ierror /= OASIS_Success ) CALL oasis_abort( ssnd(ji)%vari_id, &
      'oas_define', 'Failure in oasis_def_var for '//TRIM(ssnd(ji)%var_name) )
  ENDDO

  DO ji = 1, nfld_rcv_tot
    CALL oasis_def_var( srcv(ji)%vari_id, srcv(ji)%var_name, part_id, &
      var_nodims, OASIS_In, ishape, OASIS_REAL, ierror )
    IF ( ierror /= OASIS_Success ) CALL oasis_abort( srcv(ji)%vari_id, &
      'oas_define', 'Failure in oasis_def_var for '//TRIM(srcv(ji)%var_name) )
  ENDDO

  ! Allocate array to store received fields between two coupling steps
  ALLOCATE( frcv(iec, jec, nfld_rcv_tot), stat = ierror )
  IF ( ierror > 0 ) THEN
    CALL oasis_abort( comp_id, 'oas_define', 'Failure in allocating frcv' )
    RETURN
  ENDIF

!write(*,*), "finished oasis var definition"   

  CALL oasis_enddef(ierror)

  IF ( ierror /= OASIS_Success ) CALL oasis_abort ( comp_id, &
    'oas_define', 'Failure in oasis_enddef')

DEALLOCATE(wet,tmsk,vmsk,tmp1,tmp2,x_vert,y_vert,x_vert_tot,y_vert_tot,  &
           x_T,y_T,x_C,y_C, area_c,area_t,STAT=ierror )

CALL MPI_Barrier(local_comm, ierror)

END SUBROUTINE oas_define
#ENDIF
