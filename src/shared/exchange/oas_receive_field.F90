#ifdef COUP_OAS
SUBROUTINE receive_fld

!!---------------------------------------------------------------------
!!              ***  ROUTINE receive_fld  ***
!!
!! ** Purpose : - Prepare and receive coupling fields to OASIS
!!    
!!----------------------------------------------------------------------

USE oas_cos_vardef

USE data_modelconfig, ONLY : &
  ke, ke1,         &
  dt,              &
  ie, je,          & 
  istartpar,       &
  iendpar,         &
  jstartpar,       &
  jendpar,         &
  ie_tot,          &
  je_tot,          &
  pollat,          &
  pollon,          &
  raddeg,          &
  idt_qv

USE data_parameters,  ONLY : &
  ireals,          &
  iintegers

USE data_parallel,    ONLY : &
  isubpos,         & ! positions of the subdomains in the total domain
  nboundlines,     &     
  my_cart_id,      &
  num_compute,     &
  sendbuf,         &
  isendbuflen,     &
  imp_reals,       &
  icomm_cart,      &
  my_cart_neigh,   &
  ncomm_type,      &
  ldatatypes

USE environment,      ONLY : &
  exchg_boundaries, &
  comm_barrier,     &
  get_free_unit,    &
  release_unit,	    &
  model_abort

USE data_runcontrol , ONLY : &
  lphys,          & ! forecast with physical parametrizations
  llake,          & ! forecst with lake model FLake
  lgsp,           & ! forecast with grid-scale precipitation
  lprog_qi,       & ! if .TRUE., running with cloud ice
  ltur,           & ! forecast with turbulent diffusion
  lcpfluc,        & ! consideration of fluctuations of the heat capacity of air
  itype_turb,     & ! type of turbulent diffusion parametrization
  imode_turb,     & ! mode of turbulent diffusion parametrization
  l2tls,          & ! forecast with 2-TL integration scheme 
  nold  ,         & ! corresponds to ntstep - 1
  nnow  ,         & ! corresponds to ntstep
  nnew  ,         & ! corresponds to ntstep + 1
  ntstep,         & !
  nstart,         & ! first time step of the forecast
  ltime,          &
  nstop,          &
  nbl_exchg,      &
  nbd1, nbd2

USE data_constants,   ONLY : &
  r_d,            & ! gas constant for dry air
  rdv,            & ! r_d / r_v
  o_m_rdv,        & ! 1 - r_d/r_v
  rvd_m_o,        & ! r_v/r_d - 1
  cp_d,           & ! specific heat of dry air at constant pressure
  rdocp,          & ! r_d / cp_d
  cpdr,           & ! 1 / cp_d
  gamma,          & ! 1 / (1 - rdocp)   ( = cp_d/cv_d)
  lh_v,           & ! latent heat of vapourization
  g,              & ! acceleration due to gravity
  gq,             & ! g*g
  gr,             & ! 1 / g
  r_earth,        & ! radius of the earth/acrlat
  sigma,          & ! Boltzmann-constant
  b1,             &
  b2w,            &
  b3,             &
  b4w 

USE data_soil,        ONLY : &
  ctalb       ! thermal albedo (1-emissivity) ?!

USE data_io,          ONLY : &
  ydir_restart      ! directory for restart files

USE data_fields,      ONLY : &
  llandmask  ,    & ! landpoint mask (fr_land>=0.5)
  fr_land,        & ! land fraction
  alb_rad     ,   & ! albedo of the ground                            --
  tch         ,   & ! turbulent transfer coefficient for heat       ( -- )
  tcm         ,   & ! turbulent transfer coefficient for momentum   ( -- )
  tfh         ,   & ! factor of laminar transfer of scalars         ( -- )
  a1t, a2t  ,     & ! implicit weight of vertical diffusion
  hhl       ,     & ! geometrical height of half levels             ( m   )
  hsurf     ,     & ! height of surface topography                  ( m   )
  rlat,           & ! latitudes of the true geographical system     ( rad )
  rlon,           & ! longitudes of the true geographical system    ( rad )
    
  ! 1. constant fields for the reference atmosphere                     (unit)
  ! -----------------------------------------------
  p0         ,    & ! base state pressure                           (Pa)
    
  ! 3. prognostic variables                                             (unit)
  ! -----------------------
  u          ,    & ! zonal wind speed                              ( m/s )
  v          ,    & ! meridional wind speed                         ( m/s )
  t          ,    & ! temperature                                   (  k  )
!  qv         ,    & ! specific water vapor content                  (kg/kg)
  qc         ,    & ! cloud water content                           (kg/kg)
  qi         ,    & ! cloud ice content                             (kg/kg)
  pp         ,    & ! deviation from the reference pressure         ( pa  )
  ps        ,     & ! surface pressure                              ( pa  )

  !    fields for surface values and soil model variables               (unit )
  ! -----------------------------------------------------
  t_g       ,     & ! weighted surface temperature                  (  K  )
  t_s       ,     & ! temperature of the ground surface             (  k  )
  qv_s      ,     & ! specific water vapor content on the surface   (kg/kg)
  fr_ice    ,     & ! ice fraction 
  qvsflx    ,     & ! surface flux of water vapour                  (1/m2s)
  umfl_s    ,     & ! u-momentum flux (surface)                     ( N/m2)
  vmfl_s    ,     & ! v-momentum flux (surface)                     ( N/m2)
  shfl_s    ,     & ! sensible heat flux (surface)                  ( W/m2)
  lhfl_s    ,     & ! latent heat flux (surface)                    ( W/m2)
  t_snow    ,     & ! temperature of the snow-surface
  w_snow    ,     & ! water content of snow
  h_snow    ,     & ! snow height
  w_i       ,     & ! water content of interception water
  runoff_g  ,     & ! surface water runoff
  runoff_s  ,     & ! soil water runoff
  rho_snow  ,     & ! prognostic snow density

  ! 4. boundary fields
  ! -----------------------
  u_bd      ,     &
  v_bd      ,     & 
  t_bd      ,     &
  qv_bd     ,     &
  qc_bd     ,     &
  qi_bd     ,     &
  pp_bd     ,     &
  t_s_bd    ,     &
  t_snow_bd ,     &
  w_snow_bd ,     &
  qv_s_bd

USE src_2wn,          ONLY : &
  ps_cm , ps_gl   , & ! CM surface pressure
  fis_cm, fis_gl  , & ! CM surface geopotential
  fis_lm,           &
  qv_cm , qv_lm   , &
  qc_cm , qc_lm   , &
  qi_cm , qi_lm   , &
  t_cm  , t_lm    , &
  u_cm  , u_lm    , &
  u_prof,           &
  v_cm  , v_lm    , &
  v_prof,           &
          pp_lm   , &
          ps_lm   , &
  w_snow_cm       , &
  w_snow_lm       , &
  t_snow_lm       , &
  t_s_lm          , &
  t_s_cm,  t_s_gl , &
  qv_s_lm         , &
  ke_in,            &
  beta_cm         , &   
  org_vert_interp_2wn, &
  var_names,        &
  nc_diagnose,      &
  nrecord,          &
  rh_s_gl,          &
  read_fis_gl,      &
  write_restart_2wn,&
  read_restart_2wn, &
  convert_uv2uvrot, &
  convert_uvrot2uv, &
  write_2d_ps

USE time_utilities,   ONLY : &
  get_timings,      & 
  i_oasis_cpl,      &
  i_cpl_add_comp,   &
  i_cpl_get,        &
  i_cpl_vert_interp

USE meteo_utilities,  ONLY : &
  qsat,             &
  psat_w

USE src_tracer,      ONLY : trcr_get, trcr_errorstr

USE parallel_utilities,  ONLY:  &
  gather_field,     &
  distribute_field

USE utilities,  ONLY:  &
  uv2uvrot_vec

USE netcdf

USE src_integrals, ONLY : &
  integral_3d_total

USE data_turbulence, ONLY:  &
  vel_min

IMPLICIT NONE

!
! local parameters, variables and arrays
!

INTEGER, PARAMETER ::   jps_taux   =  1    ! zonal wind stress
INTEGER, PARAMETER ::   jps_tauy   =  2    ! meridional wind stress
INTEGER, PARAMETER ::   jps_lat    =  3    ! total latent heat flux (W/m**2)
INTEGER, PARAMETER ::   jps_sens   =  4    ! total sensible heat flux (W/m**2)
INTEGER, PARAMETER ::   jps_ir     =  5    ! emitted infrared (longwave) radiation (W/m**2)                        
INTEGER, PARAMETER ::   jps_alb    =  6    ! albedo

INTEGER, PARAMETER ::   jpss_tf    =  1    ! weighted surface/canopy temperature                                
INTEGER, PARAMETER ::   jpss_ts    =  2    ! temperature of the ground surface                                
INTEGER, PARAMETER ::   jpss_tsnow =  3    ! temperature of the snow surface                          
INTEGER, PARAMETER ::   jpss_alb   =  4    ! albedo
INTEGER, PARAMETER ::   jpss_evtp  =  5    ! surface flux of water vapour
INTEGER, PARAMETER ::   jpss_hstp  =  6    ! sensible heat flux (surface)

INTEGER(KIND=iintegers) :: &
  isec,               &
  isec_2wn,           &
  info,               &
  istatus,            &
  jn, jm, jo,         &
  i, j, k, n,         &
  ncfileid, ncvarid,  & ! NetCDF IDs
  ier,                & ! return error code
  nx,                 & ! time index
  im1, jm1,           & ! i-1, j-1
  ip1, jp1,           & ! i+1, j+1
  izerror,            &
  il_file_id,         &
  mype,               &
  ib,                 &
  npes,               &
  ierror,             &
  status,             &
  isim_time             ! total time (s) of currently simulated time slice

INTEGER(KIND=iintegers) :: &
  kzdims (24),        & 
  il_var_id (6),      & 
  dimids (2),         &
  nrcvinfo(nfld_rcv_tot)  ! OASIS info argument

INTEGER(KIND=4)         :: &
  mask(nlei-nldi+1,nlej-nldj+1) ! mask array

REAL(KIND=ireals)       :: & 
  ztvb, zvbke,        & !
  ztcm,ztch,          &
  za1t_surf,          &
  za2t_surf,          & !
  zfpi, zppa,zpp,     &
  znew, ztmcmq,       &
  dzke, tcm_epsi,     &
  coef_lim, tch_epsi, & ! by Oliver Fuhrer
  integ3d

REAL (KIND=ireals)      :: &
  za1t    (      ke1),     & !
  za2t    (      ke1),     & !
  zpia    (ie,je    ),     & !
  zpianf  (ie,je    ),     & !
  ztmp1   (ie,je    ),     &
  zlon    (ie,je    ),     &
  zlat    (ie,je    ),     &
  ps_tot  (ie_tot,je_tot)

CHARACTER(LEN=25)       :: &
  fmt='(F20.2,1X,F20.2,1X,F20.2)'  ! format descriptor

CHARACTER(LEN=200)      :: &
  yzerrmsg               ! error message for error handling

REAL(KIND=ireals) :: &
  t_s_TRIM0  (nldi:nlei,nldj:nlej) , & ! t_s from TRIM before gathering
  t_s_TRIM1  (ie_tot,je_tot)       , & ! t_s from TRIM after gathering
  t_s_CCLM   (ie_tot,je_tot)       , & ! t_s on CCLM grid
  seamask    (nldi:nlei,nldj:nlej) , & ! land/sea mask of CCLM before gathering
  seamask_tot(1:ie_tot,1:je_tot),    & ! land/sea mask of CCLM after gathering
  var_tot    (ie_tot,je_tot)           !SW: just a temporary field (will be deleted later on)

INTEGER(KIND=4)   ::  &
  trimmask (1:ie_tot,1:je_tot)         ! TRIM domain mask

!SW used for writing profiles to file
REAL (KIND=ireals)        :: &
  zxexp_tot (ie_tot,je_tot), &
  zpexp_tot (ie_tot,je_tot), &
  zp (ie,je,ke)
INTEGER (KIND=iintegers)  :: &
  igp, jgp,            &
  izerrstat,           &
  nuprofile              ! unit number for ascii files
CHARACTER (LEN=6)       :: &
  fmt_ij = '(I3.3)',   & ! format descriptor: an integer of width 3 with preceding zero(s)
  fmt_step = '(I5.5)'    ! format descriptor: an integer of width 5 with preceding zero(s)
CHARACTER (LEN=3)       :: &
  ig, jg,              &
  yflist(3)
CHARACTER (LEN=5)       :: &
  idtstep
CHARACTER(LEN=80) ::       &
  yfname                 ! filename
CHARACTER (LEN=25)       :: yzroutine

REAL (KIND=ireals), POINTER :: &
    qv_now     (:,:,:)   => NULL() !, &     ! QV at nx

!------------------------------------------------------------------------------
!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! 1. Attribute fields to received fields (independent of the specific coupling) 
!-------------------------------------------------------------------------------

! Statement function zfpi for calculation of the exner function
! where the dummy argument zppa is pressure
zfpi(zppa) = (1.E-5_ireals*zppa)**rdocp

   izerror  = 0_iintegers
   yzerrmsg = ''
   yzroutine= 'receive_fld'

   CALL trcr_get(izerror, idt_qv, ptr_tlev = nnow, ptr = qv_now)
   IF (izerror /= 0) THEN
      yzerrmsg = trcr_errorstr(izerror)
      CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
   ENDIF

! Coupling only with PEs that contain at least one land point
IF ( lpe_cpl ) THEN

  nrcvinfo (:) = OASIS_idle
  ztmp1 (:,:) = 0.0_ireals

  isec = nint ( ( ntstep - nstart ) * dt )
  isec_2wn = isec 

  IF (  debug_oasis > 15 .AND. my_cart_id == 0 ) THEN
    WRITE(nulout,*) 'in receive_fld: ntstep=',ntstep,' isec=',isec,' isec_2wn=',isec_2wn !SW, 26.08.13
  ENDIF

  IF (ltime) CALL get_timings (i_cpl_add_comp, ntstep, dt, izerror)

!------------------------------------------------------------------------------
! 2. Receive all coupling fields (independent of the specific coupling) 
!-------------------------------------------------------------------------------

  DO jn = 1, nfld_rcv_tot
    IF( srcv(jn)%laction ) THEN
      ! 1) Receive fields of a land surface model
      IF ( jn >= 1 .AND. jn <= nfld_rcv_lsm ) THEN
        CALL oas_cos_rcv( jn, isec, ztmp1(:,:), nrcvinfo(jn) )
        IF( nrcvinfo(jn) == OASIS_Rcv ) frcv(:,:,jn)=ztmp1(:,:)
      ENDIF
  
      ! 2) Receive fields of ECHAM (only at a coupling time step)
      IF ( (jn >= nfld_rcv_lsm+1) .AND. (jn <= nfld_rcv_lsm+SUM(nlev_rcv_2wn)) ) THEN
        ! only receive at coupling time step, dt_cp, and prevent last receive
        IF ( MOD(isec_2wn,dt_cp) == 0 .AND. ntstep < nstop ) THEN
          CALL oas_cos_rcv( jn, isec_2wn, ztmp1(:,:), nrcvinfo(jn) )
          IF( nrcvinfo(jn) == OASIS_Rcv ) frcv(:,:,jn) = ztmp1(:,:)
        ENDIF
      ENDIF
  
      ! 3) Receive fields of an ocean model
      IF ( jn >= (nfld_rcv_lsm+SUM(nlev_rcv_2wn)+1) ) THEN
        CALL oas_cos_rcv( jn, isec, ztmp1(:,:), nrcvinfo(jn) )
        IF( nrcvinfo(jn) == OASIS_Rcv ) frcv(:,:,jn)=ztmp1(:,:)
      ENDIF
    ENDIF
  ENDDO

ENDIF ! lpe_cpl

IF (ltime) CALL get_timings (i_cpl_get, ntstep, dt, izerror)

!----------------------------------------------------------------------------
! 2.1 clmXX-specific handling of received fields
!----------------------------------------------------------------------------

IF( ytype_lsm == 'clmXX' ) THEN 

  IF ( lpe_cpl ) THEN

    ! Land points updated at each time step
    ! frcv array updated only at coupling time step

    ! COSMO fluxes are positive downwards
    ! in CLM convention differs between 3.5 and 4.0/4.5 !!

    IF( srcv(jps_taux)%laction ) THEN
      WHERE ( llandmask )  umfl_s(:,:) = abs(frcv(:,:,jps_taux))
    ENDIF

    IF( srcv(jps_tauy)%laction ) THEN
      WHERE ( llandmask )  vmfl_s(:,:) = abs(frcv(:,:,jps_tauy))
    ENDIF

    IF( srcv(jps_lat)%laction ) THEN
      WHERE ( llandmask )  lhfl_s(:,:) = frcv(:,:,jps_lat)
    ENDIF

    IF( srcv(jps_sens)%laction ) THEN
      WHERE ( llandmask )  shfl_s(:,:) = frcv(:,:,jps_sens)
    ENDIF

    IF( srcv(jps_alb)%laction ) THEN
      WHERE ( llandmask )  alb_rad(:,:) = frcv(:,:,jps_alb)
    ENDIF

      IF( srcv(jps_ir)%laction ) THEN
         WHERE ( llandmask )
            t_g(:,:,nnew) = (frcv(:,:,jps_ir) / sigma / (1._ireals - ctalb))**0.25
            t_g(:,:,nnow) = t_g(:,:,nnew)
            t_s(:,:,nnow) = t_g(:,:,nnow)
         ENDWHERE
         IF (.NOT. l2tls) THEN
            WHERE ( llandmask )  t_g(:,:,nold) = t_g(:,:,nnew)
         ENDIF
      ENDIF

    ! derive COSMO transfer coefficients (this is a bit of reverse engineering
    ! since COSMO does not exchange fluxes but calculates them based on updated coefficients
    ! and states)

      ! select timelevel and timestep for calculations
      nx  = nnow

      do j=jstartpar,jendpar

       jm1  = MAX( 1, j-1 )

       do i=istartpar,iendpar

          if (llandmask(i,j)) then     ! for land-points only
             im1        = MAX( 1, i-1)

             ! Calculation for potential temperature gradient at surface
             zpp         = p0(i,j,ke) + pp(i,j,ke,nnow)
             dzke  = hhl(i,j,ke) - hhl(i,j,ke1)
             zpia(i,j) = zfpi( zpp )
             zpianf(i,j) = zfpi( ps(i,j,nx) )

             ! Calculation for modified transfer coefficients ztcm/ztch
             zvbke      = MAX( 0.5*SQRT ( (u(i,j,ke,nx) + u(im1,j,ke,nx))**2    &
                 + (v(i,j,ke,nx) + v(i,jm1,ke,nx))**2 ), vel_min)
             ztvb       = t_s (i,j,nx)*(1.0 + rvd_m_o*qv_s(i,j,nx))
             ! formulation below is the correct one, 
             ! but should be made consistent in slow_tendencies as well
!             ztvb       = t_g (i,j,nx)*(1.0 + rvd_m_o*qv_s(i,j,nx))
!             ztcm       = tcm(i,j)*zvbke*g*ps(i,j,nx)/(r_d*ztvb)
!             ztch       = tch(i,j)*zvbke*g*ps(i,j,nx)/(r_d*ztvb)

             !-----------------------------------------------------------------
             ! sensible heat flux inversion: derive new tch from CLM flux
             !-----------------------------------------------------------------

             tch(i,j) = -shfl_s (i,j) / ( &
                  (zvbke*ps(i,j,nx)/(r_d*ztvb))*cp_d* &
                  ( t_g(i,j,nx) - zpianf(i,j)*t(i,j,ke,nx)/zpia(i,j) ) )

             !-------------------------------------------------------------------
             ! calculate effective surface humidity based on CLM latent heat flux
             !-------------------------------------------------------------------

             ztch = tch(i,j)*zvbke*ps(i,j,nx)/(r_d*ztvb)
             IF (ABS(ztch) >= 1e-20_ireals ) THEN  ! ELD
                qv_s(i,j,nx) = qv_now(i,j,ke) - lhfl_s(i,j) / (lh_v * ztch)
             ELSE
                qv_s(i,j,nx) = qv_now(i,j,ke)
             END IF

             ! limit qv_s to positive values (no negative surface humidity ;-)
             qv_s(i,j,nx) = min(max(1.E-5_ireals,qv_s(i,j,nx)),1.E-1_ireals)
             ! set qv_s(nnew)=qv_s(nnow) similarly to what is done in TERRA
             ! this will ensure that the flux rediagnosed in COSMO is the same as here
             qv_s(i,j,nnew) = qv_s(i,j,nx)
             IF (.NOT. l2tls) qv_s(i,j,nold) = qv_s(i,j,nx)  ! ELD

             !----------------------------------------------------------------
             ! momentum flux inversion: derive new tcm from CLM momentum flux
             !----------------------------------------------------------------

             ! tcm needs to be always positive
             if (u(i,j,ke,nx) == 0) then
                tcm(i,j)=0.
             else
                tcm(i,j) = umfl_s(i,j) / ( zvbke*g*ps(i,j,nx)/(r_d*ztvb) * gr * u(i,j,ke,nx) )
             !tcm(i,j) = vmfl_s(i,j) / ( zvbke*g*ps(i,j,nx)/(r_d*ztvb) * gr * v(i,j,ke,nx) )
             endif

             !write(6,*), umfl_s(i,j)/u(i,j,ke,nx),  vmfl_s(i,j)/v(i,j,ke,nx) 

             ! implement limiter so exchange coeff. is not too large (Oliver Fuhrer, 2007)
             ! NOTE: might add a check if tcm<0 here!
             tcm_epsi=0.5   
             coef_lim=tcm_epsi*dzke/(zvbke*dt)
             if (abs(tcm(i,j))>coef_lim) then
!!                write(6,*) 'WARNING: TCM limiter!!!',tcm(i,j),coef_lim
                tcm(i,j)=sign(coef_lim,tcm(i,j))
             end if

         end if !llandmask
       end do   !i
      end do    !j
  ENDIF ! lpe_cpl

  ! Update SSTs here because tgcom is not called anymore
  WHERE ( .NOT. llandmask )   t_g(:,:,nnew) = t_s(:,:,nnew)

ENDIF ! ytype_lsm == 'clmXX'

!----------------------------------------------------------------------------
! 2.2 veg3d-specific handling of received fields
!----------------------------------------------------------------------------

IF( ytype_lsm == 'veg3d' ) THEN 

  IF ( lpe_cpl ) THEN

      ! veg3d fluxes are postive upwards
      ! COSMO fluxes are positive downwards
      IF( srcv(jpss_hstp)%laction .and. nrcvinfo(jpss_hstp) == OASIS_Rcv ) frcv(:,:,jpss_hstp) = -frcv(:,:,jpss_hstp)
      IF( srcv(jpss_evtp)%laction .and. nrcvinfo(jpss_evtp) == OASIS_Rcv ) frcv(:,:,jpss_evtp) = -frcv(:,:,jpss_evtp)


      IF( srcv(jpss_tf)%laction ) THEN
         WHERE ( fr_land >= 0.5 )  t_g(:,:,nnew) = frcv(:,:,jpss_tf)
      ENDIF

      IF( srcv(jpss_ts)%laction ) THEN
         WHERE ( fr_land >= 0.5 )  t_s(:,:,nnew) = frcv(:,:,jpss_ts)
      ENDIF

      IF( srcv(jpss_tsnow)%laction ) THEN
         WHERE ( fr_land >= 0.5 )  t_snow(:,:,nnew) = frcv(:,:,jpss_tsnow)
      ENDIF

      IF( srcv(jpss_alb)%laction ) THEN
         WHERE ( fr_land >= 0.5 )  alb_rad(:,:) = frcv(:,:,jpss_alb)
      ENDIF

      IF( srcv(jpss_evtp)%laction ) THEN
         WHERE ( fr_land >= 0.5 )  qvsflx(:,:) = frcv(:,:,jpss_evtp)
      ENDIF

      IF( srcv(jpss_hstp)%laction ) THEN
         WHERE ( fr_land >= 0.5 )  shfl_s(:,:) = frcv(:,:,jpss_hstp)
      ENDIF

    ! derive COSMO transfer coefficients (this is a bit of reverse engineering
    ! since COSMO does not exchange fluxes but calculates them based on updated coefficients
    ! and states)

    ! select timelevel and timestep for calculations
    nx  = nnow

    ! most code below comes from slow_tendencies.f90
    IF (ltur.EQV. .TRUE. .AND. (itype_turb /= 3 .OR. imode_turb < 2)) THEN
      DO k = 1, ke1
        za1t(k) = a1t(k)
        za2t(k) = a2t(k)
      ENDDO
    ELSE
      DO k = 1, ke1
        za1t(k) = 0.0_ireals
        za2t(k) = 0.0_ireals
      ENDDO
    END IF

    ! Selection of lower boundary conditions
    IF (imode_turb == 0) THEN
      ! condition of the lower boundary is the surface concentration
      za1t_surf = za1t(ke1)  ! implicit weight for ke-value
      za2t_surf = za2t(ke1)  ! explicit weight for ke-value
    ELSE
      ! condition of the lower boundary is the explicit surface mass flux
      za1t_surf = 0.0_ireals ! no implicit weight
      za2t_surf = 1.0_ireals ! full explicit weight
    ENDIF

    DO j = jstartpar, jendpar

      jm1  = MAX( 1, j-1 )
!      jp1  = MIN( je, j+1 )

      DO i = istartpar, iendpar

        IF ( llandmask(i,j) ) THEN     ! for land-points only
          im1 = MAX(  1, i-1)
!          ip1 = MIN( ie, i+1)

          t_g(i,j,nx) = t_g(i,j,nnew)
          t_s(i,j,nx) = t_s(i,j,nnew)

          ! Calculation for potential temperature gradient at surface
          zpp         = p0(i,j,ke) + pp(i,j,ke,nnow)
          dzke        = hhl(i,j,ke) - hhl(i,j,ke1)
          zpia(i,j)   = zfpi( zpp )
          zpianf(i,j) = zfpi( ps(i,j,nx) )

          ! Calculation for modified transfer coefficients ztcm/ztch
          zvbke      = MAX(0.5*SQRT ( (u(i,j,ke,nx) + u(im1,j,ke,nx))**2     &
                           +(v(i,j,ke,nx) + v(i,jm1,ke,nx))**2 ), vel_min)
          ztvb       = t_s (i,j,nx)*(1.0 + rvd_m_o*qv_s(i,j,nx))

          tch(i,j) = -shfl_s (i,j) / ( &
                     zvbke*g*ps(i,j,nx)/(r_d*ztvb)*gr*cp_d* &
                     ( t_g(i,j,nx) - zpianf(i,j)*t(i,j,ke,nx)/zpia(i,j) ) )

          ! latent heat flux inversion: derive new qv_s from CLM flux
          ztch  = tch(i,j)*zvbke*g*ps(i,j,nx)/(r_d*ztvb)

          !lhfl_s derive from qvsflx
          lhfl_s(i,j) = qvsflx(i,j) * lh_v
          qv_s(i,j,nx) = -qvsflx(i,j) / ( ztch*gr)  + qv_now(i,j,ke)

          ! limit qv_s to positive values (no negative surface humidity ;-)
          qv_s(i,j,nx) = MIN(MAX(1.E-5_ireals,qv_s(i,j,nx)),1.E-1_ireals)
          qv_s(i,j,nnew) = qv_s(i,j,nx)
          IF (.NOT. l2tls) qv_s(i,j,nold) = qv_s(i,j,nx)  ! ED

        ENDIF !llandmask

      ENDDO

   ENDDO

  END IF   ! coupling on valid PE only
 
   ! Update SSTs here because tgcom is not called anymore
  WHERE ( fr_land < 0.5 )   t_g(:,:,nnew) = t_s(:,:,nnew)
  !! ELD: replace with llandmask!

ENDIF ! ytype_lsm = veg3d

!----------------------------------------------------------------------------
! 2.3 echam-specific handling of received fields
!----------------------------------------------------------------------------

IF( ytype_2wn == 'echam' .AND. lpe_cpl ) THEN

  !----------------------------------------------------------------------------
  ! 2.3.1 Update the _cm fields only case of a coupling time step 
  !----------------------------------------------------------------------------

  ! Coupled fields in this order
  ! '_PS', 'BET', 'SGP', 'WSN', '_TS', '__T', '__U', '__V', '_QV', '_QC', '_QI'
  jn = nfld_rcv_lsm 
  DO jo = 1, nfld_rcv_2wn
    DO jm = 1, nlev_rcv_2wn(jo)
      jn = jn + 1
      IF( srcv(jn)%laction .AND. nrcvinfo(jn) == OASIS_Rcv ) THEN
        IF ( jo == 1 ) ps_cm (nldi:nlei,nldj:nlej) = frcv(nldi:nlei,nldj:nlej,jn) ! BET NUR ZURUECK
        IF ( jo == 2 .AND. isec_2wn == 0 ) fis_cm (nldi:nlei,nldj:nlej) = frcv(nldi:nlei,nldj:nlej,jn) ! BET NUR ZURUECK
        IF ( jo == 3 ) w_snow_cm (nldi:nlei,nldj:nlej) = frcv(nldi:nlei,nldj:nlej,jn) ! BET NUR ZURUECK
        IF ( jo == 4 ) t_s_cm (nldi:nlei,nldj:nlej) = frcv(nldi:nlei,nldj:nlej,jn) ! BET NUR ZURUECK
        IF ( jo == 5 ) t_cm (nldi:nlei,nldj:nlej,jm) = frcv(nldi:nlei,nldj:nlej,jn) ! BET NUR ZURUECK
        IF ( jo == 6 ) u_cm (nldi:nlei,nldj:nlej,jm) = frcv(nldi:nlei,nldj:nlej,jn) ! BET NUR ZURUECK
        IF ( jo == 7 ) v_cm (nldi:nlei,nldj:nlej,jm) = frcv(nldi:nlei,nldj:nlej,jn) ! BET NUR ZURUECK
        IF ( jo == 8 ) qv_cm (nldi:nlei,nldj:nlej,jm) = frcv(nldi:nlei,nldj:nlej,jn) ! BET NUR ZURUECK
        IF ( jo == 9 ) qc_cm (nldi:nlei,nldj:nlej,jm) = frcv(nldi:nlei,nldj:nlej,jn) ! BET NUR ZURUECK
        IF ( jo == 10 ) qi_cm (nldi:nlei,nldj:nlej,jm) = frcv(nldi:nlei,nldj:nlej,jn) ! BET NUR ZURUECK
      ENDIF
    ENDDO
  ENDDO

  !----------------------------------------------------------------------------
  ! 2.3.2 Calculate variables that are used for computing the time interpolation
  !       weights in initialize_loop.
  !----------------------------------------------------------------------------

  IF ( (ntstep+1 > nnext_cp) .AND. (ntstep < nstop) ) THEN
    nbd1_cp = 3 - nbd1_cp
    nbd2_cp = 3 - nbd2_cp

    hlast_cp = hnext_cp
    hnext_cp = hlast_cp + hinc_cp
    nlast_cp = nnext_cp
    nnext_cp = NINT ( 3600.0_ireals * hnext_cp / dt )
    ninc_cp  = nnext_cp - nlast_cp

    IF ( debug_oasis > 15 .AND. my_cart_id == 0 ) THEN
      WRITE(nulout,*) '------------------------------------'
      WRITE(nulout,*) 'update time variables in receive_fld:'
      WRITE(nulout,*) 'ntstep=',ntstep
      WRITE(nulout,*) 'nbd1_cp=',nbd1_cp,' nbd2_cp=',nbd2_cp
      WRITE(nulout,*) 'hlast_cp, hnext_cp, hinc_cp: ',hlast_cp, hnext_cp, hinc_cp
      WRITE(nulout,*) 'nlast_cp, nnext_cp, ninc_cp: ',nlast_cp, nnext_cp, ninc_cp
      WRITE(nulout,*) '------------------------------------'
      CALL flush(nulout)
    ENDIF
  ENDIF

  !----------------------------------------------------------------------------
  ! 2.3.3 Perform the vertical interpolation and update of boundary data only 
  !       in case of a coupling time step
  !----------------------------------------------------------------------------

  IF ( MOD(isec_2wn,dt_cp) == 0 .AND. ntstep < nstop ) THEN
     
    IF (  debug_oasis > 15 .AND. my_cart_id == 0) THEN
      WRITE(nulout,*) 'received coupling fields from ECHAM at ntstep=',ntstep
    ENDIF

    ! count the coupling time steps
    dt_cp_num = dt_cp_num + 1
    nrecord = nrecord + 1

    ! Check humidity quantities for unreasonable values due to horizontal
    ! interpolation
    WHERE ( qv_cm(:,:,:) < 0.0_ireals ) qv_cm(:,:,:) = 0.0_ireals
    WHERE ( qv_cm(:,:,:) > 1.0_ireals ) qv_cm(:,:,:) = 1.0_ireals
    WHERE ( qc_cm(:,:,:) < 0.0_ireals ) qc_cm(:,:,:) = 0.0_ireals
    WHERE ( qc_cm(:,:,:) > 1.0_ireals ) qc_cm(:,:,:) = 1.0_ireals
    IF ( lprog_qi ) THEN
      WHERE ( qi_cm(:,:,:) < 0.0_ireals ) qi_cm(:,:,:) = 0.0_ireals
      WHERE ( qi_cm(:,:,:) > 1.0_ireals ) qi_cm(:,:,:) = 1.0_ireals
    ENDIF

    ! copy received fields to working arrays
    t_lm  (:,:,:) = t_cm  (:,:,:)
    u_lm  (:,:,:) = u_cm  (:,:,:)
    v_lm  (:,:,:) = v_cm  (:,:,:)
    qv_lm (:,:,:) = qv_cm (:,:,:)
    qc_lm (:,:,:) = qc_cm (:,:,:)
    qi_lm (:,:,:) = qi_cm (:,:,:)
    ps_gl (:,:)   = ps_cm (:,:)
    fis_gl(:,:)   = fis_cm(:,:)
    t_s_gl(:,:)   = t_s_cm(:,:)
    w_snow_lm(:,:)= w_snow_cm(:,:)

    ! exchange boundaries of all received fields
    IF ( num_compute > 1 ) THEN
      kzdims(1:24) = (/ ke_in,ke_in,ke_in,ke_in,ke_in,ke_in,           &
        ke_in,ke_in,ke_in,ke_in,ke_in,ke_in,1,1,1,1,0,0,0,0,0,0,0,0 /)
 !! naveed check argument  
  !  CALL exchg_boundaries                                            &
  !      (0, sendbuf, isendbuflen, imp_reals, icomm_cart, ie, je,       &
  !      kzdims, jstartpar, jendpar, 3, nboundlines, my_cart_neigh,     &
  !      20000+ntstep, .FALSE., ncomm_type, izerror, yzerrmsg,          &
  !      t_cm (:,:,:),  u_cm(:,:,:),  v_cm(:,:,:),                      & 
  !      qv_cm(:,:,:), qc_cm(:,:,:), qi_cm(:,:,:),                      &
  !      t_lm (:,:,:),  u_lm(:,:,:),  v_lm(:,:,:),                      &
  !      qv_lm(:,:,:), qc_lm(:,:,:), qi_lm(:,:,:),                      &
  !      ps_gl(:,:), fis_gl(:,:), t_s_gl(:,:), w_snow_lm(:,:))
    ENDIF

CALL gather_field(fis_gl(:,:),ie,je,zxexp_tot,ie_tot,je_tot,0,nerror)
IF (my_cart_id==0) PRINT *, 'mean(fis_gl)=',SUM(zxexp_tot)/(MAX(1,SIZE(zxexp_tot)))
CALL gather_field(fis_lm(:,:),ie,je,zxexp_tot,ie_tot,je_tot,0,nerror)
IF (my_cart_id==0) PRINT *, 'mean(fis_lm)=',SUM(zxexp_tot)/(MAX(1,SIZE(zxexp_tot)))

    ! compute the relative humidity at the surface; this will be needed for
    ! computing qv_s in SUBROUTINE ground_fields
    rh_s_gl(:,:) = qv_lm(:,:,ke_in)
    DO j = 1, je
      DO i = 1, ie
        rh_s_gl(i,j) = rh_s_gl(i,j) /  &
          qsat(psat_w(t_s_gl(i,j), b1, b2w, b3, b4w), ps_gl(i,j), rdv, o_m_rdv)
      ENDDO
    ENDDO

    ! Write some "debug" output
    CALL gather_field (ps_gl(:,:),ie,je,ps_tot,ie_tot,je_tot,0,nerror)
    CALL write_2d_ps (ps_tot,'PS','ps_rcv_1.nc','write',1,nrecord,1)
    IF ( nrecord < 800 ) CALL nc_diagnose ('add_var', 'vor_vi.nc', ke_in, 1, nrecord)
    IF ( nrecord < 800 ) CALL nc_diagnose ('add_var_2d', 'vor_vi.nc', ke_in, 1, nrecord)

    IF (ltime) CALL get_timings (i_cpl_add_comp, ntstep, dt, izerror)

    !----------------------------------------------------------------------------
    ! 2.3.4 Convert the U and V velocities from the real geographical to the
    !       rotated coordinate system (see int2lm: SUBROUTINE interpol_coarse_uv)
    !----------------------------------------------------------------------------

    ! uncomment these lines to check the errors introduced by the back-and-forth
    ! conversion
    !DO k = 1, ke_in
    !  CALL gather_field (u_lm(:,:,k),ie,je,ps_tot,ie_tot,je_tot,0,nerror)
    !  CALL write_2d_ps (ps_tot,'T','ps_rcv_1.nc','write',1,nrecord,k)
    !  CALL gather_field (v_lm(:,:,k),ie,je,ps_tot,ie_tot,je_tot,0,nerror)
    !  CALL write_2d_ps (ps_tot,'T','ps_rcv_2.nc','write',1,nrecord,k)
    !ENDDO
    !CALL convert_uv2uvrot
    !CALL convert_uvrot2uv
    !DO k = 1, ke_in
    !  CALL gather_field (u_lm(:,:,k),ie,je,ps_tot,ie_tot,je_tot,0,nerror)
    !  CALL write_2d_ps (ps_tot,'T','ps_snd_1.nc','write',1,nrecord,k)
    !  CALL gather_field (v_lm(:,:,k),ie,je,ps_tot,ie_tot,je_tot,0,nerror)
    !  CALL write_2d_ps (ps_tot,'T','ps_snd_2.nc','write',1,nrecord,k)
    !ENDDO

    IF ( ltransform_uv ) THEN
      CALL convert_uv2uvrot
   
      ! Save the converted u and v velocities for later use in src_2wn for plugging
      ! back in parts of the original profiles (see SUBROUTINE vert_interpol_r2g) 
      u_prof(:,:,:) = u_lm(:,:,:)
      v_prof(:,:,:) = v_lm(:,:,:)
 
      IF (ltime) CALL get_timings (i_cpl_add_comp, ntstep, dt, izerror)
    ENDIF

    !----------------------------------------------------------------------------
    ! 2.3.5 Do the vertical interpolation from ECHAM to CCLM levels
    !----------------------------------------------------------------------------

    CALL org_vert_interp_2wn ('g2r')

    IF (ltime) CALL get_timings (i_cpl_vert_interp, ntstep, dt, izerror)

    IF ( num_compute > 1 ) THEN
      kzdims(1:24) = (/ ke_in,ke_in,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /)
  !! naveed --    
  !CALL exchg_boundaries                                            &
  !      (0, sendbuf, isendbuflen, imp_reals, icomm_cart, ie, je,       &
  !      kzdims, jstartpar, jendpar, 3, nboundlines, my_cart_neigh,     &
  !      20000+ntstep, .FALSE., ncomm_type, izerror, yzerrmsg,          &
  !      u_lm(:,:,:),  v_lm(:,:,:))
    ENDIF

    !----------------------------------------------------------------------------
    ! 2.3.6 Copy interpolated fields to the future boundary time level
    !----------------------------------------------------------------------------

    t_bd (:,:,:,nbd2_cp) = t_lm (:,:,1:ke)
    u_bd (:,:,:,nbd2_cp) = u_lm (:,:,1:ke)
    v_bd (:,:,:,nbd2_cp) = v_lm (:,:,1:ke)
    qv_bd(:,:,:,nbd2_cp) = qv_lm(:,:,1:ke)
    qc_bd(:,:,:,nbd2_cp) = qc_lm(:,:,1:ke)
    qi_bd(:,:,:,nbd2_cp) = qi_lm(:,:,1:ke)
    pp_bd(:,:,:,nbd2_cp) = pp_lm(:,:,1:ke)
    t_s_bd (:,:,nbd2_cp) = t_s_lm (:,:)
    qv_s_bd (:,:,nbd2_cp) = qv_s_lm (:,:)
    t_snow_bd (:,:,nbd2_cp) = t_snow_lm (:,:)
    w_snow_bd (:,:,nbd2_cp) = w_snow_lm (:,:)

    ! write interpolated fields to some diagnoses output for evaluation against true lbfd files
    IF ( nrecord < 100 ) CALL nc_diagnose ('add_nbd2', 'nbd2.nc', ke, 2, nrecord)

    !----------------------------------------------------------------------------
    ! 2.3.7 Calculate volume integrals
    !----------------------------------------------------------------------------

    IF ( lmeasure_mass ) THEN
      !integ3d = integral_3d_total (qv_bd(:,:,:,nbd2_cp),3)
      integ3d = integral_3d_total (qv_cm(:,:,:),3)
      IF (my_cart_id==0) THEN
        OPEN (numassqv,FILE='total_mass_qv.txt',FORM='FORMATTED',ACTION='WRITE',POSITION='APPEND',STATUS='OLD',IOSTAT=izerrstat)
        WRITE (numassqv,'(A,1X,I10,1X,I10,1X,F30.5)') 'rcv',ntstep,isec_2wn,integ3d
        CLOSE (numassqv,STATUS='KEEP',IOSTAT=izerrstat)
      ENDIF
      !integ3d = integral_3d_total (qc_bd(:,:,:,nbd2_cp),3)
      integ3d = integral_3d_total (qc_cm(:,:,:),3)
      IF (my_cart_id==0) THEN
        OPEN (numassqc,FILE='total_mass_qc.txt',FORM='FORMATTED',ACTION='WRITE',POSITION='APPEND',STATUS='OLD',IOSTAT=izerrstat)
        WRITE (numassqc,'(A,1X,I10,1X,I10,1X,F30.5)') 'rcv',ntstep,isec_2wn,integ3d
        CLOSE (numassqc,STATUS='KEEP',IOSTAT=izerrstat)
      ENDIF
      !integ3d = integral_3d_total (qi_bd(:,:,:,nbd2_cp),3)
      integ3d = integral_3d_total (qi_cm(:,:,:),3)
      IF (my_cart_id==0) THEN
        OPEN (numassqi,FILE='total_mass_qi.txt',FORM='FORMATTED',ACTION='WRITE',POSITION='APPEND',STATUS='OLD',IOSTAT=izerrstat)
        WRITE (numassqi,'(A,1X,I10,1X,I10,1X,F30.5)') 'rcv',ntstep,isec_2wn,integ3d
        CLOSE (numassqi,STATUS='KEEP',IOSTAT=izerrstat)
      ENDIF
    ENDIF

    !----------------------------------------------------------------------------
    ! 2.3.8 Special handling of the boundary time level nbd1_cp at the first 
    !       time step of a time segment, e.g. one month, in case of initialisation
    !       or restart
    !----------------------------------------------------------------------------

    IF ( ntstep == nstart .AND. ntstep /= 0 ) THEN
      ! in case of a restart read the nbd1 fields from a file
      CALL read_restart_2wn (ke, ydir_restart)
    ELSEIF ( ntstep == 0 ) THEN
      ! use the fields from the initilisation for nbd1_cp in case of a 
      ! simulation start
      t_bd (:,:,:,nbd1_cp) = t (:,:,:,nnew)
      u_bd (:,:,:,nbd1_cp) = u (:,:,:,nnew)
      v_bd (:,:,:,nbd1_cp) = v (:,:,:,nnew)
      qv_bd(:,:,:,nbd1_cp) = qv_now(:,:,:)
!jb check      qc_bd(:,:,:,nbd1_cp) = qc(:,:,:,nnew)
!jb check      qi_bd(:,:,:,nbd1_cp) = qi(:,:,:,nnew)
      pp_bd(:,:,:,nbd1_cp) = pp(:,:,:,nnew)

      t_s_bd    (:,:,nbd1_cp) = t_s    (:,:,nnew)
      qv_s_bd   (:,:,nbd1_cp) = qv_s   (:,:,nnew)
      t_snow_bd (:,:,nbd1_cp) = t_snow (:,:,nnew)
      w_snow_bd (:,:,nbd1_cp) = w_snow (:,:,nnew)

      IF ( lvert_profiles ) THEN
        ! calculate pressure on CCLM full levels
        zp(:,:,:) = p0(:,:,:) + pp (:,:,:,nnew)

        ! hard-coded list of variables for which vertical profiles will be written
        yflist = (/'t  ','u  ','v  '/)
        WRITE (ig, fmt_ij) gppi
        WRITE (jg, fmt_ij) gppj

        var_loop : DO j = 1, 3

          ! create file name
          yfname = TRIM(yflist(j))//'_'//ig//'_'//jg//'_00000_init'

          CALL get_free_unit (nuprofile)
          IF (my_cart_id == 0) THEN
            OPEN (nuprofile,FILE=yfname,FORM='FORMATTED',ACTION='WRITE',STATUS='NEW',IOSTAT=izerrstat)
          ENDIF

          DO k = 1, ke
            IF ( TRIM(yflist(j))=='t' ) THEN
              CALL gather_field(t(:,:,k,nnew),ie,je,zxexp_tot,ie_tot,je_tot,0,nerror)
            ELSEIF ( TRIM(yflist(j))=='u' ) THEN
              CALL gather_field(u(:,:,k,nnew),ie,je,zxexp_tot,ie_tot,je_tot,0,nerror)
            ELSEIF ( TRIM(yflist(j))=='v' ) THEN
              CALL gather_field(v(:,:,k,nnew),ie,je,zxexp_tot,ie_tot,je_tot,0,nerror)
            ENDIF
            CALL gather_field(zp  (:,:,k),ie,je,zpexp_tot,ie_tot,je_tot,0,nerror)
            IF ( my_cart_id == 0 ) THEN
              WRITE (nuprofile,'(I3,1X,F15.5,1X,F15.5)') k,zpexp_tot(gppi,gppj),zxexp_tot(gppi,gppj)
            ENDIF
          ENDDO

          IF (my_cart_id == 0) THEN
            CLOSE (nuprofile, STATUS='KEEP', IOSTAT=izerrstat)
            IF (izerrstat /= 0) THEN
              PRINT *, ' ERROR *** while opening or closing file ',yfname
              RETURN
            ENDIF
          ENDIF
          CALL release_unit (nuprofile)

        ENDDO var_loop
      ENDIF ! lvert_profiles
    ENDIF ! ntstep == 0

    !----------------------------------------------------------------------------
    ! 2.3.9 Write the received fields at the last coupling time step to a file;
    !       this file will be read at the restart of the next time segment
    !----------------------------------------------------------------------------

    isim_time = nint( ( nstop - nstart ) * dt )
    IF ( isec == isim_time - dt_cp ) THEN
      CALL write_restart_2wn (ke, ydir_restart)
    ENDIF 

! This is called at the end of this routine
!    IF (ltime) CALL get_timings (i_cpl_add_comp, ntstep, dt, izerror)

  ENDIF ! IF (OASIS_Rcv)
ENDIF ! ytype_2wn

!----------------------------------------------------------------------------
! 2.4 nemoO- and romsO- and trimI-specific handling of received fields
!----------------------------------------------------------------------------

IF( (ytype_oce == 'romsO' .OR. ytype_oce == 'nemoO' .OR. ytype_oce == 'trimI'  .OR. ytype_oce == 'nemoD' .OR. ytype_oce == 'momBS' ) & ! sandra 02022018
    .AND. lpe_cpl ) THEN
! UOI, HaHo 2013-09-04 }

  istatus=nf90_open('masks.nc', NF90_NOWRITE, ncfileid)
  IF ( ytype_oce == 'nemoO' .OR. ytype_oce == 'nemoD' ) istatus=nf90_inq_varid(ncfileid, 'coms.msk' , ncvarid)
  IF ( ytype_oce == 'trimI' ) istatus=nf90_inq_varid(ncfileid, 'cobs.msk' , ncvarid)
  IF ( ytype_oce == 'romsO' ) istatus=nf90_inq_varid(ncfileid, 'coas.msk' , ncvarid)
  istatus=nf90_get_var(ncfileid, ncvarid, mask, &
         (/ isubpos(my_cart_id,1)-iboundstart,         &
            isubpos(my_cart_id,2)-jboundstart /),      &
          (/ nlei-nldi+1,nlej-nldj+1 /))
  istatus=nf90_close(ncfileid)

  ! Update non-masked t_s with frcv array
  ! SST index
  jn = nfld_rcv_lsm + SUM ( nlev_rcv_2wn )
  jn = jn + 1
! IF( srcv(jn)%laction .AND. nrcvinfo(jn) == OASIS_Rcv ) THEN
IF( srcv(jn)%laction ) THEN 
   IF( ytype_oce == 'trimI') THEN
      IF(CPL_FLG .gt. 0) THEN     ! since CURRENT_DATE == CPL_START, otherwise SST is from ERA-int/40 ...
        IF ( my_cart_id == 0 ) &
          write(0,*) "+++cclm_4.8_clm19_uoi_SW_20131007_14Jan14_varsCCLM: CCLM uses TRIM's SST"
        WHERE (mask == 0) t_s(nldi:nlei,nldj:nlej,nnow) = frcv(nldi:nlei,nldj:nlej,jn)
        WHERE (mask == 0) t_s(nldi:nlei,nldj:nlej,nnew) = frcv(nldi:nlei,nldj:nlej,jn)
      ELSE
        IF ( my_cart_id == 0 ) &
          write(0,*) "+++cclm_4.8_clm19_uoi_SW_20131007_14Jan14_varsCCLM: CCLM uses ERAint's SST"
      ENDIF
    ELSE
      WHERE (mask == 0) t_s(nldi:nlei,nldj:nlej,nnow) = frcv(nldi:nlei,nldj:nlej,jn)
      WHERE (mask == 0) t_s(nldi:nlei,nldj:nlej,nnew) = frcv(nldi:nlei,nldj:nlej,jn)
    ENDIF ! ytype_oce == 'trimI'
  ENDIF  ! srcv(jn)%laction
ENDIF ! ytype_oce == 'nemoO' .OR. 'romsO' .OR. 'trimI'

IF( ( ytype_oce == 'nemoI' .OR. ytype_oce == 'nemoD' .OR. ytype_oce == 'momBS') &  ! sandra 02022018
    .AND. lpe_cpl ) THEN
! UOI, HaHo 2013-09-04 }

  istatus=nf90_open('masks.nc', NF90_NOWRITE, ncfileid)
  istatus=nf90_inq_varid(ncfileid, 'cobs.msk' , ncvarid)
  istatus=nf90_get_var(ncfileid, ncvarid, mask, &
         (/ isubpos(my_cart_id,1)-iboundstart,         &
            isubpos(my_cart_id,2)-jboundstart /),      &
          (/ nlei-nldi+1,nlej-nldj+1 /))
  istatus=nf90_close(ncfileid)

  ! Update non-masked t_s with frcv array
  ! SST index
  jn = nfld_rcv_lsm + SUM ( nlev_rcv_2wn )
  jn = jn + 1
!IF( srcv(jn)%laction .AND. nrcvinfo(jn) == OASIS_Rcv ) THEN
IF( srcv(jn)%laction ) THEN
    WHERE (mask == 0) t_s(nldi:nlei,nldj:nlej,nnow) = frcv(nldi:nlei,nldj:nlej,jn)
    WHERE (mask == 0) t_s(nldi:nlei,nldj:nlej,nnew) = frcv(nldi:nlei,nldj:nlej,jn)
  ENDIF
  jn = jn + 1 
 ! IF( srcv(jn)%laction .AND. nrcvinfo(jn) == OASIS_Rcv ) &
  IF( srcv(jn)%laction ) &
    WHERE (mask == 0) fr_ice(nldi:nlei,nldj:nlej) =  frcv(nldi:nlei,nldj:nlej,jn)
ENDIF ! ytype_oce == 'nemoI' ! and 'momBS' sandra 02022018

IF (ltime) CALL get_timings (i_cpl_add_comp, ntstep, dt, izerror)

END SUBROUTINE receive_fld
#endif
