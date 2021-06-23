#ifdef COUP_OAS
SUBROUTINE oas_send_field

!!---------------------------------------------------------------------
!!              ***  ROUTINE send_fld  ***
!!
!! ** Purpose : Prepare and send coupling fields to OASIS
!!      
!!----------------------------------------------------------------------

USE oas_vardef

USE mpp_mod,         only: mpp_pe, mpp_root_pe

USE data_fields,      ONLY : &
  hhl,           & ! geometrical height of half levels             ( m   )
  hsurf    ,     & ! height of surface topography                  ( m   )
  ! 1. constant fields for the reference atmosphere                     (unit)
  ! -----------------------------------------------
  p0        ,    & ! base state pressure                           (Pa)
  ! 3. prognostic variables                                             (unit)
  ! -----------------------
  u         ,    & ! zonal wind speed                              ( m/s )
  v         ,    & ! meridional wind speed                         ( m/s )
  t         ,    & ! temperature                                   (  k  )
!  qv        ,    & ! specific water vapor content                  (kg/kg)
  qi        ,    & ! specific cloud water content                  (kg/kg)
  qc        ,    & ! specific cloud ice content                    (kg/kg)
  qr        ,    & ! specific rain content                         (kg/kg)
  qs        ,    & ! specific snow content                         (kg/kg)
  pp        ,    & ! deviation from the reference pressure         ( pa  )
  ! 6. fields that are computed in the parametrization and dynamics     (unit )
  ! ---------------------------------------------------------------
  !   fields of convective and grid-scale precipitation
  prr_con    ,   & ! precipitation rate of rain, convective        (kg/m2*s)
  prs_con    ,   & ! precipitation rate of snow, convective        (kg/m2*s)
  prr_gsp    ,   & ! precipitation rate of rain, grid-scale        (kg/m2*s)
  prs_gsp    ,   & ! precipitation rate of snow, grid-scale        (kg/m2*s)
  prg_gsp    ,   & ! precipitation rate of graupel, grid-scale     (kg/m2*s)
 !   fields of the radiation
 swdir_s     ,   & ! direct shortwave downward radiation           ( W/m2)
 swdifd_s    ,   & ! diffuse shortwave downward radiation          ( W/m2)
 lwd_s       ,   & ! thermal downward radiation at the ground      ( W/m2)
 lhfl_s      ,   & ! latent heat flux (surface)                    ( W/m2)
 shfl_s      ,   & ! sensible heat flux (surface)                  ( W/m2)
 sobs        ,   & ! solar radiation at the ground                 ( W/m2)
 thbs        ,   & ! thermal radiation at the ground               ( W/m2)
 !   momentum
 umfl_s      ,   & ! u-momentum flux (surface)                     ( N/m2)
 vmfl_s      ,   & ! v-momentum flux (surface)                     ( N/m2)
 ! variables needed by the coupling with veg3d
 lai         ,   & ! leaf area index of plants
 plcov       ,   & ! fraction of plant cover
 ps          ,   & ! surface pressure
 rlon        ,   & ! geographical longitude
 rlat        ,   & ! geographical latitude
 soiltyp     ,   &
 ! variables needed by the coupling with trimI
 u_10m       ,   & ! zonal wind at 10m                      ( m/s )
 v_10m       ,   & ! meridional wind at 10m                 ( m/s )
 t_2m        ,   & ! temperature at 2m                      ( K   )
 rh_2m       ,   & ! relative humidity at 2m                ( %   )
 clct        ,   & ! total cloud cover                      ( --  )
 hhl         ,   & ! geometrical hgt of half levs           ( m   )
 rhoA        ,   & ! low level total density of air         (kg/m3)
 Tpot        ,   & ! low level potential air temperature    ( K   )
 ! UOI, HaHo 2013-09-04 {
 pmsl1       ,   & ! mean sea level pressure                ( Pa  )
 ! UOI, HaHo 2013-09-04 }
 ! variables additionally needed by the coupling with momBS  ! sandra 02022018
 qv_2m             ! specific humidity at 2m                (kg/kg)
 
USE netcdf

USE src_2wn,          ONLY : &
  ke_in,               &
  ke_out,              &
  t_lm, t_cm,          &
  u_lm, u_cm,          &  
  v_lm, v_cm,          &
  qv_lm, qv_cm,        &
  qc_lm, qc_cm,        &
  qi_lm, qi_cm,        & 
  pp_lm, ps_lm,        & 
  org_vert_interp_2wn, &
  ps_cm, ps_gl,        & ! CM surface pressure
  beta_cm,             &
  nc_diagnose,         &
  nrecord,             &
  gradients_2wn,       &
  dtdx, dtdy,          &
  dpsdx, dpsdy,        &
  dudx, dudy,          &
  dvdx, dvdy,          &
  dqvdx, dqvdy,        &
  dqcdx, dqcdy,        &
  dqidx, dqidy,        &
!  u_lm_rot,            &
!  v_lm_rot,            &
  convert_uvrot2uv,    &
  convert_uv2uvrot,    &
  print_profile,       &
  akh_in, bkh_in,      &
  write_2d_ps

USE time_utilities,   ONLY : & 
  get_timings,       & 
  i_cpl_put,         &
  i_cpl_vert_interp, &
  i_cpl_add_comp,    &
  i_cpl_grad

USE src_tracer,      ONLY : trcr_get, trcr_errorstr

USE environment,      ONLY : model_abort

USE utilities,  ONLY : & 
  phirot2phi,          &
  uvrot2uv_vec,        &
  uv2uvrot_vec

USE parallel_utilities, ONLY:  &
  gather_field,        &
  distribute_field

USE src_integrals, ONLY: &
  integral_3d_total

USE data_turbulence, ONLY:  &
  vel_min

IMPLICIT NONE

!
! local parameters, variables and arrays
!

INTEGER, PARAMETER ::   jps_t   =  1      ! temperature
INTEGER, PARAMETER ::   jps_u   =  2      ! u wind
INTEGER, PARAMETER ::   jps_v   =  3      ! v wind
INTEGER, PARAMETER ::   jps_q   =  4      ! specific water vapor content
INTEGER, PARAMETER ::   jps_th  =  5      ! thickness of lowest level (m)
INTEGER, PARAMETER ::   jps_pr  =  6      ! surface pressure (Pa)
INTEGER, PARAMETER ::   jps_rs  =  7      ! direct shortwave downward radiation (W/m2)
INTEGER, PARAMETER ::   jps_fs  =  8      ! diffuse shortwave downward radiation (W/m2)
INTEGER, PARAMETER ::   jps_lw  =  9      ! longwave downward radiation (W/m2)
INTEGER, PARAMETER ::   jps_cr  = 10      ! convective rain precipitation      (kg/m2*s)
INTEGER, PARAMETER ::   jps_cs  = 11      ! convective snow precipitation      (kg/m2*s)
INTEGER, PARAMETER ::   jps_gr  = 12      ! gridscale rain precipitation
INTEGER, PARAMETER ::   jps_gs  = 13      ! gridscale snow precipitation
INTEGER, PARAMETER ::   jps_gg  = 14      ! gridscale graupel precipitation
INTEGER, PARAMETER ::   jps_cp  = 15      ! total convective precipitation
INTEGER, PARAMETER ::   jps_gp  = 16      ! total gridscale precipitation

INTEGER, PARAMETER ::   jpss_lai   =  1   ! LAI
INTEGER, PARAMETER ::   jpss_plcov =  2   ! PLCOV
INTEGER, PARAMETER ::   jpss_xveg  =  3   ! vegetation factor
INTEGER, PARAMETER ::   jpss_sw    =  4   ! short-wave downward radiation
INTEGER, PARAMETER ::   jpss_lw    =  5   ! long-wave downward radiation
INTEGER, PARAMETER ::   jpss_sn    =  6   ! snow precipitation
INTEGER, PARAMETER ::   jpss_pr    =  7   ! rain precipitation
INTEGER, PARAMETER ::   jpss_t     =  8   ! temperature lowest model layer
INTEGER, PARAMETER ::   jpss_w     =  9   ! wind velocitiy calculated out of u and v
INTEGER, PARAMETER ::   jpss_qv    = 10   ! specific water vapor content lowest model layer
INTEGER, PARAMETER ::   jpss_ps    = 11   ! surface pressure 
INTEGER, PARAMETER ::   jpss_pa    = 12   ! air pressure lowest model layer

INTEGER(KIND=iintegers) :: &
  isec,            &
  isec_2wn,        &
  kinfo,           &
  il_var_id(16),   &  
  dimids(2),       &
  il_file_id,      &
  mype,            &
  ib,              &
  npes,            &
  ierror,          &
  status,          &
  ier,             &
  jn, jm, jo,      &
  i, j, k,         &
  ii, jj,          &
  ztl,             &
  ingp,            &
  izerror,         &
  izerrstat,       &
  kzdims(24)         ! vertical dimensions for boundary exchange

REAL(KIND=ireals)      :: &
  ztmp1  (ie,je),         &
  ztprec (ie,je),         &
  ztevap (ie,je),         &
  xveg   (ie,je),         &
  soiltyp_real (ie,je),   &
  sodown (ie,je),         &
  prec_s(ie,je),          &
  prec_r(ie,je),          &
  u_mitte(ie,je),         &
  v_mitte(ie,je),         &
  wind(ie,je),            &
  p_air(ie,je),           &
  zpart,                  &
  zlev,                   &
  integ3d,                &
  t_tot  (ie_tot,je_tot), &
  t_lm_tot  (ie_tot,je_tot), &
  t_cm_tot  (ie_tot,je_tot), &
  zps_tot (ie_tot,je_tot), &
  zlat   (ie,je),         &
  zlon   (ie,je),         &
  zp     (ie,je,ke_in)

CHARACTER (LEN=80)      :: &
  yzerrmsg

CHARACTER(LEN=29)       :: &
  fmt='(F20.2,1X,F20.2,1X,F20.2)'  ! format descriptor

CHARACTER (LEN=25)       :: yzroutine

REAL (KIND=ireals), POINTER :: &
    qv_now     (:,:,:)   => NULL()! , &     ! QV at nx

!------------------------------------------------------------------------------
!- End of header
!------------------------------------------------------------------------------

  izerror = 0_iintegers   
  yzerrmsg = ''
  yzroutine= 'send_fld'

 ! Retrieve the required microphysics tracers
  CALL trcr_get(izerror, idt_qv,  ptr_tlev = nnow, ptr = qv_now)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF

!----------------------------------------------------------------------------
! 1 coupling-specific handling of sent fields
!----------------------------------------------------------------------------

! Coupling only on PE with at least one unmasked grid point
IF ( lpe_cpl ) THEN

  isec     = nint( (ntstep - nstart ) * dt )
  isec_2wn = nint( (ntstep - nstart ) * dt - dt_cp ) ! start sending at 0 s

  IF ( debug_oasis > 15 .AND. my_cart_id == 0 ) THEN
    WRITE(nulout,*) 'in send_fld: ntstep=',ntstep,' isec=',isec,' isec_2wn=',isec_2wn !SW, 26.08.13
  ENDIF

!----------------------------------------------------------------------------
! 1.1 clmXX-specific handling of sent fields
!----------------------------------------------------------------------------

  IF ( ytype_lsm == 'clmXX' ) THEN

    IF( ssnd(jps_t)%laction )   CALL oas_cos_snd( jps_t, isec, t(:,:,ke,nnow), kinfo )
    IF( ssnd(jps_u)%laction )   CALL oas_cos_snd( jps_u, isec, u(:,:,ke,nnow), kinfo )
    IF( ssnd(jps_v)%laction )   CALL oas_cos_snd( jps_v, isec, v(:,:,ke,nnow), kinfo )
    IF( ssnd(jps_q)%laction )   CALL oas_cos_snd( jps_q, isec, qv_now(:,:,ke), kinfo)
    ! Send only total height
    ztmp1(:,:) = hhl(:,:,ke) - hsurf(:,:)
    IF( ssnd(jps_th)%laction )   CALL oas_cos_snd( jps_th, isec, ztmp1(:,:), kinfo )

    ! Send only total pressure
    ztmp1(:,:) = p0(:,:,ke) + pp(:,:,ke,nnow)
    IF( ssnd(jps_pr)%laction )   CALL oas_cos_snd( jps_pr, isec, ztmp1(:,:), kinfo )

    IF( ssnd(jps_rs)%laction )   CALL oas_cos_snd( jps_rs, isec, swdir_s(:,:), kinfo )
    IF( ssnd(jps_fs)%laction )   CALL oas_cos_snd( jps_fs, isec, swdifd_s(:,:), kinfo )
    IF( ssnd(jps_lw)%laction )   CALL oas_cos_snd( jps_lw, isec, lwd_s(:,:), kinfo )
    IF( ssnd(jps_cr)%laction )   CALL oas_cos_snd( jps_cr, isec, prr_con(:,:), kinfo )
    IF( ssnd(jps_cs)%laction )   CALL oas_cos_snd( jps_cs, isec, prs_con(:,:), kinfo )
    IF( ssnd(jps_gr)%laction )   CALL oas_cos_snd( jps_gr, isec, prr_gsp(:,:), kinfo )
    IF( ssnd(jps_gs)%laction )   CALL oas_cos_snd( jps_gs, isec, prs_gsp(:,:), kinfo )

    ! Send graupel if any, 0 otherwise
    ztmp1(:,:) = 0.
    IF ( itype_gscp == 4 ) ztmp1(:,:) = prg_gsp(:,:)
    IF( ssnd(jps_gg)%laction )   CALL oas_cos_snd( jps_gg, isec, ztmp1(:,:), kinfo )

    ! Send convective and gridscale total precipitations only
    ztmp1(:,:) = prr_con(:,:) + prs_con(:,:)
    IF( ssnd(jps_cp)%laction )   CALL oas_cos_snd( jps_cp, isec, ztmp1(:,:), kinfo )
    ztmp1(:,:) = prr_gsp(:,:) + prs_gsp(:,:)
    IF ( itype_gscp == 4 ) ztmp1(:,:) = ztmp1(:,:) + prg_gsp(:,:)
    IF( ssnd(jps_gp)%laction )   CALL oas_cos_snd( jps_gp, isec, ztmp1(:,:), kinfo )

  ENDIF ! ytype_lsm == 'clmXX'

!----------------------------------------------------------------------------
! 1.2 veg3d-specific handling of sent fields
!----------------------------------------------------------------------------

  IF ( ytype_lsm == 'veg3d' ) THEN

     CALL vegfaktor (istartpar, iendpar, jstartpar, jendpar, xveg)
     sodown(:,:)              = swdir_s(:,:) + swdifd_s(:,:)
     prec_s(:,:)              = prs_con(:,:) + prs_gsp(:,:)
     prec_r(:,:)              = prr_con(:,:) + prr_gsp(:,:)
     p_air(:,:)               = pp(:, :, ke, nnow) + p0(:, :, ke)

     DO j = jstartpar, jendpar
          jj  = MAX( j-1, 1 )
        DO i = istartpar, iendpar
           ii  = MAX( i-1, istartpar )
           u_mitte(i,j) = 0.5_ireals * ( u(i,j,ke,nnow) + u(ii,j,ke,nnow) )
           v_mitte(i,j) = 0.5_ireals * ( v(i,j,ke,nnow) + v(i,jj,ke,nnow) )
           wind(i,j)  = SQRT(u_mitte(i,j) ** 2.0_ireals + v_mitte(i,j) ** 2.0_ireals)
           wind(i,j)  = MAX(wind(i,j), vel_min)
        ENDDO
     ENDDO

     IF( ssnd(jpss_lai)%laction )   CALL oas_cos_snd( jpss_lai, isec, lai(:,:), kinfo )
     IF( ssnd(jpss_plcov)%laction )   CALL oas_cos_snd( jpss_plcov, isec, plcov(:,:), kinfo )
     IF( ssnd(jpss_xveg)%laction )   CALL oas_cos_snd( jpss_xveg, isec, xveg(:,:), kinfo )
     IF( ssnd(jpss_sw)%laction )   CALL oas_cos_snd( jpss_sw, isec, sodown(:,:), kinfo )
     IF( ssnd(jpss_lw)%laction )   CALL oas_cos_snd( jpss_lw, isec, lwd_s(:,:), kinfo )
     IF( ssnd(jpss_sn)%laction )   CALL oas_cos_snd( jpss_sn, isec, prec_s(:,:), kinfo )
     IF( ssnd(jpss_pr)%laction )   CALL oas_cos_snd( jpss_pr, isec, prec_r(:,:), kinfo )
     IF( ssnd(jpss_t)%laction )   CALL oas_cos_snd( jpss_t, isec, t(:,:,ke,nnow), kinfo )
     IF( ssnd(jpss_w)%laction )   CALL oas_cos_snd( jpss_w, isec, wind(:,:), kinfo )
     IF( ssnd(jpss_qv)%laction )   CALL oas_cos_snd( jpss_qv, isec, qv_now(:,:,ke), kinfo )
     IF( ssnd(jpss_ps)%laction )   CALL oas_cos_snd( jpss_ps, isec, ps(:,:,nnow), kinfo )
     IF( ssnd(jpss_pa)%laction )   CALL oas_cos_snd( jpss_pa, isec, p_air(:,:), kinfo )
     
  ENDIF ! ytype_lsm == 'veg3d'


!------------------------------------------------------------------------------
! 1.3 echam-specific handling of sent fields
!------------------------------------------------------------------------------

  IF ( (ytype_2wn == 'echam') .AND. l2wn ) THEN
    IF ( ntstep > nstart .AND. MOD(isec_2wn,dt_cp) == 0 ) THEN   ! start at time=0s

      ztl = nnew
      IF ( lreal_2wn ) THEN
        t_lm (:,:,:) = 0.0_ireals 
        u_lm (:,:,:) = 0.0_ireals
        v_lm (:,:,:) = 0.0_ireals
        qv_lm(:,:,:) = 0.0_ireals
        qc_lm(:,:,:) = 0.0_ireals
        qi_lm(:,:,:) = 0.0_ireals

        ! copy current, i.e. at the end of a time step, (nnew) fields to working arrays
        t_lm (:,:,1:ke) = t (:,:,:,ztl) 
        u_lm (:,:,1:ke) = u (:,:,:,ztl)
        v_lm (:,:,1:ke) = v (:,:,:,ztl)
!!!! ELD WARNING this has to be changed to nnew !!!!!
        qv_lm(:,:,1:ke) = qv_now(:,:,:)
!        qv_lm(:,:,1:ke) = qv(:,:,:,ztl)
        qc_lm(:,:,1:ke) = qc(:,:,:,ztl)
        qi_lm(:,:,1:ke) = qi(:,:,:,ztl)
        pp_lm(:,:,1:ke) = pp(:,:,:,ztl)
        ps_lm(:,:)      = ps(:,:,ztl) ! ps(nnew) is calculated in "near_surface" 
                                      ! at the end of lmorg
      ENDIF ! lreal_2wn

      IF (ltime) CALL get_timings (i_cpl_put, ntstep, dt, ierror)

      !------------------------------------------------------------------------
      ! 1.3.1 vertical interpolation to coarse-model levels
      !------------------------------------------------------------------------

      ! some debug output
      CALL gather_field (ps_lm(:,:),ie,je,zps_tot,ie_tot,je_tot,0,nerror)
      CALL write_2d_ps (zps_tot,'PS','ps_snd_1.nc','write',5,nrecord,1)
      !DO k = 1, ke_in
      !  CALL gather_field (t_lm(:,:,k),ie,je,zps_tot,ie_tot,je_tot,0,nerror)
      !  CALL write_2d_ps (zps_tot,'T','ps_snd_1.nc','write',5,nrecord,k)
      !ENDDO

      CALL org_vert_interp_2wn ('r2g')

!IF ( debug_oasis > 15 .AND. my_cart_id==0) WRITE(nulout,*) 'Mean temperature difference mean(t_lm(k))-mean(t_cm(k)) in send_fld'
!DO k = ke_out, ke_in
!  CALL gather_field(t_cm(:,:,k),ie,je,t_cm_tot,ie_tot,je_tot,0,nerror)
!  CALL gather_field(t_lm(:,:,k),ie,je,t_lm_tot,ie_tot,je_tot,0,nerror)
!  IF (my_cart_id==0) THEN
!    WRITE(nulout,'(I2,1X,F10.5)') k, (SUM(t_lm_tot)/(MAX(1,SIZE(t_lm_tot))))-(SUM(t_cm_tot)/(MAX(1,SIZE(t_cm_tot))))
!    CALL flush(nulout)
!  ENDIF
!ENDDO

      IF (ltime) CALL get_timings (i_cpl_vert_interp, ntstep, dt, ierror)

      !----------------------------------------------------------------------------
      ! 1.3.2 Convert the U and V velocities from the rotated to the real
      !       geographical coordinate system (see int2lm: SUBROUTINE interpol_coarse_uv)
      !----------------------------------------------------------------------------

      IF ( ltransform_uv ) THEN

        !u_cm(:,:,:) = 10.0_ireals
        !v_cm(:,:,:) = 5.0_ireals
        !u_lm(:,:,:) = 10.0_ireals
        !v_lm(:,:,:) = 5.0_ireals
        !CALL convert_uv2uvrot

        CALL convert_uvrot2uv

        IF (ltime) CALL get_timings (i_cpl_add_comp, ntstep, dt, izerror)
      ENDIF

      !------------------------------------------------------------------------
      ! 1.3.3 apply vertical relaxation of COSMO to ECHAM fields below the 
      !       Rayleigh damping (rdheight); the relaxation function is a cosine 
      !       function; u and v velocities are all in the true geographical 
      !       lat/lon coordinate system
      !------------------------------------------------------------------------

      IF ( lmerge_profiles ) THEN
        zlev = -1.0_ireals
        DO k = ke_out, ke_out+nlev_merge
          zlev = zlev + 1.0_ireals
          zpart = COS((pi/2.0_ireals) * (zlev/REAL(nlev_merge,ireals)))**2
          !IF ( my_cart_id == 0 ) PRINT *, 'actual level k=',k
          !IF ( my_cart_id == 0 ) PRINT *, 'zlev=',zlev
          !IF ( my_cart_id == 0 ) PRINT *, 'zpart=',zpart
          IF ( zpart < 0.0_ireals ) THEN
            zpart = 0.0_ireals
            !IF ( my_cart_id == 0 ) PRINT *, 'zpart=',zpart
          ENDIF
          t_lm  (:,:,k) = t_lm  (:,:,k)*(1.0_ireals-zpart) + t_cm  (:,:,k)*zpart
          u_lm  (:,:,k) = u_lm  (:,:,k)*(1.0_ireals-zpart) + u_cm  (:,:,k)*zpart
          v_lm  (:,:,k) = v_lm  (:,:,k)*(1.0_ireals-zpart) + v_cm  (:,:,k)*zpart
          qv_lm (:,:,k) = qv_lm (:,:,k)*(1.0_ireals-zpart) + qv_cm (:,:,k)*zpart
          qc_lm (:,:,k) = qc_lm (:,:,k)*(1.0_ireals-zpart) + qc_cm (:,:,k)*zpart
          qi_lm (:,:,k) = qi_lm (:,:,k)*(1.0_ireals-zpart) + qi_cm (:,:,k)*zpart
        ENDDO
      ENDIF

      !----------------------------------------------------------------------------
      ! 1.3.4 Overwrite original ECHAM fields in parts of all levels by the 
      !       CCLM solution
      !----------------------------------------------------------------------------

      ! Only the levels ke_out to ke_in will be sent back to ECHAM
      IF ( loverwrite_cm ) THEN
        t_cm  (:,:,ke_out:ke_in) = t_lm  (:,:,ke_out:ke_in)
        u_cm  (:,:,ke_out:ke_in) = u_lm  (:,:,ke_out:ke_in)
        v_cm  (:,:,ke_out:ke_in) = v_lm  (:,:,ke_out:ke_in)
        qv_cm (:,:,ke_out:ke_in) = qv_lm (:,:,ke_out:ke_in)
        qc_cm (:,:,ke_out:ke_in) = qc_lm (:,:,ke_out:ke_in)
        qi_cm (:,:,ke_out:ke_in) = qi_lm (:,:,ke_out:ke_in)
        !ps_cm (:,:)              = ps_lm (:,:) !SW 18.07.14
        ps_cm (:,:)              = ps_gl (:,:)
      ENDIF

      IF ( ke_out > 1 ) THEN
        t_lm  (:,:,1:ke_out-1) = t_cm  (:,:,1:ke_out-1)
        u_lm  (:,:,1:ke_out-1) = u_cm  (:,:,1:ke_out-1)
        v_lm  (:,:,1:ke_out-1) = v_cm  (:,:,1:ke_out-1)
        qv_lm (:,:,1:ke_out-1) = qv_cm (:,:,1:ke_out-1)
        qc_lm (:,:,1:ke_out-1) = qc_cm (:,:,1:ke_out-1)
        qi_lm (:,:,1:ke_out-1) = qi_cm (:,:,1:ke_out-1)
      ENDIF

      IF (ltime) CALL get_timings (i_cpl_add_comp, ntstep, dt, ierror)

      !------------------------------------------------------------------------
      ! 1.3.5 Calculation of horizontal gradients (based on *_cm fields)
      !------------------------------------------------------------------------

      IF ( lhoriz_grad ) THEN

        ! exchange boundaries of fields for which horizontal gradients will be
        ! computed
        IF ( num_compute > 1 ) THEN
          kzdims(1:24) = (/ ke_in,ke_in,ke_in,ke_in,ke_in,ke_in,1,         &
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /)
     !! naveed --    
     ! CALL exchg_boundaries                                            &
     !       (0, sendbuf, isendbuflen, imp_reals, icomm_cart, ie, je,       &
     !       kzdims, jstartpar, jendpar, 3, nboundlines, my_cart_neigh,     &
     !       20000+ntstep, .FALSE., ncomm_type, izerror, yzerrmsg,          &
     !       t_cm(:,:,:),   u_cm(:,:,:),  v_cm(:,:,:),                      &
     !       qv_cm(:,:,:), qc_cm(:,:,:), qi_cm(:,:,:), ps_cm(:,:))
        ENDIF

        CALL gradients_2wn (ke_out, ke_in)

        IF (ltime) CALL get_timings (i_cpl_grad, ntstep, dt, ierror)

      ENDIF

      !------------------------------------------------------------------------
      ! 1.3.6 Calculation of a special mask
      !------------------------------------------------------------------------

      ! This section should only be activiated if new masks need to be created (see the 
      ! documentation on the namelist parameters and the description of how to create new
      ! masks).
      IF ( lcreate_mask ) THEN
        ! 1) Create a mask for ECHAM to COSMO grid:
        ps_cm(:,:) = 1.0_ireals
        ! 2) Create a mask for COSMO to ECHAM grid:
        CALL gather_field(ps_cm(:,:),ie,je,t_tot,ie_tot,je_tot,0,nerror)
        IF ( my_cart_id == 0 ) THEN
          t_tot (1:ngp_mask_frame,:) = 0.0_ireals               !left edge
          t_tot (ie_tot-ngp_mask_frame+1:ie_tot,:) = 0.0_ireals !right edge
          t_tot (:,1:ngp_mask_frame) = 0.0_ireals               !lower edge
          t_tot (:,je_tot-ngp_mask_frame+1:je_tot) = 0.0_ireals !upper edge
          ! Print some output:
          WRITE(nulout,*) 'Create mask'
          WRITE(nulout,*) 'Number of grid point rows/columns set to zero: ',ngp_mask_frame
          WRITE(nulout,*) 'Left edge : from 1 to ',ngp_mask_frame
          WRITE(nulout,*) 'Right edge: from ',ie_tot-ngp_mask_frame+1,' to ',ie_tot
          WRITE(nulout,*) 'Lower edge: from 1 to ',ngp_mask_frame
          WRITE(nulout,*) 'Upper edge: from ',je_tot-ngp_mask_frame+1,' to ',je_tot
          CALL flush(nulout)
        ENDIF
!! naveed --       
! CALL distribute_field (t_tot(:,:),ie_tot,je_tot,ps_cm(:,:),ie,je,nerror)
      ENDIF

      !----------------------------------------------------------------------------
      ! 1.3.7 Calculate volume integrals of all humidity variables and write to file
      !----------------------------------------------------------------------------

      IF ( lmeasure_mass ) THEN
        integ3d = integral_3d_total (qv_lm(:,:,:),3)
        IF (my_cart_id==0) THEN
          OPEN (numassqv,FILE='total_mass_qv.txt',FORM='FORMATTED',ACTION='WRITE',&
            POSITION='APPEND',STATUS='OLD',IOSTAT=izerrstat)
          WRITE (numassqv,'(A,1X,I10,1X,I10,1X,F30.5)') 'snd',ntstep,isec_2wn,integ3d
          CLOSE (numassqv,STATUS='KEEP',IOSTAT=izerrstat)
        ENDIF
        integ3d = integral_3d_total (qc_lm(:,:,:),3)
        IF (my_cart_id==0) THEN
          OPEN (numassqc,FILE='total_mass_qc.txt',FORM='FORMATTED',ACTION='WRITE',&
            POSITION='APPEND',STATUS='OLD',IOSTAT=izerrstat)
          WRITE (numassqc,'(A,1X,I10,1X,I10,1X,F30.5)') 'snd',ntstep,isec_2wn,integ3d
          CLOSE (numassqc,STATUS='KEEP',IOSTAT=izerrstat)
        ENDIF
        integ3d = integral_3d_total (qi_lm(:,:,:),3)
        IF (my_cart_id==0) THEN
          OPEN (numassqi,FILE='total_mass_qi.txt',FORM='FORMATTED',ACTION='WRITE',&
            POSITION='APPEND',STATUS='OLD',IOSTAT=izerrstat)
          WRITE (numassqi,'(A,1X,I10,1X,I10,1X,F30.5)') 'snd',ntstep,isec_2wn,integ3d
          CLOSE (numassqi,STATUS='KEEP',IOSTAT=izerrstat)
        ENDIF
        integ3d = integral_3d_total (qr(:,:,:,ztl),3)
        IF (my_cart_id==0) THEN
          OPEN (numassqr,FILE='total_mass_qr.txt',FORM='FORMATTED',ACTION='WRITE',&
            POSITION='APPEND',STATUS='OLD',IOSTAT=izerrstat)
          WRITE (numassqr,'(I10,1X,I10,1X,F30.5)') ntstep,isec_2wn,integ3d
          CLOSE (numassqr,STATUS='KEEP',IOSTAT=izerrstat)
        ENDIF
        integ3d = integral_3d_total (qs(:,:,:,ztl),3)
        IF (my_cart_id==0) THEN
          OPEN (numassqs,FILE='total_mass_qs.txt',FORM='FORMATTED',ACTION='WRITE',&
            POSITION='APPEND',STATUS='OLD',IOSTAT=izerrstat)
          WRITE (numassqs,'(I10,1X,I10,1X,F30.5)') ntstep,isec_2wn,integ3d
          CLOSE (numassqs,STATUS='KEEP',IOSTAT=izerrstat)
        ENDIF

        IF (ltime) CALL get_timings (i_cpl_add_comp, ntstep, dt, ierror)
      ENDIF

      !------------------------------------------------------------------------
      ! 1.3.8 send the fields
      !------------------------------------------------------------------------

      IF ( nrecord < 800 ) CALL nc_diagnose ('add_var', 'nach_vi.nc', ke_in, 2, nrecord)
      IF ( nrecord < 800 ) CALL nc_diagnose ('add_var_2d', 'nach_vi.nc', ke_in, 2, nrecord)

      jn = nfld_snd_lsm 
      DO jo = 1, nfld_snd_2wn
        DO jm = ke_out, nlev_snd_2wn(jo)+ke_out-1
          jn = jn + 1
          IF( ssnd(jn)%laction ) THEN

            ! distinguish between first- and second-order conservative interpolation
            ! via OASIS3-MCT
            IF ( lhoriz_grad .AND. lconserv_2o ) THEN ! second-order
              IF ( jo == 1 ) CALL oas_cos_snd(jn, isec_2wn, ps_cm(:,:), kinfo,    &
                                              dpsdx(:,:), dpsdy(:,:))
              !IF ( jo == 2 ) CALL oas_cos_snd(jn, isec_2wn, beta_cm(:,:), kinfo)
              IF ( jo == 2 .AND. isec_2wn == 0 ) CALL oas_cos_snd(jn, isec_2wn, beta_cm(:,:), kinfo) ! BET NUR ZURUECK
              IF ( jo == 3 ) CALL oas_cos_snd(jn, isec_2wn, t_cm(:,:,jm), kinfo,  &
                                              dtdx(:,:,jm), dtdy(:,:,jm))
              IF ( jo == 4 ) CALL oas_cos_snd(jn, isec_2wn, u_cm(:,:,jm), kinfo,  &
                                              dudx(:,:,jm), dudy(:,:,jm))
              IF ( jo == 5 ) CALL oas_cos_snd(jn, isec_2wn, v_cm(:,:,jm), kinfo,  &
                                              dvdx(:,:,jm), dvdy(:,:,jm))
              IF ( jo == 6 ) CALL oas_cos_snd(jn, isec_2wn, qv_cm(:,:,jm), kinfo, &
                                              dqvdx(:,:,jm), dqvdy(:,:,jm))
              IF ( jo == 7 ) CALL oas_cos_snd(jn, isec_2wn, qc_cm(:,:,jm), kinfo, &
                                              dqcdx(:,:,jm), dqcdy(:,:,jm))
              IF ( jo == 8 ) CALL oas_cos_snd(jn, isec_2wn, qi_cm(:,:,jm), kinfo, &
                                              dqidx(:,:,jm), dqidy(:,:,jm))
            ELSE ! first-order (if specified in the namcouple file)
              IF ( jo == 1 ) CALL oas_cos_snd(jn, isec_2wn, ps_cm(:,:), kinfo )
              !IF ( jo == 2 ) CALL oas_cos_snd(jn, isec_2wn, beta_cm(:,:), kinfo )
              IF ( jo == 2 .AND. isec_2wn == 0 ) CALL oas_cos_snd(jn, isec_2wn, beta_cm(:,:), kinfo ) ! BET NUR ZURUECK
              IF ( jo == 3 ) CALL oas_cos_snd(jn, isec_2wn, t_cm(:,:,jm), kinfo )
              IF ( jo == 4 ) CALL oas_cos_snd(jn, isec_2wn, u_cm(:,:,jm), kinfo )
              IF ( jo == 5 ) CALL oas_cos_snd(jn, isec_2wn, v_cm(:,:,jm), kinfo )
              IF ( jo == 6 ) CALL oas_cos_snd(jn, isec_2wn, qv_cm(:,:,jm), kinfo )
              IF ( jo == 7 ) CALL oas_cos_snd(jn, isec_2wn, qc_cm(:,:,jm), kinfo )
              IF ( jo == 8 ) CALL oas_cos_snd(jn, isec_2wn, qi_cm(:,:,jm), kinfo )
            ENDIF

            ! send gradients themselves (these are interpolated first-order only) 
            IF ( lhoriz_grad ) THEN
              IF ( jo ==  9 ) CALL oas_cos_snd( jn, isec_2wn, dtdx (:,:,jm), kinfo )
              IF ( jo == 10 ) CALL oas_cos_snd( jn, isec_2wn, dtdy (:,:,jm), kinfo )
              IF ( jo == 11 ) CALL oas_cos_snd( jn, isec_2wn, dpsdx(:,:   ), kinfo )
              IF ( jo == 12 ) CALL oas_cos_snd( jn, isec_2wn, dpsdy(:,:   ), kinfo )
              IF ( jo == 13 ) CALL oas_cos_snd( jn, isec_2wn, dudx (:,:,jm), kinfo )
              IF ( jo == 14 ) CALL oas_cos_snd( jn, isec_2wn, dudy (:,:,jm), kinfo )
              IF ( jo == 15 ) CALL oas_cos_snd( jn, isec_2wn, dvdx (:,:,jm), kinfo )
              IF ( jo == 16 ) CALL oas_cos_snd( jn, isec_2wn, dvdy (:,:,jm), kinfo )
            ENDIF
          ENDIF
        ENDDO ! jm / level
      ENDDO ! jo / field

      IF (  debug_oasis > 15 .AND. my_cart_id == 0 ) THEN
        WRITE(nulout,*) 'Sent fields at ntstep=',ntstep,' isec_2wn=',isec_2wn
        CALL flush(nulout)
      ENDIF

    ENDIF ! MOD
  ENDIF ! ytype_2wn .AND. l2wn

!----------------------------------------------------------------------------
! 1.4 nemoO-nemoI-nemoD and romsO-specific handling of sent fields
!----------------------------------------------------------------------------

  IF (  ytype_oce /= 'nooce' ) THEN
    ! Index of the first ocean coupled field
    jn = nfld_snd_lsm + SUM ( nlev_snd_2wn )  
  ENDIF

  IF (ytype_oce == 'romsO') THEN
     jn = jn + 1
     IF( ssnd(jn)%laction ) CALL oas_cos_snd (jn, isec, umfl_s, kinfo )
     jn = jn + 1
     IF( ssnd(jn)%laction ) CALL oas_cos_snd (jn, isec, vmfl_s, kinfo )
     jn = jn + 1
     IF( ssnd(jn)%laction ) THEN
       ztmp1 = sobs + thbs + lhfl_s + shfl_s
       CALL oas_cos_snd (jn, isec, ztmp1, kinfo )
     ENDIF
     jn = jn + 1
     IF( ssnd(jn)%laction ) THEN
       ztmp1 = swdir_s +swdifd_s
       CALL oas_cos_snd (jn, isec,ztmp1, kinfo )
     ENDIF
     jn = jn + 1
     IF( ssnd(jn)%laction ) THEN
       ztevap = -lhfl_s/lh_v
       ztmp1 = ztevap - prr_con - prs_con - prr_gsp - prs_gsp
       CALL oas_cos_snd (jn, isec, ztmp1, kinfo )
    ENDIF
  ENDIF

  IF ( ytype_oce == 'nemoO' .OR. ytype_oce == 'nemoD' ) THEN
    jn = jn + 1
    IF( ssnd(jn)%laction ) THEN
      ztprec = prr_con + prs_con + prr_gsp + prs_gsp
      IF (itype_gscp == 4 ) ztprec = ztprec + prg_gsp 
        !   the surface freshwater flux is calculated, should be positiv downward
        !   for nemo, so the sign is changed
      ztevap = lhfl_s/lh_v
      ztmp1 = -ztevap - ztprec
      CALL oas_cos_snd (jn, isec, ztmp1, kinfo )
    ENDIF
    jn = jn + 1
    IF( ssnd(jn)%laction ) CALL oas_cos_snd (jn, isec, umfl_s, kinfo )
    jn = jn + 1
    IF( ssnd(jn)%laction ) CALL oas_cos_snd (jn, isec, vmfl_s, kinfo ) 
    jn = jn + 1
    IF( ssnd(jn)%laction ) CALL oas_cos_snd (jn, isec, sobs, kinfo )
    jn = jn + 1
    IF( ssnd(jn)%laction ) THEN
      ztmp1 = thbs + lhfl_s + shfl_s
      CALL oas_cos_snd (jn, isec, ztmp1, kinfo )
    ENDIF
  ENDIF !! ytype_oce  

  IF ( ytype_oce == 'nemoI' .OR. ytype_oce == 'nemoD' ) THEN 
    jn = jn + 1
    IF( ssnd(jn)%laction ) THEN
      ztprec = prr_con + prs_con + prr_gsp + prs_gsp
      IF (itype_gscp == 4 ) ztprec = ztprec + prg_gsp
      !   the surface freshwater flux is calculated, should be positiv downward
      !   for nemo, so the sign is changed
      ztevap = lhfl_s/lh_v
      ztmp1 = -ztevap - ztprec
      CALL oas_cos_snd (jn, isec, ztmp1, kinfo )
    ENDIF
    jn = jn + 1 
    IF( ssnd(jn)%laction ) CALL oas_cos_snd (jn, isec, umfl_s, kinfo )
    jn = jn + 1
    IF( ssnd(jn)%laction ) CALL oas_cos_snd (jn, isec, vmfl_s, kinfo )
    jn = jn + 1
    IF( ssnd(jn)%laction ) CALL oas_cos_snd (jn, isec, sobs, kinfo )
    jn = jn + 1
    IF( ssnd(jn)%laction ) THEN
      ztmp1 = thbs + lhfl_s + shfl_s
      CALL oas_cos_snd (jn, isec, ztmp1, kinfo )
    ENDIF
    jn = jn + 1 
    IF( ssnd(jn)%laction ) CALL oas_cos_snd (jn, isec, ps, kinfo ) 
  ENDIF ! ytype_oce

!----------------------------------------------------------------------------
! 1.5 TRIMNP+CICE-specific handling of sent fields
!----------------------------------------------------------------------------
! uoi HaHo 2013-09-04 {
  IF ( ytype_oce == 'trimI' ) THEN

    jn = jn + 1
    IF( ssnd(jn)%laction ) CALL oas_cos_snd (jn, isec, pmsl1, kinfo )  ! 1) PMSLtTB [Pa]
    jn = jn + 1
    IF( ssnd(jn)%laction ) CALL oas_cos_snd (jn, isec, u_10m, kinfo )  ! 2) U_10MtTB [m/s]
    jn = jn + 1
    IF( ssnd(jn)%laction ) CALL oas_cos_snd (jn, isec, v_10m, kinfo )  ! 3) V_10MtTB [m/s]
    jn = jn + 1
    IF( ssnd(jn)%laction ) THEN
      ztmp1 = prr_con + prs_con + prr_gsp + prs_gsp
      IF (itype_gscp == 4 ) ztmp1 = ztmp1 + prg_gsp
      CALL oas_cos_snd (jn, isec, ztmp1, kinfo )                       ! 4) PR_tTB [kg/m2/s]
    ENDIF
    jn = jn + 1
    IF( ssnd(jn)%laction ) CALL oas_cos_snd (jn, isec, lhfl_s, kinfo ) ! 5) LHFL_StTB [W/m2]
    jn = jn + 1
    IF( ssnd(jn)%laction ) CALL oas_cos_snd (jn, isec, sobs, kinfo )   ! 6) SOBS_RADtTB [W/m2]
    jn = jn + 1
    IF( ssnd(jn)%laction ) THEN
      ztmp1 = sobs + lwd_s + shfl_s + lhfl_s
      CALL oas_cos_snd (jn, isec, ztmp1, kinfo )                       ! 7) HFLtTB [W/m2]
    ENDIF
!!!!! use state vars from CCLM to calculate HFLX in TRIM {
    jn = jn + 1
    IF( ssnd(jn)%laction ) CALL oas_cos_snd (jn, isec, t_2m, kinfo )   ! 8) T_2MtTB [K]

    jn = jn + 1
    IF( ssnd(jn)%laction ) CALL oas_cos_snd (jn, isec, rh_2m, kinfo )  ! 9) RH_2MtTB [%]

    jn = jn + 1
    IF( ssnd(jn)%laction ) CALL oas_cos_snd (jn, isec, clct, kinfo )   ! 10) CLCTtTB [1]
!!!!! use state vars from CCLM to calculate HFLX in TRIM }

    jn = jn + 1
    IF( ssnd(jn)%laction ) THEN
      ztmp1 (nldi:nlei,nldj:nlej) = 0.5* (hhl(nldi:nlei,nldj:nlej,ke) &
                                        + hhl(nldi:nlei,nldj:nlej,ke+1))
      CALL oas_cos_snd (jn, isec, ztmp1, kinfo )                       ! 11) HML_KEtTB [m]
    ENDIF
    jn = jn + 1
    IF( ssnd(jn)%laction ) THEN
!      ztmp1(nldi:nlei,nldj:nlej) = qv (nldi:nlei,nldj:nlej,ke,nnew)
!!!! ELD WARNING change to nnew !!!!!
      ztmp1(nldi:nlei,nldj:nlej) = qv_now (nldi:nlei,nldj:nlej,ke)
      CALL oas_cos_snd (jn, isec, ztmp1, kinfo )                       ! 12) QV_KEtTB [kg/kg]
    ENDIF
    jn = jn + 1
    IF( ssnd(jn)%laction ) CALL oas_cos_snd (jn, isec, rhoA, kinfo )   ! 13) RHOA_KEtTB [kg/m3]
    jn = jn + 1
    IF( ssnd(jn)%laction ) CALL oas_cos_snd (jn, isec, Tpot, kinfo )   ! 14) Tpot_KEtTB [K]
    jn = jn + 1
    IF( ssnd(jn)%laction ) THEN
      ztmp1(nldi:nlei,nldj:nlej) = t (nldi:nlei,nldj:nlej,ke,nnew)     ! 15) T_KEtTB [K]
      CALL oas_cos_snd (jn, isec, ztmp1, kinfo )
    ENDIF
    jn = jn + 1
    IF( ssnd(jn)%laction ) THEN
      ztmp1 = prr_con + prr_gsp   ! kg/m2/s = mm/s
      CALL oas_cos_snd (jn, isec, ztmp1, kinfo )                       ! 16) PRRtTB [kg/m2/s]
    ENDIF
    jn = jn + 1
    IF( ssnd(jn)%laction ) THEN
      ztmp1 = prs_con + prs_gsp   ! kg/m2/s = mm/s
      CALL oas_cos_snd (jn, isec, ztmp1, kinfo )                       ! 17) PRStTB [kg/m2/s]
    ENDIF
    jn = jn + 1
    IF( ssnd(jn)%laction ) CALL oas_cos_snd (jn, isec, lwd_s, kinfo )  ! 18) LWD_StTB [W/m2]

  ENDIF ! ytype_oce == 'trimI'
! uoi HaHo 2013-09-04 }

!----------------------------------------------------------------------------
! 1.6 MOM5-specific handling of sent fields
!----------------------------------------------------------------------------

  IF ( ytype_oce == 'momBS' ) THEN  ! sandra 02022018
    jn = jn + 1
    IF( ssnd(jn)%laction ) THEN                                         ! precipitation 
      ztprec = prr_con + prr_gsp 
      IF (itype_gscp == 4 ) ztprec = ztprec + prg_gsp
      !   the surface freshwater flux is calculated, should be positiv downward
      !   for MOM5, so the sign is changed
      !ztevap = lhfl_s/lh_v
      !ztmp1 = -ztevap - ztprec
      CALL oas_cos_snd (jn, isec, -ztprec, kinfo )
    ENDIF
    jn = jn + 1
    IF( ssnd(jn)%laction ) THEN                                         ! evaporation
      ztevap = lhfl_s/lh_v
      CALL oas_cos_snd (jn, isec, -ztevap, kinfo )
    ENDIF
    jn = jn + 1
    IF( ssnd(jn)%laction ) THEN                                         ! snow
      ztmp1 = prs_con + prs_gsp
      CALL oas_cos_snd (jn, isec, -ztmp1, kinfo )
    ENDIF
    jn = jn + 1
    IF( ssnd(jn)%laction ) CALL oas_cos_snd (jn, isec, ps, kinfo )      ! pressure at surface
    jn = jn + 1
    IF( ssnd(jn)%laction ) CALL oas_cos_snd (jn, isec, t_2m, kinfo )    ! 2m temerature
    jn = jn + 1
    IF( ssnd(jn)%laction ) CALL oas_cos_snd (jn, isec, qv_2m, kinfo )   ! specific humidity
    jn = jn + 1
    IF( ssnd(jn)%laction ) CALL oas_cos_snd (jn, isec, u_10m, kinfo )   ! 10m u wind 
    jn = jn + 1
    IF( ssnd(jn)%laction ) CALL oas_cos_snd (jn, isec, v_10m, kinfo )   ! 10m v wind
    jn = jn + 1
    IF( ssnd(jn)%laction ) CALL oas_cos_snd (jn, isec, thbs, kinfo )    ! net long wave radiation
    jn = jn + 1
    IF( ssnd(jn)%laction ) CALL oas_cos_snd (jn, isec, sobs, kinfo )    ! net short wave radiation
    jn = jn + 1
    IF( ssnd(jn)%laction ) CALL oas_cos_snd (jn, isec, lhfl_s, kinfo )  ! latent heat flux
    jn = jn + 1
    IF( ssnd(jn)%laction ) CALL oas_cos_snd (jn, isec, shfl_s, kinfo )  ! sensible heat flux

  ENDIF ! ytype_oce  ! sandra 02022018



ENDIF ! lpe_cpl

IF (ltime) CALL get_timings (i_cpl_put, ntstep, dt, ierror)

END SUBROUTINE send_fld


!!---------------------------------------------------------------------
!!              ***  ROUTINE vegfaktor  ***
!!
!! ** Purpose : compute a vegetation factor depending on the circle
!!              of latitude and the season (adopted from INT2LM)
!!
!!----------------------------------------------------------------------

SUBROUTINE vegfaktor (istartpar, iendpar, jstartpar, jendpar, xveg)

USE data_parameters,   ONLY: &
  ireals,                    &
  iintegers
USE data_modelconfig,  ONLY: &
 ie,                         &
 je,                         &
 pollon,                     &
 pollat,                     &
 polgam
USE data_runcontrol,   ONLY: &
 nzjulianday
USE data_constants,    ONLY: &
 g,                          &
 pi
USE data_fields,       ONLY :&
 hsurf,                      &
 rlon,                       &
 rlat
USE utilities,         ONLY: &
 phirot2phi

  IMPLICIT NONE

  INTEGER, INTENT(in)   ::   &
 istartpar,                  &
 iendpar,                    &
 jstartpar,                  &
 jendpar
  REAL (KIND=ireals), DIMENSION(ie, je), INTENT(out)  :: &
   xveg                        !vegetation factor

  ! local variables
  INTEGER (KIND=iintegers)  :: &
  i,                           &
  j,                           &
 nzday
!day of the year
  REAL (KIND=ireals)                                              ::         &!
 zlat                , &
 zhred               , & !height reduction factor
 zbvp                , & !begin of the vegetation period (in day of theyear) as a
                               !function of the latitude
 zdvp                    !length of the vegetation period (in days) as a
                               !function of the latitude
!------------------------------------------------------------------------------
! -End of header
!======================================================================================
!------------------------------------------------------------------------------
! Begin Subroutine vegfaktor
!------------------------------------------------------------------------------

     DO i = istartpar, iendpar
        DO j = jstartpar, jendpar
        !use rotated coordinates
           zlat = phirot2phi (rlat(i,j), rlon(i,j), pollat, pollon, polgam)
           zhred = EXP(-5.0E-9_ireals * g * hsurf(i,j)**2)
           zbvp  = MAX(1.0_ireals, 3.0_ireals*(ABS(zlat) - 20.0_ireals))
           zdvp  = MIN(365.0_ireals,                                           &
                   345.0_ireals - 4.5_ireals*(ABS(zlat) - 20.0_ireals))

           IF ( zlat < 0.0_ireals ) THEN
              nzday = MOD(nzjulianday + 180_iintegers, 365_iintegers)
           ELSE
              nzday = nzjulianday
           ENDIF

           IF ( zdvp >= 345.0_ireals ) THEN
              xveg(i,j) = zhred
           ELSE IF ( nzday < NINT(zbvp) ) THEN
              xveg(i,j) = 0.0_ireals
           ELSE IF ( nzday > NINT(zbvp+zdvp) ) THEN
              xveg(i,j) = 0.0_ireals
           ELSE
              xveg(i,j) = MAX(0.0_ireals, MIN ( 1.0_ireals, 1.12_ireals*      &
                          SIN( pi*MAX(0.0_ireals, (nzday-zbvp))/zdvp ) ) )*zhred
           ENDIF

        END DO
     END DO

END SUBROUTINE vegfaktor
#endif
