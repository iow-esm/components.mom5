!
! Macros used from fms_platform.h:
! __ALLOCATABLE maps to either POINTER  or _ALLOCATABLE
! _NULL        maps to either =>NULL() or "nothing"
#include <fms_platform.h>

module ocean_drifters_mod

!!#define _DEBUG
!
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Matt Harrison 
!</CONTACT>
!<CONTACT EMAIL="torsten.seifert@io-warnemuende.de"> Torsten Seifert 
!</CONTACT>
!<CONTACT EMAIL="klaus-ketelsen@t-online.de"> Klaus Ketelsen 
!</CONTACT>
!
!<OVERVIEW>
!  Advect USER supplied drifters using the shared/drifters package.
!</OVERVIEW>
!
!<DESCRIPTION>
!  Advect USER supplied drifters using the shared/drifters package.
!  IOW version 3.0 from 2014/11/28
!  code changes:
!  synchronize drifter velocities to avoid dependence on domain decomposition
!  define x1d,y1d on xu,yu grid to find t-cells with u,v on its edges for B-grid
!  unique drifters confined to compute domain (data domain provides velocity halo)
!  explicit wind and Stokes drift added for 2D (surface) drifters
!  drifter MPI communication reorganized by Klaus Ketelsen (!kk)
!  get a large number (10000) of drifters run in 2D and 3D by adapted MPI buffer
!</DESCRIPTION>
!
!<NAMELIST NAME="ocean_drifters_nml">
!  <DATA NAME="use_this_module" TYPE="logical">
!  Must be true to run this module. Default use_this_module=.false.
!  </DATA> 
!  <DATA NAME="output_interval" TYPE="integer">
!  Interval in timesteps between drifter writes.
!  Drifter time steps are two ocean time steps because 2 cycles are needed for Runge-Kutta scheme of 4th order.
!  </DATA> 
!  <DATA NAME="h_drift" TYPE="real">
!  Factor to enhance drift by current, default h_drift=1e0.
!  Drifter time steps are two ocean time steps because 2 cycles are needed for Runge-Kutta scheme of 4th order.
!  </DATA> 
!  <DATA NAME="debug_this_module" TYPE="logical">
!  Must be true for debugging, e.g. list of initial drifters. Default debug_this_module=.false.
!  </DATA> 
!  <DATA NAME="debug_drift" TYPE="logical">
!  Must be true for debugging with constant drift velocities. Default debug_drift=.false.
!  </DATA> 
!  <DATA NAME="u_debug" TYPE="real">
!  Constant zonal velocity for debugging, default u_debug=0.1 [m/s]
!  </DATA> 
!  <DATA NAME="v_debug" TYPE="real">
!  Constant meridional velocity for debugging, default v_debug=0.1 [m/s]
!  </DATA> 
!  <DATA NAME="w_debug" TYPE="real">
!  Constant vertical velocity for debugging, default w_debug=1e-5 [m/s]
!  </DATA> 
!  <DATA NAME="use_noise" TYPE="logical">
!  Must be true to add random noise to horizontal drifter positions after each time step. Default use_noise=.false.
!  Noise is characterized by velocity scales which adept to the drifter time step. It is applied as random
!  displacements to new positions of drifters in drifters_comm_update.
!  </DATA> 
!  <DATA NAME="a_noise" TYPE="real">
!  Amplitude scale of random noise yielding displacements between [-1,+1]*a_noise*timestep. Default 0.01 [m/s].
!  </DATA> 
!  <DATA NAME="v_noise" TYPE="real">
!  Cutoff velocity for random noise. Default 2*a_noise [m/s].
!  Suppressing random displacements if drifters are moving slower than v_noise. This avoids crossing into land cells, 
!  but there is also no random dispersion for resting drifters in open sea.
!  </DATA> 
!  <DATA NAME="use_wind_drift" TYPE="logical">
!  Must be true to add wind drift to horizontal surface velocities for 2D drifters. 
!  Requires the wind stress.
!  Default use_wind_drift=.false.
!  </DATA> 
!  <DATA NAME="c_wind" TYPE="real">
!  Factor to scale wind stress [N/m²] to wind drift velocity [m/s]. Default 1.0
!  </DATA> 
!  <DATA NAME="use_stokes_drift" TYPE="logical">
!  Must be true to add Stokes drift to horizontal surface velocities for 2D drifters at each time step. 
!  Requires a wave model providing the Stokes drift components.
!  Default use_stokes_drift=.false.
!  </DATA> 
!  <DATA NAME="c_stokes" TYPE="real">
!  Factor to tune influence of Stokes drift on surface drifters. Default 1.0
!  </DATA> 
!</NAMELIST>

  use drifters_mod
  use drifters_comm_mod,only       : id_trans
  use fms_mod
  use mpp_domains_mod
  use mpp_mod,          only : input_nml_file, mpp_pe, mpp_error, FATAL, stdout, stdlog ! iow !
  use mpp_mod,          only : mpp_clock_id, mpp_clock_begin, mpp_clock_end
  use time_manager_mod, only : time_type, get_time

  use ocean_domains_mod,    only : get_local_indices
  use ocean_parameters_mod, only : rho0r
  use ocean_types_mod,      only : ocean_grid_type, ocean_domain_type
  use ocean_types_mod,      only : ocean_time_type, ocean_velocity_type
  use ocean_types_mod,      only : ocean_adv_vel_type, ocean_prog_tracer_type
  use wave_types_mod,       only : ocean_wave_type ! iow ! for Stokes drift

  implicit none

  type(drifters_type),save :: drfts

  character(len=128) :: ermesg
  integer :: isc,iec,jsc,jec,isd,ied,jsd,jed,nk
  integer :: dt ! ocean time step in seconds
  integer :: id_init, id_upds
  real :: xcmin, xcmax,ycmin,ycmax,xdmin,xdmax,ydmin,ydmax,xmin,xmax
  real, dimension(:), allocatable :: x1d, y1d, z1d, dlon, dlat
  real, dimension(:,:,:), allocatable, save :: u, v, w, lon3d, lat3d, depth3d

  ! nml variables 
  logical :: use_this_module=.false.
  integer :: output_interval=1
  real    :: h_drift=1e0  ! iow ! factor to enhance drift by currents
! iow ! set debugging option
  logical :: debug_this_module=.false.
! iow ! constant velocities for debugging
  logical :: debug_drift=.false.
  real    :: u_debug=0.1  ! iow ! zonal velocity for debugging [m/s]
  real    :: v_debug=0.1  ! iow ! meridional velocity for debugging [m/s]
  real    :: w_debug=1e-5 ! iow ! vertical velocity for debugging [m/s]
! iow ! default velocity scale factors for noise
  logical :: use_noise=.false.
  real    :: a_noise=1e-2 ! iow ! amplitude scale [m/s]
  real    :: v_noise=2e-2 ! iow ! cutoff scale [m/s]
! iow ! explicit 2D wind drift (derived from wind stress)
  logical :: use_wind_drift=.false.
  real    :: c_wind=1.0   ! iow ! scale factor wind stress [N/m²] to wind drift [m/s]
! iow ! Stokes drift (imported from wave model)
  logical :: use_stokes_drift=.false.
  real    :: c_stokes=1.0 ! iow ! tuning factor for Stokes drift [m/s]
  namelist /ocean_drifters_nml/ use_this_module, output_interval, h_drift &
          , debug_this_module, debug_drift, u_debug, v_debug, w_debug &
          , use_noise, a_noise, v_noise &
          , use_wind_drift, c_wind, use_stokes_drift, c_stokes

! iow ! copy ocean_wave_nml from ocean_wave.F90 for local check of calc_stokes_drift
  logical  :: use_wave_module   = .false.
  logical  :: write_a_restart   = .true.
  logical  :: use_TMA           = .true.
  logical  :: filter_wave_mom   = .true.
  logical  :: damp_where_ice    = .true.
  real     :: wavedamp = -10.
  logical  :: calc_stokes_drift = .false.         
  real     :: stokes_fac   = 1.0                  
  real     :: spread_loss  = 1.0                  
  real     :: crit_ice     = 50.0                 
  real     :: stokes_p_lim = 1.0                  
  real     :: stokes_a_lim = 0.005365             
  namelist /ocean_wave_nml/ use_this_module, debug_this_module, write_a_restart &
           , use_TMA, filter_wave_mom, damp_where_ice, wavedamp &
           , calc_stokes_drift, stokes_fac, spread_loss, crit_ice, stokes_p_lim, stokes_a_lim

! iow ! copy ocean_sbc_nml from ocean_sbc.F90 for local check of read_stokes_drift
  logical :: use_waterflux                  =.true.
  logical :: waterflux_tavg                 =.false.
  logical :: use_waterflux_override_calving =.false.
  logical :: use_waterflux_override_fprec   =.false.
  logical :: use_waterflux_override_evap    =.false.
  logical :: rotate_winds                   =.false.
  logical :: taux_sinx                      =.false.
  logical :: tauy_siny                      =.false.
  logical :: runoffspread                   =.false.
  logical :: calvingspread                  =.false.
  logical :: salt_restore_under_ice         =.true.
  logical :: salt_restore_as_salt_flux      =.true.
  logical :: zero_net_salt_restore          =.false.
  logical :: zero_net_salt_correction       =.false.
  logical :: zero_net_water_restore         =.false.
  logical :: zero_net_water_correction      =.false.
  logical :: zero_net_water_coupler         =.false.
  logical :: zero_net_water_couple_restore  =.false.
  logical :: zero_net_pme_eta_restore       =.false.
  logical :: debug_water_fluxes             =.false.
  logical :: zero_water_fluxes              =.false. 
  logical :: zero_pme_fluxes                =.false. 
  logical :: zero_calving_fluxes            =.false. 
  logical :: zero_runoff_fluxes             =.false. 
  logical :: zero_river_fluxes              =.false. 
  logical :: convert_river_to_pme           =.false.
  logical :: zero_heat_fluxes               =.false. 
  logical :: zero_surface_stress            =.false.
  logical :: read_restore_mask              =.false. 
  logical :: restore_mask_gfdl              =.false.
  logical :: land_model_heat_fluxes         =.false. 
  logical :: do_flux_correction             =.false.
  logical :: sbc_heat_fluxes_const          =.false.
  logical :: sbc_heat_fluxes_const_seasonal =.false.
  logical :: use_constant_sss_for_restore   =.false.
  logical :: use_constant_sst_for_restore   =.false.
  logical :: use_ideal_runoff               =.false.
  logical :: use_ideal_calving              =.false.
  logical :: read_stokes_drift              =.false.
  logical :: do_langmuir                    =.false.
  real    :: constant_sss_for_restore       = 35.0
  real    :: constant_sst_for_restore       = 12.0
  real    :: sbc_heat_fluxes_const_value    = 0.0    ! W/m2
  real    :: ice_salt_concentration         = 0.005  ! kg/kg
  real    :: runoff_salinity                = 0.0    ! psu
  real    :: runoff_temp_min                = 0.0    ! degC
  real    :: temp_restore_tscale            = -30.
  real    :: salt_restore_tscale            = -30.
  real    :: eta_restore_tscale             = -30.
  real    :: max_ice_thickness              = 5.0 
  real    :: salinity_ref                   = 35.0
  real    :: max_delta_salinity_restore     = -0.5
  real    :: temp_correction_scale          = 0.0
  real    :: salt_correction_scale          = 0.0
  real    :: tau_x_correction_scale         = 0.0
  real    :: tau_y_correction_scale         = 0.0
  logical :: constant_hlf                   = .true.
  logical :: constant_hlv                   = .true.
  logical :: avg_sfc_velocity               = .true.
  logical :: avg_sfc_temp_salt_eta          = .true.
  logical :: use_full_patm_for_sea_level    = .false. 
  logical :: do_bitwise_exact_sum           = .true.
  namelist /ocean_sbc_nml/ temp_restore_tscale, salt_restore_tscale, salt_restore_under_ice, salt_restore_as_salt_flux,        &   
           eta_restore_tscale, zero_net_pme_eta_restore,                                                                       &  
           rotate_winds, taux_sinx, tauy_siny, use_waterflux, waterflux_tavg, max_ice_thickness, runoffspread, calvingspread,  &
           use_waterflux_override_calving, use_waterflux_override_evap, use_waterflux_override_fprec,                          &
           salinity_ref, zero_net_salt_restore, zero_net_water_restore, zero_net_water_coupler, zero_net_water_couple_restore, &
           zero_net_salt_correction, zero_net_water_correction,                                                                &
           debug_water_fluxes, zero_water_fluxes, zero_calving_fluxes, zero_pme_fluxes, zero_runoff_fluxes, zero_river_fluxes, &
           convert_river_to_pme, zero_heat_fluxes, zero_surface_stress, avg_sfc_velocity, avg_sfc_temp_salt_eta,               &
           ice_salt_concentration, runoff_salinity, runoff_temp_min, read_restore_mask, restore_mask_gfdl,                     &
           land_model_heat_fluxes, use_full_patm_for_sea_level, max_delta_salinity_restore, do_flux_correction,                &
           temp_correction_scale, salt_correction_scale, tau_x_correction_scale, tau_y_correction_scale, do_bitwise_exact_sum, &
           sbc_heat_fluxes_const, sbc_heat_fluxes_const_value, sbc_heat_fluxes_const_seasonal,                                 &
           use_constant_sss_for_restore, constant_sss_for_restore, use_constant_sst_for_restore, constant_sst_for_restore,     &
           use_ideal_calving, use_ideal_runoff, constant_hlf, constant_hlv, read_stokes_drift, do_langmuir                      

contains

  subroutine ocean_drifters_init(Domain,Grid, Time, T_prog, Velocity, Adv_vel)

    type(ocean_domain_type),   intent(in)  :: Domain
    type(ocean_grid_type),     intent(in)  :: Grid
    type(ocean_time_type),     intent(in)  :: Time
    type(ocean_prog_tracer_type), intent(in), dimension(:) :: T_prog
    type(ocean_velocity_type), intent(in)  :: Velocity
    type(ocean_adv_vel_type),  intent(in)  :: Adv_vel
    
    integer :: i,j,k, days, npes, ioun, ierr, io_status
    integer, dimension(:), allocatable :: pelist

# ifdef _DEBUG
    write (stdout(),*) 'ocean_drifters_init starts on PE',mpp_pe()
# endif
    id_init   = mpp_clock_id('(Drifters initialization) ',grain=CLOCK_ROUTINE)
    id_upds   = mpp_clock_id('(Drifters updates) ',grain=CLOCK_ROUTINE)
!kk Additional timings
    id_trans  = mpp_clock_id('(Drifters comm_update) ',grain=CLOCK_ROUTINE)
    id_fields = mpp_clock_id('(Drifters comm_set_field) ',grain=CLOCK_ROUTINE)
    id_compk  = mpp_clock_id('(Drifters compute_k) ',grain=CLOCK_ROUTINE)
  
    call mpp_clock_begin(id_init)

! iow ! read calc|read_stokes_drift from "imported" namelists
! iow ! must be done first not to overwrite settings from ocean_drifters_nml
#ifdef INTERNAL_FILE_NML
    read (input_nml_file, nml=ocean_wave_nml, iostat=io_status)
    ierr = check_nml_error(io_status,'ocean_wave_nml')
#else
    ioun = open_namelist_file()
    read (ioun, ocean_wave_nml,iostat=io_status)
    call close_file (ioun)
#endif
    use_wave_module=use_this_module ! iow ! keep for wave module
    use_this_module=.false.         ! iow ! reset for drifter module
#ifdef INTERNAL_FILE_NML
    read (input_nml_file, nml=ocean_sbc_nml, iostat=io_status)
    ierr = check_nml_error(io_status,'ocean_sbc_nml')
#else
    ioun = open_namelist_file()
    read (ioun, ocean_sbc_nml,iostat=io_status)
    call close_file (ioun)
#endif

    ! provide for namelist over-ride of defaults 
#ifdef INTERNAL_FILE_NML
    read (input_nml_file, nml=ocean_drifters_nml, iostat=io_status)
    ierr = check_nml_error(io_status,'ocean_drifters_nml')
#else
    ioun = open_namelist_file()
    read (ioun, ocean_drifters_nml,iostat=io_status)
    ierr = check_nml_error(io_status,'ocean_drifters_nml')
    call close_file (ioun)
#endif
# ifdef _DEBUG
    write (stdout(),*) 'ocean_drifters_init ocean_drifters_nml read on PE',mpp_pe()
! iow ! the following shows that PEs read input output_interval=1.0 arbitrarily as 1 or 0 !!!
    write (stdout(),*) 'ocean_drifters_init output_interval on PE',mpp_pe(),output_interval
# endif
    write (stdout(),'(/)')
    write (stdout(), ocean_drifters_nml)

! iow ! check output_interval to avoid arbitrary hangup by zero input
    if (output_interval.lt.1) then
      write (stdout(),*) '==> WARNING: output_interval =',output_interval,' reset to 1'
      output_interval=1
    endif
! iow ! check drift enhancement 
    if (h_drift.ne.1e0) then
      write (stdout(),*) '--> drift by currents is enhanced by h_drift =',h_drift
      if (h_drift.lt.1e0) then
        write (stdout(),*) '==> WARNING: drift by currents is diminished'
      endif
      if (h_drift.le.0e0) then
        write (stdout(),*) '==> WARNING: drift factor reset to h_drift = 1e0'
        h_drift=1e0
      endif
    endif

    call get_local_indices(Domain,isd,ied,jsd,jed,isc,iec,jsc,jec)
    nk=Grid%nk

# ifdef _DEBUG
    write (stdout(),*) 'ocean_drifters_init get_local_indices on PE',mpp_pe(),isd,ied,jsd,jed,isc,iec,jsc,jec,nk
# endif

    if(.not. use_this_module) return 

    ! iow ! save debugging option to drifters%core
    drfts%core%debug = debug_this_module
    ! iow ! debugging with constant velocities
    if (debug_drift .and. mpp_pe().eq.0) then
      write (stdout(),*)
      write (stdout(),*) '--> ocean_drifters_init on',mpp_pe()
      write (stdout(),*) '--> debugging with constant velocities'
      write (stdout(),*) '--> u_debug =',u_debug,' [m/s]'
      write (stdout(),*) '--> v_debug =',v_debug,' [m/s]'
      write (stdout(),*) '--> w_debug =',w_debug,' [m/s]'
    endif

    ! Determine local compute and data domain boundaries
    ! Note: This is only appropriate for spherical grids
    ! will not work for the tripolar grid region, for instance.

    allocate(x1d(isd:ied),y1d(jsd:jed),z1d(nk))
    allocate(u(isd:ied,jsd:jed,nk),v(isd:ied,jsd:jed,nk),w(isd:ied,jsd:jed,nk))
    allocate(lon3d(isd:ied,jsd:jed,nk),lat3d(isd:ied,jsd:jed,nk),depth3d(isd:ied,jsd:jed,nk))    
    allocate(dlon(isd:ied),dlat(jsd:jed))

    u=0.0;v=0.0;w=0.0
!!     
!!     do i=isc,iec
!!        x1d(i)=Grid%grid_x_t(i)
!!     enddo
!!     
!!     if (isc.eq.1) then
!!         x1d(isd)=Grid%grid_x_t(Grid%ni)-360.0
!!     else
!!         x1d(isd)=Grid%grid_x_t(isd)
!!     endif
!! 
!!     if (iec.eq.Grid%ni) then
!!         x1d(ied)=Grid%grid_x_t(1)+360.0
!!     else
!!         x1d(ied)=Grid%grid_x_t(ied)
!!     endif
!!     
!!     do j=jsc,jec
!!        y1d(j)=Grid%grid_y_t(j)
!!     enddo
!! 
!!     if (jsc.eq.1) then
!!         y1d(jsd) = Grid%grid_y_t(1)-1.e-5
!!     else
!!         y1d(jsd) = Grid%grid_y_t(jsd)
!!     endif
!! 
!!     if (jec.eq.Grid%nj) then
!!         y1d(jed) = Grid%grid_y_t(Grid%nj)+1.e-5
!!     else
!!         y1d(jed) = Grid%grid_y_t(jed)
!!     endif
    
! iow ! set x1d, y1d, z1d used as velocity axes below also for defining drifter domain boundaries
! iow ! this applies for a B-grid to identify t-grid cell (i,j) to interpolate (u,v) on its edges
! iow ! moreover halo=1 is assumed here, i.e isd=isc-1, ied=iec+1, jsd=jsc+1, jed=jec+1
! iow ! otherwise isd+1,...,isc-1 and jsd+1,...,jsc-1 remain undefined
    do i=isc,iec
       x1d(i)=Grid%grid_x_u(i)
    enddo
    
    if (isc.eq.1) then
        x1d(isd)=Grid%grid_x_u(Grid%ni)-360.0
    else
        x1d(isd)=Grid%grid_x_u(isd)
    endif

    if (iec.eq.Grid%ni) then
        x1d(ied)=Grid%grid_x_u(1)+360.0
    else
        x1d(ied)=Grid%grid_x_u(ied)
    endif
    
    do j=jsc,jec
       y1d(j)=Grid%grid_y_u(j)
    enddo

    if (jsc.eq.1) then
        y1d(jsd) = Grid%grid_y_u(1)-1.e-5
    else
        y1d(jsd) = Grid%grid_y_u(jsd)
    endif

    if (jec.eq.Grid%nj) then
        y1d(jed) = Grid%grid_y_u(Grid%nj)+1.e-5
    else
        y1d(jed) = Grid%grid_y_u(jed)
    endif

! iow ! not clear what to do between surface & zt(1), resp. zt(nk) & sea bottom
! iow ! use z1d(0:nk+1) ?
! iow ! how to tackle with changing zstar depths & for partial bottom cells ?
    do k=1,nk
       z1d(k)=Grid%zt(k)
    enddo
    
! iow ! compute domain bounded by u-grid coordinates
   if (isc.eq.1) then
     xcmin = 2*x1d(isc)-x1d(isc+1)
   else
     xcmin = x1d(isc-1)
   endif
   xcmax = x1d(iec)
   if (isc.eq.1) then
     ycmin = 2*y1d(isc)-y1d(isc+1)
   else
     ycmin = y1d(jsc-1)
   endif
   ycmax = y1d(jec)    

    xdmin = x1d(isd)
    xdmax = x1d(ied)
    ydmin = y1d(jsd)
    ydmax = y1d(jed)    

# ifdef _DEBUG
    write (stdout(),*) 'ocean_drifters_init on PE xdmin,xdmax,ydmin,ydmax',mpp_pe(),xdmin,xdmax,ydmin,ydmax
    write (stdout(),*) 'ocean_drifters_init on PE xcmin,xcmax,ycmin,ycmax',mpp_pe(),xcmin,xcmax,ycmin,ycmax
# endif

!!    xmin = Grid%grid_x_t(Grid%ni)-360.; xmax = Grid%grid_x_t(1)+360.
    xmin = Grid%grid_x_u(Grid%ni)-360.; xmax = Grid%grid_x_u(1)+360. ! iow !
    
    call get_time(Time%Time_Step,dt,days)
    
    call drifters_new(drfts, &
       & input_file ='INPUT/drifters_inp.nc'  ,  &
       & output_file='DRIFTERS/drifters_out.nc', &
       & ermesg=ermesg)
    if (ermesg/='') then
      call mpp_error(FATAL, ermesg)
    endif
! (iow) check if drifter fields correspond to T_prog%name(1|2)
    if (size(drfts%input%field_names).ne.2) then
      ermesg='this version of drifter code works only for first 2 tracer fields'
      call mpp_error(FATAL, ermesg)
    else
      write (stdout(),*) '==> check that drifter fields correspond to first 2 tracer fields'
      write (stdout(),*) '--> number tracer_field T_prog%name'      
      do i = 1, size(drfts%input%field_names)
        write(stdout(),'(i4,32a,1x,32a)') i,trim(drfts%input%field_names(i)),trim(T_prog(i)%name)
      enddo
    endif

! (iow) check if wind drift is set for 2D drifters
    if (use_wind_drift) then
      if (drfts%core%nd.eq.2) then
        if (mpp_pe().eq.0) then
          write (stdout(),*) '--> Note: use_wind_drift for 2D drifters'
          write (stdout(),*) '--> Wind drift scaled from wind stress by c_wind =',c_wind
        endif
      else
        call mpp_error(FATAL, 'use_wind_drift only for 2D drifters')
      endif
    endif

! (iow) check if Stokes drift is set for 2D drifters
    if (use_stokes_drift) then
      if (drfts%core%nd.eq.2) then
        if (mpp_pe().eq.0) then
          write (stdout(),*) '--> Note: use_stokes_drift for 2D drifters'
          write (stdout(),*) '--> Stokes drift modified by c_stokes =',c_stokes
! iow ! check if wave model is active
          if (.not.use_wave_module) then
            call mpp_error(FATAL, 'use_stokes_drift=.true. in ocean_drifters_nml needs use_this_module=.true. in ocean_wave_nml')
          endif
! iow ! check if calc_stokes_drift is set
          if (calc_stokes_drift) then
            write (stdout(),*) '--> found calc_stokes_drift=.true. in ocean_wave_nml'
          else
            call mpp_error(FATAL, 'use_stokes_drift=.true. in ocean_drifters_nml needs calc_stokes_drift=.true. in ocean_wave_nml')
          endif
! iow ! check if read_stokes_drift is set
          if (read_stokes_drift) then
            write (stdout(),*) '--> found read_stokes_drift=.true. in ocean_sbc_nml'
          else
            call mpp_error(FATAL, 'use_stokes_drift=.true. in ocean_drifters_nml needs read_stokes_drift=.true. in ocean_sbc_nml')
          endif
        endif
      else
        call mpp_error(FATAL, 'use_stokes_drift only for 2D drifters')
      endif
    endif

! (iow) check drifter input for debugging, drfts%core%nf, drfts%core%np not yet defined
!!    if (debug_this_module .and. mpp_pe().eq.0) then
    if (debug_this_module) then
      write (stdout(),*)
      write (stdout(),*) '--> ocean_drifters_init on PE',mpp_pe()
      write (stdout(),*) '--> nd, nf, np, npdim =',drfts%core%nd,size(drfts%input%field_names),size(drfts%input%ids),drfts%core%npdim 
      write (stdout(),*) '--> title = '//trim(drfts%input%title)
      write (stdout(),*) '--> j, position_names, position_units'
      do i = 1, size(drfts%input%position_names)
        write (stdout(),*) i,trim(drfts%input%position_names(i))//' '//trim(drfts%input%position_units(i))
      enddo
      write (stdout(),*) '--> j, field_names, field_units'
      do i = 1, size(drfts%input%field_names)
        write (stdout(),*) i,trim(drfts%input%field_names(i))//' '//trim(drfts%input%field_units(i))
      enddo
      write (stdout(),*) '--> j, velocity_names'
      do i = 1, size(drfts%input%velocity_names)
        write (stdout(),*) i,trim(drfts%input%velocity_names(i))
      enddo
! (iow) show input drifters later after distribution to PEs
    endif
 
    npes=mpp_npes()
    allocate(pelist(npes))
    call mpp_get_pelist(Domain%domain2d,pelist)

    drfts%comm%pe_beg=pelist(1)
    drfts%comm%pe_end=pelist(npes)
    
  ! set the initial time and dt
    drfts%time = 0.0
  ! (iow) drifter movement corresponds to 2*dt (checked with u,v=const.)
    drfts%dt   = 2.0*float(dt)
    drfts%core%add_noise = use_noise
    if (use_noise) then
  ! rescale to displacements & set factor of 2 to scale (random-0.5) to [-1:+1]
      drfts%core%anoise = 2.0*a_noise*drfts%dt ! max. displacements [m]
      drfts%core%vnoise = v_noise*drfts%dt ! cutoff for smaller displacements [m]
      if (debug_this_module .and. mpp_pe().eq.0) then
        write (stdout(),*) '--> add isotropic noise to drifter positions with scales:'
        write (stdout(),*) '--> amplitude scale',a_noise,' [m/s], yielding max. displacements of -/+',a_noise*drfts%dt,' [m]'
        write (stdout(),*) '--> cutoff velocity',v_noise,' [m/s], suppressing displacements below',v_noise*drfts%dt,' [m]'
      endif
    endif

  ! set the PE domain boundaries. xmin_comp/ymin_comp, xmax_comp/ymax_comp
  ! refer to the "compute" domain, which should cover densily the domain:
  ! xcmax[pe] = xcmin[pe_east]
  ! ycmax[pe] = ycmin[pe_north]
  ! Xmin_data/ymin_data, xmax_data/ymax_data refer to the "data" domain, which
  ! should be larger than the compute domain and therefore overlap:
  ! xdmax[pe] > xdmin[pe_east]
  ! ydmax[pe] > ydmin[pe_north]
  ! Particles in the overlap regions are tracked by several PEs. 

! iow ! debugging
    if (debug_this_module) then
      write(stdout(),*) 'ocean_drifters set comp domain on PE',mpp_pe(),isc,iec,jsc,jec,xcmin,xcmax,ycmin,ycmax
    endif

    call drifters_set_domain(drfts, &
       & xmin_comp=xcmin, xmax_comp=xcmax, &
       & ymin_comp=ycmin, ymax_comp=ycmax, &
       & xmin_data=xdmin, xmax_data=xdmax, &
       & ymin_data=ydmin, ymax_data=ydmax, &
       & xmin_glob=xmin , xmax_glob=xmax,  & ! this will set periodicity in x 
       & ermesg=ermesg)

! iow ! debugging
    if (debug_this_module) then
      write(stdout(),*) 'drifters_set_domain comp domain on PE',drfts%comm%pe,isc,iec,jsc,jec,drfts%comm%xcmin,drfts%comm%xcmax,drfts%comm%ycmin,drfts%comm%ycmax
    endif
  
  ! set neighboring PEs [domain2d is of type(domain2d)]
    call drifters_set_pe_neighbors(drfts, domain=Domain%domain2d, ermesg=ermesg)
    drfts%comm%pe = mpp_pe() ! save pe (iow)

  ! set the velocities axes. Each velocity can have different axes.
    call drifters_set_v_axes(drfts, component='u', &
       & x=x1d, y=y1d, z=z1d, ermesg=ermesg)

    call drifters_set_v_axes(drfts, component='v', &
       & x=x1d, y=y1d, z=z1d, ermesg=ermesg)

    call drifters_set_v_axes(drfts, component='w', &
       & x=x1d, y=y1d, z=z1d, ermesg=ermesg)  

! iow ! debugging
    if (debug_this_module) then
      write (stdout(),*) 'ocean_drifters set xu axes:',size(drfts%xu),drfts%xu(1:3),' ...',maxval(drfts%xu)
      write (stdout(),*) 'ocean_drifters set yu axes:',size(drfts%yu),drfts%yu(1:3),' ...',maxval(drfts%yu)
      write (stdout(),*) 'ocean_drifters set zu axes:',size(drfts%zu),drfts%zu(1:3),' ...',maxval(drfts%zu)
      write (stdout(),*) 'ocean_drifters set xv axes:',size(drfts%xv),drfts%xv(1:3),' ...',maxval(drfts%xv)
      write (stdout(),*) 'ocean_drifters set yv axes:',size(drfts%yv),drfts%yv(1:3),' ...',maxval(drfts%yv)
      write (stdout(),*) 'ocean_drifters set zv axes:',size(drfts%zv),drfts%zv(1:3),' ...',maxval(drfts%zv)
      write (stdout(),*) 'ocean_drifters set xw axes:',size(drfts%xw),drfts%xw(1:3),' ...',maxval(drfts%xw)
      write (stdout(),*) 'ocean_drifters set yw axes:',size(drfts%yw),drfts%yw(1:3),' ...',maxval(drfts%yw)
      write (stdout(),*) 'ocean_drifters set zw axes:',size(drfts%zw),drfts%zw(1:3),' ...',maxval(drfts%zw)
   endif

  ! Distribute the drifters across PEs
    call drifters_distribute(drfts, ermesg)

  ! (iow) show distribution of input drifters to PEs, drfts%core%np is now defined
    if (debug_this_module) then
      if (drfts%core%np.ge.1) then
        write(stdout(),*) '--> PE',mpp_pe(),' has',drfts%core%np,' initial drifters'
        write (stdout(),*) '--> id, input positions'
        if (drfts%core%nd.eq.2) then
          write(stdout(),'(i8,1p2e18.10)') (drfts%input%ids(i),drfts%input%positions(:,i),i=1,drfts%core%np)
        else
          write(stdout(),'(i8,1p3e18.10)') (drfts%input%ids(i),drfts%input%positions(:,i),i=1,drfts%core%np) 
        endif
      else
        write(stdout(),*) '--> PE',mpp_pe(),' has no initial drifters'
      endif
    endif

! iow ! save initial drifters
      call drifters_save(drfts, ermesg=ermesg)      

  ! Push the drifters. u_comp, v_comp etc are provided by the host code

! iow ! this applies only for regular grid steps & halo=1
!!     do j=jsc,jec
!!        dlat(j) = 0.5*(y1d(j+1)-y1d(j-1))
!!     enddo
! iow ! x1d & y1f refer to velocity coordinates bounding the t-cells
    do j=jsc,jec
      dlat(j) = y1d(j)-y1d(j-1)
    enddo
    if (jsc.eq.1) dlat(1) = y1d(2)-y1d(1)
    dlat(jsd)=dlat(jsc)
    dlat(jed)=dlat(jec)

!!     do i=isc,iec
!!        dlon(i) = 0.5*(x1d(i+1)-x1d(i-1))
!!     enddo
    do i=isc,iec
       dlon(i) = x1d(i)-x1d(i-1)
    enddo
    if (isc.eq.1) dlon(1) = x1d(2)-x1d(1)
    dlon(isd)=dlon(isc)
    dlon(ied)=dlon(iec)
    
    do k=1,nk
       do j=jsd,jed
          do i=isd,ied
             u(i,j,k)=h_drift*Velocity%u(i,j,k,1,Time%taup1)*Grid%dxtr(i,j)*dlon(i)
             v(i,j,k)=h_drift*Velocity%u(i,j,k,2,Time%taup1)*Grid%dytr(i,j)*dlat(j)
             w(i,j,k)=-1.0*Adv_vel%wrho_bt(i,j,k)*rho0r
             lon3d(i,j,k) = Grid%xt(i,j)
             lat3d(i,j,k) = Grid%yt(i,j)
             depth3d(i,j,k) = Grid%zt(k)                      
          enddo
       enddo
    enddo
! (iow) u,v,w=const for debugging (show that drifters move with 2*dt)
    if (debug_drift) then
       do k=1,nk
          do j=jsd,jed
             do i=isd,ied
                u(i,j,k)=u_debug*Grid%umask(i,j,k)*Grid%dxtr(i,j)*dlon(i) ! (iow) ! deg/s
                v(i,j,k)=v_debug*Grid%umask(i,j,k)*Grid%dytr(i,j)*dlat(j) ! (iow) ! deg/s
                w(i,j,k)=w_debug ! (iow) ! m/s
             enddo
          enddo
       enddo
    endif
! (iow) add explicit wind drift (only for 2D surface drifters)
    if (use_wind_drift) then
       do j=jsd,jed
          do i=isd,ied
             u(i,j,1)=u(i,j,1)+c_wind*Velocity%smf(i,j,1)*Grid%dxtr(i,j)*dlon(i)
             v(i,j,1)=v(i,j,1)+c_wind*Velocity%smf(i,j,2)*Grid%dytr(i,j)*dlat(j)
          enddo
       enddo
    endif
! (iow) add Stokes drift (only for 2D surface drifters)
    if (use_stokes_drift) then
       do j=jsd,jed
          do i=isd,ied
             u(i,j,1)=u(i,j,1)+c_stokes*Velocity%stokes_drift(i,j,0,1)*Grid%dxtr(i,j)*dlon(i)
             v(i,j,1)=v(i,j,1)+c_stokes*Velocity%stokes_drift(i,j,0,2)*Grid%dytr(i,j)*dlat(j)
          enddo
       enddo
    endif
  
  ! push the drifters ! (iow) for 2D or 3D !
# ifdef _DEBUG
    write (stdout(),*) 'ocean_drifters_init drifters_push on PE',mpp_pe()
# endif
    if (drfts%core%nd .eq. 2) then
      call drifters_push(drfts, u=u(:,:,1), v=v(:,:,1), ermesg=ermesg)
    else
      call drifters_push(drfts, u=u, v=v, w=w, ermesg=ermesg)
    endif

  ! check if RK4 integration is complete
    if(drfts%rk4_completed .and. mod(Time%itt0,output_interval) .eq. 0 ) then
      
      ! interpolate fields ! (iow) for 2D or 3D !      
      if (drfts%core%nd .eq. 2) then

        call drifters_set_field(drfts, index_field=1, x=x1d, y=y1d, &
           &    data=T_prog(1)%field(:,:,1,Time%taup1), ermesg=ermesg)

        call drifters_set_field(drfts, index_field=2, x=x1d, y=y1d, &
           &    data=T_prog(2)%field(:,:,1,Time%taup1), ermesg=ermesg)      
      
      else

        call drifters_set_field(drfts, index_field=1, x=x1d, y=y1d, z=z1d, &
           &    data=T_prog(1)%field(:,:,:,Time%taup1), ermesg=ermesg)

        call drifters_set_field(drfts, index_field=2, x=x1d, y=y1d, z=z1d, &
           &    data=T_prog(2)%field(:,:,:,Time%taup1), ermesg=ermesg)      

      endif
      
      ! save data 
      call drifters_save(drfts, ermesg=ermesg)
      
    endif

# ifdef _DEBUG
    write (stdout(),*) 'ocean_drifters_init ready on PE',mpp_pe()
# endif
  
    call mpp_clock_end(id_init)

  end subroutine ocean_drifters_init


  subroutine update_ocean_drifters(Velocity, Adv_vel, T_prog, Grid, Time, Domain) ! iow !

    type(ocean_velocity_type),    intent(in) :: Velocity
    type(ocean_adv_vel_type),     intent(in) :: Adv_vel
    type(ocean_prog_tracer_type), intent(in), dimension(:) :: T_prog
    type(ocean_grid_type),        intent(in) :: Grid
    type(ocean_time_type),        intent(in) :: Time
! iow ! needed for halo update
    type(ocean_domain_type)      :: Domain

    integer :: i,j,k

    if(.not. use_this_module) return 

# ifdef _DEBUG
    write (stdout(),*) 'update_ocean_drifters:',drfts%core%np,' drifters on PE',drfts%comm%pe,' at itt',Time%itt0
    if (drfts%core%np.ge.1) then
      write (stdout(),*) 'update_ocean_drifters: drifter ids ',minval(drfts%core%ids(1:drfts%core%np)),'  ...',maxval(drfts%core%ids(1:drfts%core%np))
    endif
# endif
    call mpp_clock_begin(id_upds)

    do k=1,nk
       do j=jsd,jed
          do i=isd,ied
!!# ifdef _DEBUG
!!    if (Grid%umask(i,j,k).gt.0.and.k.eq.1) & ! iow ! 2D drifters ! 
!!    write (stdout(),*) 'velo_drift',Time%itt0,mpp_pe(),i,j,k &
!!             ,h_drift*Velocity%u(i,j,k,1,Time%taup1)*Grid%dxtr(i,j)*dlon(i) &
!!             ,h_drift*Velocity%u(i,j,k,2,Time%taup1)*Grid%dytr(i,j)*dlat(j) &
!!             ,-1.0*Adv_vel%wrho_bt(i,j,k)*rho0r
!!# endif
             u(i,j,k)=h_drift*Velocity%u(i,j,k,1,Time%taup1)*Grid%dxtr(i,j)*dlon(i)
             v(i,j,k)=h_drift*Velocity%u(i,j,k,2,Time%taup1)*Grid%dytr(i,j)*dlat(j)
             w(i,j,k)=-1.0*Adv_vel%wrho_bt(i,j,k)*rho0r
          enddo
       enddo
    enddo
! (iow) u,v,w=const for debugging (show that drifters move with 2*dt)
    if (debug_drift) then
       do k=1,nk
          do j=jsd,jed
             do i=isd,ied
                u(i,j,k)=u_debug*Grid%umask(i,j,k)*Grid%dxtr(i,j)*dlon(i) ! (iow) ! deg/s
                v(i,j,k)=v_debug*Grid%umask(i,j,k)*Grid%dytr(i,j)*dlat(j) ! (iow) ! deg/s
                w(i,j,k)=w_debug ! (iow) ! m/s
             enddo
          enddo
       enddo
    endif
! (iow) add explicit wind drift (only for 2D surface drifters)
    if (use_wind_drift) then
       do j=jsd,jed
          do i=isd,ied
!!# ifdef _DEBUG
!!    if (Grid%umask(i,j,1).gt.0) & 
!!    write (stdout(),*) 'wind_drift',Time%itt0,mpp_pe(),i,j,1 &
!!             ,c_wind*Velocity%smf(i,j,1)*Grid%dxtr(i,j)*dlon(i) &
!!             ,c_wind*Velocity%smf(i,j,2)*Grid%dytr(i,j)*dlat(j)
!!# endif
             u(i,j,1)=u(i,j,1)+c_wind*Velocity%smf(i,j,1)*Grid%dxtr(i,j)*dlon(i)
             v(i,j,1)=v(i,j,1)+c_wind*Velocity%smf(i,j,2)*Grid%dytr(i,j)*dlat(j)
          enddo
       enddo
    endif
! (iow) add Stokes drift (only for 2D surface drifters)
    if (use_stokes_drift) then
       do j=jsd,jed
          do i=isd,ied
             u(i,j,1)=u(i,j,1)+c_stokes*Velocity%stokes_drift(i,j,0,1)*Grid%dxtr(i,j)*dlon(i)
             v(i,j,1)=v(i,j,1)+c_stokes*Velocity%stokes_drift(i,j,0,2)*Grid%dytr(i,j)*dlat(j)
!!# ifdef _DEBUG
!!    if (Grid%umask(i,j,1).gt.0) then 
!!    write (stdout(),*) 'stokes_drift',Time%itt0,mpp_pe(),i,j,1 &
!!             ,c_stokes*Velocity%stokes_drift(i,j,0,1)*Grid%dxtr(i,j)*dlon(i) &
!!             ,c_stokes*Velocity%stokes_drift(i,j,0,2)*Grid%dytr(i,j)*dlat(j)
!!    write (stdout(),*) 'full_drift',Time%itt0,mpp_pe(),i,j,k,u(i,j,1),v(i,j,1)
!!    endif
!!# endif
          enddo
       enddo
    endif

! iow ! synchronize drifter velocities
! iow ! otherwise 4th Runge-Kutta time step will differ at domain boundaries
! iow ! leading to different drifter tracks in dependence on domain decomposition
    if (drfts%core%nd.eq.2) then
      call mpp_update_domains(u(:,:,1), v(:,:,1), Domain%domain2d,gridtype=BGRID_NE)
    else
      call mpp_update_domains(u(:,:,:), v(:,:,:), Domain%domain2d,gridtype=BGRID_NE)
    endif


    ! push the drifters ! (iow) for 2D or 3D !
# ifdef _DEBUG
    write (stdout(),*) 'update_ocean_drifters: push',drfts%core%np,' drifters on PE',drfts%comm%pe,' at itt',Time%itt0
    if (drfts%core%np.ge.1) then
! iow ! minval(drfts%core%positions(1,1:drfts%core%np)) to avoid blank drifters drfts%core%np+1 ... drfts%core%npdim
      write (stdout(),*) 'update_ocean_drifters: itt xlon',Time%itt0,minval(drfts%core%positions(1,1:drfts%core%np)),' ...',maxval(drfts%core%positions(1,:))
      write (stdout(),*) 'update_ocean_drifters: itt ylat',Time%itt0,minval(drfts%core%positions(2,1:drfts%core%np)),' ...',maxval(drfts%core%positions(2,:))
    endif
# endif
    if (drfts%core%nd .eq. 2) then
      call drifters_push(drfts, u=u(:,:,1), v=v(:,:,1), ermesg=ermesg)
    else
      call drifters_push(drfts, u=u, v=v, w=w, ermesg=ermesg)
    endif
 
    ! check if RK4 integration is complete
# ifdef _DEBUG
    write (stdout(),*) 'update_ocean_drifters: check if RK4 integration is complete on PE',drfts%comm%pe,' at itt',Time%itt0
# endif
    if(drfts%rk4_completed .and. mod(Time%itt0,output_interval) .eq. 0 ) then
# ifdef _DEBUG
      write (stdout(),*) 'update_ocean_drifters: interpolate fields on PE',drfts%comm%pe,' at itt',Time%itt0
# endif
      ! interpolate fields ! (iow) for 2D or 3D ! T_prog(1|2) = temp|salt !
      if (drfts%core%nd .eq. 2) then

        call drifters_set_field(drfts, index_field=1, x=x1d, y=y1d, &
           &    data=T_prog(1)%field(:,:,1,Time%taup1), ermesg=ermesg)

        call drifters_set_field(drfts, index_field=2, x=x1d, y=y1d, &
           &    data=T_prog(2)%field(:,:,1,Time%taup1), ermesg=ermesg)      
      
      else

        call drifters_set_field(drfts, index_field=1, x=x1d, y=y1d, z=z1d, &
           &    data=T_prog(1)%field(:,:,:,Time%taup1), ermesg=ermesg)

        call drifters_set_field(drfts, index_field=2, x=x1d, y=y1d, z=z1d, &
           &    data=T_prog(2)%field(:,:,:,Time%taup1), ermesg=ermesg)      

      endif
!!        ! interpolate fields
!!        call drifters_set_field(drfts, index_field=1, x=x1d, y=y1d, z=z1d, &
!!             &    data=lon3d, ermesg=ermesg)
!!
!!        call drifters_set_field(drfts, index_field=2, x=x1d, y=y1d, z=z1d, &
!!             &    data=lat3d, ermesg=ermesg)             
!!
!!        call drifters_set_field(drfts, index_field=3, x=x1d, y=y1d, z=z1d, &
!!             &    data=depth3d, ermesg=ermesg)
!!
!!        call drifters_set_field(drfts, index_field=4, x=x1d, y=y1d, z=z1d, &
!!	     &    data=T_prog(1)%field(:,:,:,Time%taup1), ermesg=ermesg)
!!
!!        call drifters_set_field(drfts, index_field=5, x=x1d, y=y1d, z=z1d, &
!!             &    data=T_prog(2)%field(:,:,:,Time%taup1), ermesg=ermesg)

# ifdef _DEBUG
        write (stdout(),*) 'update_ocean_drifters: save',drfts%core%np,' drifters on PE',drfts%comm%pe,' at itt',Time%itt0
# endif
        call drifters_save(drfts, ermesg=ermesg)

    endif 

    call mpp_clock_end(id_upds)
    
  end subroutine update_ocean_drifters


  subroutine ocean_drifters_end(Grid)
    type(ocean_grid_type) :: Grid

    if(.not. use_this_module) return 

    ! write restart file, optionally with lon/lat data coordinates (under development)
    call drifters_write_restart(drfts, filename='RESTART/drifters_inp.nc', &
!       & x1=Grid%grid_x_u, y1=Grid%grid_y_u, geolon1=Grid%xt, &
!       & x2=Grid%grid_x_u, y2=Grid%grid_y_u, geolat2=Grid%yt, &
         ermesg=ermesg)  

    ! destroy
    call drifters_del(drfts, ermesg=ermesg)

    deallocate(x1d,y1d,z1d,u,v,w,lon3d,lat3d,depth3d,dlon,dlat)

  end subroutine ocean_drifters_end


end module ocean_drifters_mod
