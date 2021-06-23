#ifdef COUP_OAS
MODULE oas_vardef

!Controls, definitions and variables
!                   for OASIS communications
!

USE mpi
USE mod_oasis                    ! OASIS3-MCT module

IMPLICIT NONE

INCLUDE 'mpif.h'

  ! By default OASIS3-MCT exchanges data in double precision,
  ! but conversion to or from single precision data is supported in the interface
#ifdef NO_USE_DOUBLE_PRECISION
  INTEGER, PARAMETER :: wp = SELECTED_REAL_KIND(6,37)   ! real
#elif defined USE_DOUBLE_PRECISION
  INTEGER, PARAMETER :: wp = SELECTED_REAL_KIND(12,307) ! double
#endif

  !
  !CHARACTER(len=30), PARAMETER   :: data_gridname='grids.nc' ! file with the grids
  !CHARACTER(len=30), PARAMETER   :: data_maskname='masks.nc' ! file with the masks
  !CHARACTER(len=30)              :: data_filename, field_name

  ! Component name (6 characters) same as in the namcouple
  CHARACTER(len=6)   :: comp_name = 'momBS'
  CHARACTER(len=128) :: comp_out ! name of the output log file
  CHARACTER(len=3)   :: chout

  TYPE :: CPL_FLD
    CHARACTER(LEN = 16)     :: var_name    ! Name of the coupling field
    INTEGER                 :: vari_id     ! Id of the field
  END TYPE CPL_FLD

  TYPE(CPL_FLD), ALLOCATABLE :: &
    srcv(:),                 & ! All fields to be received
    ssnd(:)                    ! All fields to be sent

  INTEGER                 :: var_nodims(2)
  CHARACTER(LEN = 4)      :: grid_type   ! Grid type
  INTEGER                 :: var_type
  INTEGER               ::  ib
  INTEGER, PARAMETER    ::  il_nb_time_steps = 1 ! number of time steps
! get model time step!!!!  INTEGER, PARAMETER    ::  delta_t = 3600     ! time step ! 


  REAL (kind=wp), PARAMETER     :: field_ini = -1. ! initialisation of received fields
  INTEGER,        PARAMETER     :: err_msk = -10000

  INTEGER :: &
    comp_id,                & ! id returned by oasis_init_comp   
    ishape(2,2),             & ! shape of arrays passed to PSMILe
    local_comm,                 & ! Local communicator
    ierror,                  & ! return error code
    w_unit

  ! Grid parameters definition
  INTEGER                       :: part_id  ! use to connect the partition to the variables
                                            ! in oasis_def_var
  INTEGER                       :: var_sh(4) ! local dimensions of the arrays to the pe
                                             ! 2 x field rank (= 4 because fields are of rank = 2)
  INTEGER                       :: il_flag  ! Flag for grid writing by proc 0
  INTEGER                       :: itap_sec ! Time used in oasis_put/get
  INTEGER                       :: il_paral(5) ! Decomposition for each proc ; box partition

INTEGER                   :: &
  nfld,                    & !
  nfld_snd_tot,            & ! total number of fields to be sent
  nfld_rcv_tot,            & ! total number of fields to be received
  nfld_snd_oce,            & !
  nfld_rcv_oce

INTEGER(KIND=iintegers), ALLOCATABLE :: &
  nlev_snd_oce(:),         & !
  nlev_rcv_oce(:)

CHARACTER(LEN=16), ALLOCATABLE :: &
  nam_snd_oce(:),          & !
  nam_rcv_oce(:)             !

REAL (KIND=wp), ALLOCATABLE :: &
  exfld (:,:),             & ! Temporary buffer for receiving
  frcv  (:,:,:)              ! all fields recieved from coupled model

END MODULE oas_vardef
#endif
