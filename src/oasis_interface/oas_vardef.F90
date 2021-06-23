#IFDEF COUP_OAS
MODULE oas_vardef

!Controls, definitions and variables
!                   for OASIS communications
!

USE mpi
USE mod_oasis                    ! OASIS3-MCT module
!USE mod_oasis_grid

IMPLICIT NONE

!INCLUDE 'mpif.h'

  ! By default OASIS3-MCT exchanges data in double precision,
  ! but conversion to or from single precision data is supported in the interface
!#IFDEF NO_USE_DOUBLE_PRECISION
!  INTEGER, PARAMETER :: wp = SELECTED_REAL_KIND(6,37)   ! real
!#ELIF defined USE_DOUBLE_PRECISION
!  INTEGER, PARAMETER :: wp = SELECTED_REAL_KIND(12,307) ! double
!#ENDIF

  !
  !CHARACTER(len=30), PARAMETER   :: data_gridname='grids.nc' ! file with the grids
  !CHARACTER(len=30), PARAMETER   :: data_maskname='masks.nc' ! file with the masks
  !CHARACTER(len=30)              :: data_filename, field_name

  ! Component name (6 characters) same as in the namcouple
  CHARACTER(len=1)   :: MOM5_instance_letter = 'M' ! starting letter for OASIS3 variable names
  CHARACTER(len=6)   :: comp_name = 'momBS'
  CHARACTER(len=128) :: comp_out ! name of the output log file
  CHARACTER(len=3)   :: chout

  TYPE :: CPL_FLD
    CHARACTER(LEN = 1)      :: clgrid    ! Grid type
    LOGICAL                 :: laction = .TRUE.   ! To be coupled or not
    CHARACTER(LEN = 16)     :: var_name    ! Name of the coupling field
    INTEGER                 :: vari_id     ! Id of the field
  END TYPE CPL_FLD

#IFDEF OASIS_IOW_ESM
  INTEGER, PARAMETER :: nfld_snd_tot=18
  INTEGER, PARAMETER :: nfld_rcv_tot=44
#ELSE
  INTEGER, PARAMETER :: nfld_snd_tot=2
  INTEGER, PARAMETER :: nfld_rcv_tot=14
#ENDIF  

  TYPE(CPL_FLD), ALLOCATABLE :: &
    srcv(:),                 & ! All fields to be received
    ssnd(:)                    ! All fields to be sent

  INTEGER                 :: var_nodims(2)
  CHARACTER(LEN = 4)      :: grid_type   ! Grid type
  INTEGER                 :: var_type
  INTEGER               ::  ib
  INTEGER, PARAMETER    ::  il_nb_time_steps = 1 ! number of time steps
! get model time step!!!!  INTEGER, PARAMETER    ::  delta_t = 3600     ! time step ! 


  REAL (kind=8), PARAMETER     :: field_ini = -1. ! initialisation of received fields
  INTEGER,        PARAMETER     :: err_msk = -10000

  INTEGER :: &
    comp_id,                & ! id returned by oasis_init_comp   
    ishape(2,2),             & ! shape of arrays passed to PSMILe
    !local_comm,                 & ! Local communicator
    ierror,                  & ! return error code
    istatus,              &
    w_unit

  ! Grid parameters definition
  INTEGER                       :: part_id  ! use to connect the partition to the variables
                                            ! in oasis_def_var
  INTEGER                       :: var_sh(4) ! local dimensions of the arrays to the pe
                                             ! 2 x field rank (= 4 because fields are of rank = 2)
  INTEGER                       :: il_flag  ! Flag for grid writing by proc 0
  INTEGER                       :: itap_sec ! Time used in oasis_put/get
  INTEGER                       :: il_paral(5) ! Decomposition for each proc ; box partition

REAL (KIND=8), ALLOCATABLE :: &
  exfld1 (:,:),             & ! Temporary buffer for receiving
  frcv  (:,:,:)              ! all fields recieved from coupled model

INTEGER :: &
  OASIS_Rcv  = 1,          & ! return code if received field
  OASIS_idle = 0,          & ! return code if nothing was done by OASIS
  OASIS_Success = 0          ! return code if no error in OASIS


END MODULE oas_vardef
#ENDIF
