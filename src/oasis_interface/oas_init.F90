#IFDEF OASIS_IOW_ESM

SUBROUTINE oas_init(local_comm) 

USE oas_vardef

IMPLICIT NONE

CHARACTER(LEN=6)     :: modname = 'mom5BS'    ! Name of the model
INTEGER,INTENT(OUT)  :: local_comm
!------------------------------------------------------------------
! Initialize OASIS 
!------------------------------------------------------------------

CALL oasis_init_comp( comp_id, modname, ierror )
IF( ierror /= OASIS_Success ) THEN
  CALL oasis_abort( comp_id, 'oas_init', 'Failure in oasis_init_comp' )
ENDIF

!------------------------------------------------------------------
! Get MPI communicator for local communication
!------------------------------------------------------------------

CALL oasis_get_localcomm( local_comm, ierror )
IF( ierror /= OASIS_Success ) THEN
  CALL oasis_abort( comp_id, 'oas_cos_init', 'Failure in oasis_get_localcomm' )
ENDIF

END SUBROUTINE oas_init
#ENDIF
