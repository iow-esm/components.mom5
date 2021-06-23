#ifdef COUP_OAS
SUBROUTINE oas_init 

USE oas_vardef

IMPLICIT NONE
   
CHARACTER(LEN=6)   :: modname = 'mom5BS'    ! Name of the model

w_unit = 100 + rank
  WRITE(chout,'(I3)') w_unit
  comp_out=comp_name//'.out_'//chout
  !
  OPEN(w_unit,file=TRIM(comp_out),form='formatted')
  WRITE (w_unit,*) '-----------------------------------------------------------'
  WRITE (w_unit,*) TRIM(comp_name), ' Running with reals compiled as kind =',wp
  WRITE (w_unit,*) 'I am component ', TRIM(comp_name), ' rank :',rank
  WRITE (w_unit,*) '----------------------------------------------------------'
  CALL flush(w_unit)
      
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
  CALL prism_abort_proto( comp_id, 'oas_cos_init', 'Failure in oasis_get_localcomm' )
ENDIF

END SUBROUTINE oas_init
#endif
