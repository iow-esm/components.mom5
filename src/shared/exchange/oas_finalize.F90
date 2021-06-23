#ifdef COUP_OAS
SUBROUTINE oas_finalize

! organization of MPI ending by OASIS 

USE oas_vardef

IMPLICIT NONE

! local variables
INTEGER :: ierror

DEALLOCATE( exfld, frcv, ssnd, srcv, STAT=ierror )

CALL oasis_terminate ( ierror )         

END SUBROUTINE oas_finalize
#endif
