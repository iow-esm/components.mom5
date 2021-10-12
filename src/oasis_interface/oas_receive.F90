#IFDEF OASIS_IOW_ESM
SUBROUTINE oas_recieve( kid, kstep, pdata, kinfo )

!!---------------------------------------------------------------------
!!              ***  ROUTINE oas_cos_rcv  ***
!!
!! ** Purpose : - At each coupling time-step,this routine call fields
!!                from the coupler or remote application.
!!----------------------------------------------------------------------

USE ice_grid_mod
USE mpp_domains_mod, only: mpp_get_compute_domain, mpp_get_compute_domains, &
                           Domain2d, mpp_get_global_domain
USE oas_vardef

IMPLICIT NONE

INTEGER, INTENT(IN)   :: &
  kid    ! variable index in the array

INTEGER, INTENT(IN)   :: &
  kstep  ! time-step in seconds

REAL(kind=8), INTENT(OUT) :: & !, ALLOCATABLE :: &
  pdata(iec,jec)

INTEGER, INTENT(OUT)  :: &
  kinfo  ! OASIS info argument

!!
!! local variables
!!
LOGICAL :: &
  llaction
!!--------------------------------------------------------------------

!ALLOCATE(pdata(isc:iec, jsc:jec), stat = ierror )


!  Masked points are not modified by OASIS
!  = buffer set to 0 before calling OASIS
 exfld1(:,:) = 0.0
 pdata(:,:) = 0.0

  CALL oasis_get(srcv(kid)%vari_id, kstep, exfld1(isc:iec, jsc:jec), kinfo )         
 !  CALL oasis_get(srcv(kid)%vari_id, kstep, pdata(isc:iec, jsc:jec), kinfo ) 
  
  IF ( ierror /= OASIS_Success) THEN
    CALL oasis_abort(srcv(kid)%vari_id, 'oas_cos_rcv', &
      'Failure in oasis_get for '//TRIM(srcv(kid)%var_name) )
  ENDIF

  llaction = .FALSE.
  IF( kinfo == OASIS_Recvd   .OR. kinfo == OASIS_FromRest .OR.   &
      kinfo == OASIS_RecvOut .OR. kinfo == OASIS_FromRestOut .OR. kinfo == OASIS_Input)  llaction = .TRUE.
      
 
  ! If coupling time step
  IF ( llaction ) THEN

    ! Declare to calling routine that OASIS provided coupling field
    kinfo = OASIS_Rcv

    ! Update array which contains coupling field (only on valid shape)
     pdata(isc:iec, jsc:jec) = exfld1(isc:iec, jsc:jec)

  ELSE
    ! Declare to calling routine that OASIS did not provide coupling field
    kinfo = OASIS_idle     
  ENDIF

!DEALLOCATE( pdata, STAT=ierror )

END SUBROUTINE oas_recieve

#ENDIF
