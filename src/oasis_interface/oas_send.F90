#IFDEF OASIS_IOW_ESM
SUBROUTINE oas_send( kid, kstep, pdata1, kinfo)

!!---------------------------------------------------------------------
!!              ***  ROUTINE oas_cos_snd  ***
!!
!! ** Purpose : - At each coupling time-step,this routine sends fields
!!                to the coupler or remote application.
!!----------------------------------------------------------------------

USE ice_grid_mod
USE mpp_domains_mod, only: mpp_get_compute_domain, Domain2d
USE oas_vardef

IMPLICIT NONE

!! * Arguments
!!
INTEGER, INTENT(IN)    ::  &
  kid    ! variable index in the array

INTEGER, INTENT(OUT)   ::  &
  kinfo  ! OASIS info argument

INTEGER, INTENT(IN)    ::  &
  kstep  ! time-step in seconds

REAL, INTENT(IN)          ::  &
  pdata1(iec,jec)

! Call OASIS at each time step but field is sent to other model only at coupling time step
! (accumulation otherwise, if asked in the namcouple configuration file)
!
CALL oasis_put(ssnd(kid)%vari_id, kstep, pdata1(isc:iec, jsc:jec), kinfo)

END SUBROUTINE oas_send
#ENDIF
