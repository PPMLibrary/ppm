MODULE ppm_module_lbfgs
  USE ppm_module_typedef

!   ----------------------------------------------------------          
!     DATA                                                              
!   ----------------------------------------------------------          
!                                                                       
IMPLICIT NONE 

PRIVATE

INTEGER, PARAMETER :: MP = 6, LP = 6
REAL(8),PARAMETER  :: GTOL   = 9.0D-01
REAL(8),PARAMETER  :: STPMIN = 1.0D-20
REAL(8),PARAMETER  :: STPMAX = 0.3D+0                                                          

!INTEGER :: LP, MP 
!REAL(8) :: GTOL, STPMIN, STPMAX 
!DATA MP, LP, GTOL, STPMIN, STPMAX / 6, 6, 9.0D-01, 1.0D-20,       &
    !0.3D+0 /                                                         
    !10.0D+20 /                                                         
!                                                                       

PUBLIC LBFGS

CONTAINS

#include "lbfgs/ppm_lbfgs.f"

END MODULE ppm_module_lbfgs

