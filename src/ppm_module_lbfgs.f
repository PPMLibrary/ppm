module ppm_module_lbfgs
  USE ppm_module_typedef

!   ----------------------------------------------------------          
!     DATA                                                              
!   ----------------------------------------------------------          
!                                                                       
PRIVATE

INTEGER LP, MP 
REAL(8) GTOL, STPMIN, STPMAX 
COMMON / LB3 / MP, LP, GTOL, STPMIN, STPMAX 
DATA MP, LP, GTOL, STPMIN, STPMAX / 6, 6, 9.0D-01, 1.0D-20,       &
    0.3D+0 /                                                         
    !10.0D+20 /                                                         
!                                                                       

PUBLIC LBFGS

contains

include 'lbfgs/ppm_lbfgs.f'

end module ppm_module_lbfgs

