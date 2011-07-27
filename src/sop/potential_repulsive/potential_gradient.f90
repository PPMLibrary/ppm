!------------------------------------------------------------------------------!
! Purely repulsive potential 
! Psi = gamma * meanD * (meanD/rr)**3
! gamma = 0.00001/3 
! gradPsi = -3*gamma*(meanD/rr)**4
!press = 1.4_MK
!------------------------------------------------------------------------------!

IF (rd .GT. attractive_radius) THEN
    !Psi_part = Psi_part + 0.5_MK * meanD**2 /rd**2
    !gradPsi = meanD / rd**3
    !Psi_part = Psi_part + 0.5_MK / rd**2 + 1._MK/6._MK /rd**6
    !gradPsi = 1._MK / (meanD * rd**3) + 1._MK / (meanD * rd**7)
    Psi_part = Psi_part + meanD**2* &
        (0.5_MK / rd**2 + 1._MK/6._MK /rd**6   + 0.0_mk * rc**6)
    gradPsi = (meanD / rd**3) + (meanD / rd**7) - meanD*0.0_mk*rc**5
ELSE
    !Psi_part = Psi_part + meanD**2 * (-50._MK + 0.1_MK * rd)
    Psi_part = Psi_part + meanD**2 * (-100._MK / rd)
    gradPsi = - meanD * 100._MK / rd**2
    !Psi_part = Psi_part -1._MK / rd
    !gradPsi = - 1._MK / (meanD*rd**2)
ENDIF

Psi_max = MAX(Psi_max, 1._MK / (rd+1E-10))
