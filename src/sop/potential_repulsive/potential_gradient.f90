!------------------------------------------------------------------------------!
! Purely repulsive potential 
! Psi = gamma * meanD * (meanD/rr)**3
! gamma = 0.00001/3 
! gradPsi = -3*gamma*(meanD/rr)**4
!press = 1.4_MK
!------------------------------------------------------------------------------!

!confining potential
 !Psi_part = Psi_part + meanD**2 * &
     !(-log(rd) + 0.5_mk*rd **2)
 !gradPsi = meanD * (1._mk/rd - rd)

!Morse potential, with parameters
!such that it is h-stable

rho = opts%param_morse

IF (rd .GT. attractive_radius .or. no_fusion) THEN

 Psi_part = Psi_part + meanD**2 * coeff * &
     !(exp(-2._mk*rho*rd)-2._mk*exp(-rho*(rd+1._mk)))
     !(-exp(-rd)+rho*exp(-rd/0.8_mk))
     (-(rho)**(-4._mk*rd)+0.8_mk*rho**(1._mk-5._mk*rd))
 gradPsi = meanD * coeff * 0.25_mk/log(rho) * (-(rho)**(-4._mk*rd)+rho**(1._mk-5._mk*rd))
 !gradPsi = meanD * (-exp(-rd)+rho/0.8_mk*exp(-rd/0.8_mk))
 !gradPsi = meanD * (-2._mk*rho*exp(-rho*(rd+1._mk))*(1._mk-exp(rho*(1._mk-rd))))

ELSE
    Psi_part = Psi_part + meanD**2 * (-10._MK / sqrt(rd+0.1_mk)) &
        + meanD**2 *(-(rho)**(-4._mk*attractive_radius)+&
        0.8_mk*rho**(1._mk-5._mk*attractive_radius) + &
        10._MK/sqrt(attractive_radius+0.1_mk))
    gradPsi = - meanD * 5._MK / rd**1.5_mk
    adaptation_ok = .false.
ENDIF

!IF (rd .GT. attractive_radius) THEN
    !!Psi_part = Psi_part + 0.5_MK * meanD**2 /rd**2
    !!gradPsi = meanD / rd**3
    !!Psi_part = Psi_part + 0.5_MK / rd**2 + 1._MK/6._MK /rd**6
    !!gradPsi = 1._MK / (meanD * rd**3) + 1._MK / (meanD * rd**7)
    !Psi_part = Psi_part + meanD**2* &
        !(0.5_MK / rd**2 + 1._MK/6._MK /rd**6   + 0.0_mk * rc**6)
    !gradPsi = (meanD / rd**3) + (meanD / rd**7) - meanD*0.0_mk*rc**5
!ELSE
    !!Psi_part = Psi_part + meanD**2 * (-50._MK + 0.1_MK * rd)
    !Psi_part = Psi_part + meanD**2 * (-100._MK / rd)
    !gradPsi = - meanD * 100._MK / rd**2
    !!Psi_part = Psi_part -1._MK / rd
    !!gradPsi = - 1._MK / (meanD*rd**2)
!ENDIF

Psi_max = MAX(Psi_max, 1._MK / (rd+1E-10))
