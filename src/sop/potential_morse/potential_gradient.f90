!------------------------------------------------------------------------------!
!Morse potential, with parameters
!such that it is h-stable

!confining potential
 !Psi_part = Psi_part + meanD**2 * &
     !(-log(rd) + 0.5_MK*rd **2)
 !gradPsi = meanD * (1.0_MK/rd - rd)


IF (rd .GT. attractive_radius .OR. no_fusion) THEN

    Psi_part = Psi_part + meanD**2 * coeff * &
        (-rho**(-4._MK*rd) + 0.8_MK*rho**(1.0_MK-5._MK*rd) &
        -Psi_at_cutoff)
    gradPsi = 2._MK* meanD * coeff * &
        4._MK*LOG(rho) * (rho**(-4._MK*rd) - rho**(1.0_MK-5._MK*rd))

ELSE

    Psi_part = Psi_part + meanD**2 * (-10._MK / SQRT(rd+0.1_MK)) &
        + meanD**2 *(-rho**(-4._MK*attractive_radius)+&
        0.8_MK*rho**(1.0_MK-5._MK*attractive_radius) + &
        10._MK/SQRT(attractive_radius+0.1_MK))
    gradPsi =  2._MK*meanD * 5._MK / (rd+0.1_MK)**1.5_MK

    !adaptation_ok = .FALSE.

ENDIF

Psi_max = MAX(Psi_max, 1._MK / (rd+1E-10))
