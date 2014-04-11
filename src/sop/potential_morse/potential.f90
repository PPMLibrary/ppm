!------------------------------------------------------------------------------!
!Morse potential, with parameters
!such that it is h-stable
!
!confining potential
 !Psi_part = Psi_part + meanD**2 * &
     !(-log(rd) + 0.5_mk*rd **2)



IF (rd .GT. attractive_radius .or. no_fusion) THEN

 Psi_part = Psi_part + meanD**2 * coeff *  &
     (-rho**(-4._mk*rd) + 0.8_mk*rho**(1._mk-5._mk*rd) &
      -Psi_at_cutoff)
ELSE
    Psi_part = Psi_part + meanD**2 * (-10._MK / SQRT(rd+0.1_mk)) &
     + meanD**2 *(-rho**(-4._mk*attractive_radius)+&
     0.8_mk*rho**(1._mk-5._mk*attractive_radius) + &
     10._MK/SQRT(attractive_radius+0.1_mk))
ENDIF

Psi_max = MAX(Psi_max, 1._MK / (rd+1E-10))
