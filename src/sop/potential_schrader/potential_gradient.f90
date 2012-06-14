!---------------------------------------------------------------------------------
! Attraction-repulsion potential Psi with minimum for rr = meanD, quadratic
! Psi = gamma * (2*meanD - 2*rr) * (meanD/rr)**2
! gradPsi = -2*gamma*(1._MK - meanD/rr)*(meanD/rr)**3

IF (rd .GT. attractive_radius) THEN
    Psi_part = Psi_part + meanD**2 * &
        ((1._MK / rd - 1._MK ) **2 - 0.00_MK*rd **2)
    gradPsi = 2._MK * meanD*&
        (-(1._MK - rd) / rd**3 + 0.0_MK * rd)
ELSE
    Psi_part = Psi_part + meanD**2 * (-5._MK + rd)
    gradPsi = -1._MK * meanD
ENDIF

Psi_max = MAX(Psi_max, 1._MK / (rd+1E-10))
