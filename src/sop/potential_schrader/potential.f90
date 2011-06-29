!---------------------------------------------------------------------------------
! Attraction-repulsion potential Psi with minimum for rr = meanD, quadratic

IF (rd .GT. attractive_radius) THEN
    Psi_part = Psi_part + param_a * meanD**2 * &
        ((1._MK / rd - 1._MK ) **2 - 0.00_MK*rd **2)
ELSE
    Psi_part = Psi_part + meanD**2 * (-5._MK + rd)
ENDIF

Psi_max = MAX(Psi_max, 1._MK / (rd+1E-10))
