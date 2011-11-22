  !----------------------------------------------------------
  ! symmetry ghosts in the z-direction.
  !----------------------------------------------------------
  
  !----------------------------------------------------------
  ! those at down boundary.
  ! No '=' in the inequality checking.
  !----------------------------------------------------------
  
  zminf = min_phys(3) - ghostsize
  zmini = min_phys(3) 
  
  DO i = 1, Nall
     
     IF ( xp(3,i) > zminf .AND. xp(3,i) < zmini ) THEN
        xp(3,i)   = 2.0_MK*zminf - xp(3,i)
     END IF
  END DO
