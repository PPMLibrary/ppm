  !----------------------------------------------------------
  ! those at the south boundary.
  ! No '=' in the inequality checking.
  !----------------------------------------------------------
  
  yminf = min_phys(2) - ghostsize
  ymini = min_phys(2) 
  
  DO i = 1, Nall
     IF ( xp(2,i) > yminf .AND. xp(2,i) < ymini ) THEN
        xp(2,i)   = 2.0_MK*yminf - xp(2,i)
     END IF
  END DO
  
  
