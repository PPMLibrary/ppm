  !----------------------------------------------------------
  ! symmetry ghosts in the x-direction.
  !----------------------------------------------------------

  !----------------------------------------------------------
  ! the east bc, even we are using symmetry.
  ! No '=' in the inequality checking.
  !----------------------------------------------------------
  
  xmaxf = max_phys(1) + ghostsize
  xmaxi = max_phys(1) 
  
  
  DO i = 1, Nall
     
     IF ( xp(1,i) > xmaxi .AND. xp(1,i) < xmaxf) THEN
        xp(1,i)   = 2.0_MK*xmaxf - xp(1,i)
     END IF
  END DO
