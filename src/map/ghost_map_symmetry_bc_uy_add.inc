  !----------------------------------------------------------
  ! symmetry ghosts in the y-direction.
  !----------------------------------------------------------
  
  !----------------------------------------------------------
  ! those at north boundary, 
  ! even using symmetry inter-comunication.
  ! No '=' in the inequality checking.
  !----------------------------------------------------------
  
  ymaxf = max_phys(2) + ghostsize
  ymaxi = max_phys(2)
  
  DO i = 1, Nall
     IF ( xp(2,i) > ymaxi .AND. xp(2,i) < ymaxf ) THEN
        xp(2,i)   = 2.0_MK*ymaxf - xp(2,i)
     END IF
  END DO
