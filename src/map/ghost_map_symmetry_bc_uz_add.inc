  !----------------------------------------------------------
  ! Symmetry ghosts in the z-direction.
  !----------------------------------------------------------

  !----------------------------------------------------------
  ! those at upper boundary, 
  ! even using symmetry inter-communication.
  ! No '=' in the inequality checking.
  !----------------------------------------------------------
  
  zmaxf = max_phys(3) + ghostsize
  zmaxi = max_phys(3) 
  
  DO i = 1, Nall
     
     IF ( xp(3,i) > zmaxi .AND. xp(3,i) < zmaxf ) THEN
        xp(3,i)   = 2.0_MK*zmaxf - xp(3,i)
     END IF
  END DO
