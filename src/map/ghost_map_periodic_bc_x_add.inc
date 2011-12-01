         !-------------------------------------------------------------------
         !  copy periodic ghosts in the x-direction
         !-------------------------------------------------------------------
         xminf = min_phys(1) - ghostsize
         xmini = min_phys(1)

         DO i=1,Nall
            !----------------------------------------------------------------
            !  first those at the west boundary 
            !----------------------------------------------------------------
            IF (xp(1,i).GE.xminf.AND.xp(1,i).LT.xmini) THEN
               xp(1,i)   = xp(1,i) + len_phys(1)
            ENDIF
         ENDDO
         IF (isymm.EQ.0) THEN
            !----------------------------------------------------------------
            !  then the east bc, but only if we are not using symmetry
            !----------------------------------------------------------------
            xmaxf = max_phys(1) + ghostsize
            xmaxi = max_phys(1) 
            
            DO i=1,Nall
               IF  (xp(1,i).GE.xmaxi.AND.xp(1,i).LT.xmaxf) THEN
                  xp(1,i)   = xp(1,i) - len_phys(1)
               ENDIF
            ENDDO
         ENDIF 
