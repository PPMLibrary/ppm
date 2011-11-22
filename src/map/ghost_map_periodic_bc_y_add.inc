         !-------------------------------------------------------------------
         !  copy periodic ghosts in the y-direction
         !-------------------------------------------------------------------
         yminf = min_phys(2) + ghostsize
         ymini = min_phys(2)
         
         DO i=1,Nall
            !----------------------------------------------------------------
            !  first those at the south boundary 
            !----------------------------------------------------------------
            IF (xp(2,i).GE.yminf.AND.xp(2,i).LT.ymini) THEN
               xp(2,i)   = xp(2,i) + len_phys(2)
            ENDIF
         ENDDO
         IF (isymm.EQ.0) THEN
            !----------------------------------------------------------------
            !  then the north bc, but only if we are not using symmetry
            !----------------------------------------------------------------
            ymaxf = max_phys(2) + ghostsize
            ymaxi = max_phys(2)
            
            DO i=1,Nall
               IF  (xp(2,i).GT.ymaxi.AND.xp(2,i).LT.ymaxf) THEN
                  xp(2,i)   = xp(2,i) - len_phys(2)
               ENDIF
            ENDDO
         ENDIF 
