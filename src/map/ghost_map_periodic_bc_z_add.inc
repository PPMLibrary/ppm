            !----------------------------------------------------------------
            !  copy periodic ghosts in the z-direction
            !----------------------------------------------------------------
            zminf = min_phys(3) + ghostsize
            zmini = min_phys(3)
            
            DO i=1,Nall
               !-------------------------------------------------------------
               !  first those at the south boundary 
               !-------------------------------------------------------------
               IF (xp(3,i).GE.zminf.AND.xp(3,i).LT.zmini) THEN
                  xp(3,i)   = xp(3,i) + len_phys(3)
               ENDIF
            ENDDO
            IF (isymm.EQ.0) THEN
               !-------------------------------------------------------------
               !  then the north bc, but only if we are not using symmetry
               !-------------------------------------------------------------
               zmaxf = max_phys(3) + ghostsize
               zmaxi = max_phys(3)
               
               DO i=1,Nall
                  IF  (xp(3,i).GT.zmaxi.AND.xp(3,i).LT.zmaxf) THEN
                     xp(3,i)   = xp(3,i) - len_phys(3)
                  ENDIF
               ENDDO
            ENDIF 
