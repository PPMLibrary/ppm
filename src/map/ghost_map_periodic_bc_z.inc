
            !----------------------------------------------------------------
            !  (Re)allocate memory for the periodic ghosts
            !----------------------------------------------------------------
            iopt   = ppm_param_alloc_grow_preserve
            ldu(1) = ppm_dim
            ldu(2) = 2*nghostplus
            CALL ppm_alloc(xt,ldu,iopt,info) 
            CALL ppm_alloc(xt_offset,ldu,iopt,info) 
            IF (info.NE.0) THEN
               info = ppm_error_fatal
               CALL ppm_error(ppm_err_alloc,'ppm_map_part_ghost_get',     &
  &                'xt',__LINE__,info)
               GOTO 9999
            ENDIF

            ldu(1) = ldu(2)
            CALL ppm_alloc(ighost,ldu,iopt,info) 
            IF (info.NE.0) THEN
                info = ppm_error_fatal
                CALL ppm_error(ppm_err_alloc,'ppm_map_part_ghost_get',     &
  &                 'ighost',__LINE__,info)
                GOTO 9999
            ENDIF

            !----------------------------------------------------------------
            !  copy periodic ghosts in the z-direction
            !----------------------------------------------------------------
            zminf = min_phys(3)
            zmini = min_phys(3) + ghostsize
            
            k     = nghostplus 
            DO i=1,nghostplus
               !-------------------------------------------------------------
               !  first those at the south boundary 
               !-------------------------------------------------------------
               IF (xt(3,i).GE.zminf.AND.xt(3,i).LT.zmini) THEN
                  k         = k + 1
                  ighost(k) = ighost(i)
                  xt(1,k)   = xt(1,i) 
                  xt(2,k)   = xt(2,i)
                  xt(3,k)   = xt(3,i) + len_phys(3)

                  xt_offset(1,k) = 0.0_MK
                  xt_offset(2,k) = 0.0_MK
                  xt_offset(3,k) = len_phys(3)
               ENDIF
            ENDDO
            IF (isymm.EQ.0) THEN
               !-------------------------------------------------------------
               !  then the north bc, but only if we are not using symmetry
               !-------------------------------------------------------------
               zmaxf = max_phys(3)
               zmaxi = max_phys(3) - ghostsize
               
               DO i=1,nghostplus
                  IF  (xt(3,i).GT.zmaxi.AND.xt(3,i).LT.zmaxf) THEN
                     k         = k + 1
                     ighost(k) = ighost(i)
                     xt(1,k)   = xt(1,i)
                     xt(2,k)   = xt(2,i) 
                     xt(3,k)   = xt(3,i) - len_phys(3)
                  ENDIF
               ENDDO
            ENDIF 

            !----------------------------------------------------------------
            !  update the ghost counter
            !----------------------------------------------------------------
            nghostplus = k
