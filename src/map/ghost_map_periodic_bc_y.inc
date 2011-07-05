
         !-------------------------------------------------------------------
         !  (Re)allocate memory for the periodic ghosts
         !-------------------------------------------------------------------
         iopt   = ppm_param_alloc_grow_preserve
         ldu(1) = ppm_dim
         ldu(2) = 2*nghostplus
         CALL ppm_alloc(xt,ldu,iopt,info) 
         CALL ppm_alloc(xt_offset,ldu,iopt,info) 
         IF (info.NE.0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_alloc,'ppm_map_part_ghost_get',     &
  &             'xt',__LINE__,info)
            GOTO 9999
         ENDIF

         ldu(1) = ldu(2)
         CALL ppm_alloc(ighost,ldu,iopt,info) 
         IF (info.NE.0) THEN
             info = ppm_error_fatal
             CALL ppm_error(ppm_err_alloc,'ppm_map_part_ghost_get',     &
  &              'ighost',__LINE__,info)
             GOTO 9999
         ENDIF

         !-------------------------------------------------------------------
         !  copy periodic ghosts in the y-direction
         !-------------------------------------------------------------------
         yminf = min_phys(2)
         ymini = min_phys(2) + ghostsize
         
         k     = nghostplus
         DO i=1,nghostplus
            !----------------------------------------------------------------
            !  first those at the south boundary 
            !----------------------------------------------------------------
            IF (xt(2,i).GE.yminf.AND.xt(2,i).LT.ymini) THEN
               k         = k + 1
               ighost(k) = ighost(i)
               xt(1,k)   = xt(1,i) 
               xt(2,k)   = xt(2,i) + len_phys(2)
               xt_offset(1,k) = 0.0_MK
               xt_offset(2,k) = len_phys(2)
               IF (ppm_dim.EQ.3) THEN
                  xt(3,k)   = xt(3,i)
                  xt_offset(3,k) = 0.0_MK
               ENDIF 
            ENDIF
         ENDDO
         IF (isymm.EQ.0) THEN
            !----------------------------------------------------------------
            !  then the north bc, but only if we are not using symmetry
            !----------------------------------------------------------------
            ymaxf = max_phys(2)
            ymaxi = max_phys(2) - ghostsize
            
            DO i=1,nghostplus
               IF  (xt(2,i).GT.ymaxi.AND.xt(2,i).LT.ymaxf) THEN
                  k         = k + 1
                  ighost(k) = ighost(i)
                  xt(1,k)   = xt(1,i)
                  xt(2,k)   = xt(2,i) - len_phys(2)
                  xt_offset(1,k) = 0.0_MK
                  xt_offset(2,k) = -len_phys(2)
                  IF (ppm_dim.EQ.3) THEN
                     xt(3,k)   = xt(3,i)
                     xt_offset(3,k) = 0.0_MK
                  ENDIF 
               ENDIF
            ENDDO
         ENDIF 

         !-------------------------------------------------------------------
         !  update the ghost counter
         !-------------------------------------------------------------------
         nghostplus = k
