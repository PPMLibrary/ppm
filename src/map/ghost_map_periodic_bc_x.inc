
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
         !  copy periodic ghosts in the x-direction
         !-------------------------------------------------------------------
         xminf = min_phys(1)
#if    __MODE == __HOM
         xmini = min_phys(1) + ghostsize
#elif  __MODE == __INHOM
            !----------------------------------------------------------------
            !  For inhomogenous ghostlayers, we use the entire max
            !----------------------------------------------------------------
#if    __KIND == __SINGLE_PRECISION
            xmini = min_phys(1) + topo%ghostsizes(1)
#else
            xmini = min_phys(1) + topo%ghostsized(1)
#endif            
#endif 


         k     = nghostplus
         DO i=1,nghostplus
            !----------------------------------------------------------------
            !  first those at the west boundary 
            !----------------------------------------------------------------
            IF (xt(1,i).GE.xminf.AND.xt(1,i).LT.xmini) THEN
               k         = k + 1
               ighost(k) = ighost(i)
               xt(1,k)   = xt(1,i) + len_phys(1)
               xt(2,k)   = xt(2,i)
               xt_offset(1,k) = len_phys(1)
               xt_offset(2,k) = xt_offset(2,i)
               IF (ppm_dim.EQ.3) THEN
                  xt(3,k)   = xt(3,i)
                  xt_offset(3,k) = xt_offset(3,i)
               ENDIF 
            ENDIF
         ENDDO
         IF (isymm.EQ.0) THEN
            !----------------------------------------------------------------
            !  then the east bc, but only if we are not using symmetry
            !----------------------------------------------------------------
            xmaxf = max_phys(1)
#if    __MODE == __HOM
            xmaxi = max_phys(1) - ghostsize
#elif  __MODE == __INHOM
            !----------------------------------------------------------------
            !  For inhomogenous ghostlayers, we use the entire max
            !----------------------------------------------------------------
#if    __KIND == __SINGLE_PRECISION
            xmaxi = max_phys(1) - topo%ghostsizes(1)
#else
            xmaxi = max_phys(1) - topo%ghostsized(1)
#endif
            
#endif 
            
            DO i=1,nghostplus
               IF  (xt(1,i).GT.xmaxi.AND.xt(1,i).LT.xmaxf) THEN
                  k         = k + 1
                  ighost(k) = ighost(i)
                  xt(1,k)   = xt(1,i) - len_phys(1)
                  xt(2,k)   = xt(2,i)
                  xt_offset(1,k) = -len_phys(1)
                  xt_offset(2,k) = xt_offset(2,i)
                  IF (ppm_dim.EQ.3) THEN
                     xt(3,k)   = xt(3,i)
                     xt_offset(3,k) = xt_offset(3,i)
                  ENDIF 
               ENDIF
            ENDDO
         ENDIF 

         !-------------------------------------------------------------------
         !  update the ghost counter
         !-------------------------------------------------------------------
         nghostplus = k
