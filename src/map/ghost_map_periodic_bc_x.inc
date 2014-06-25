
         !-------------------------------------------------------------------
         !  (Re)allocate memory for the periodic ghosts
         !-------------------------------------------------------------------
         iopt   = ppm_param_alloc_grow_preserve
         ldu(1) = ppm_dim
         ldu(2) = 2*nghostplus
         CALL ppm_alloc(xt,ldu,iopt,info)
         or_fail_alloc('xt',ppm_error=ppm_error_fatal)

         CALL ppm_alloc(xt_offset,ldu,iopt,info)
         or_fail_alloc('xt_offset',ppm_error=ppm_error_fatal)

         CALL ppm_alloc(xt_off_fac,ldu,iopt,info)
         or_fail_alloc('xt_off_fac',ppm_error=ppm_error_fatal)

         ldu(1) = ldu(2)
         CALL ppm_alloc(ighost,ldu,iopt,info)
         or_fail_alloc('ighost',ppm_error=ppm_error_fatal)

         !-------------------------------------------------------------------
         !  copy periodic ghosts in the x-direction
         !-------------------------------------------------------------------
         xminf = min_phys(1)
         xmini = min_phys(1) + ghostsize

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
               xt_off_fac(1,k) = xt_off_fac(1,i)
               xt_off_fac(2,k) = xt_off_fac(2,i)
               IF (ppm_dim.EQ.3) THEN
                  xt(3,k)   = xt(3,i)
                  xt_offset(3,k) = xt_offset(3,i)
                  xt_off_fac(3,k) = xt_off_fac(3,i)
               ENDIF
            ENDIF
         ENDDO
         IF (isymm.EQ.0) THEN
            !----------------------------------------------------------------
            !  then the east bc, but only if we are not using symmetry
            !----------------------------------------------------------------
            xmaxf = max_phys(1)
            xmaxi = max_phys(1) - ghostsize

            DO i=1,nghostplus
               IF  (xt(1,i).GT.xmaxi.AND.xt(1,i).LT.xmaxf) THEN
                  k         = k + 1
                  ighost(k) = ighost(i)
                  xt(1,k)   = xt(1,i) - len_phys(1)
                  xt(2,k)   = xt(2,i)
                  xt_offset(1,k) = -len_phys(1)
                  xt_offset(2,k) = xt_offset(2,i)
                  xt_off_fac(1,k) = xt_off_fac(1,i)
                  xt_off_fac(2,k) = xt_off_fac(2,i)
                  IF (ppm_dim.EQ.3) THEN
                     xt(3,k)   = xt(3,i)
                     xt_offset(3,k) = xt_offset(3,i)
                     xt_off_fac(3,k) = xt_off_fac(3,i)
                  ENDIF
               ENDIF
            ENDDO
         ENDIF

         !-------------------------------------------------------------------
         !  update the ghost counter
         !-------------------------------------------------------------------
         nghostplus = k
