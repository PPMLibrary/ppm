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
         !  clear out the new segment of xt_offset to be able to copy the
         !  correct values from the previous loop through xt_offset
         !-------------------------------------------------------------------
         FORALL (i=1:ppm_dim,j=nghostplus+1:2*nghostplus) xt_offset(i,j)=0.0_MK

         !-------------------------------------------------------------------
         !  copy periodic ghosts in the y-direction
         !-------------------------------------------------------------------
         yminf = min_phys(2)
         ymini = min_phys(2) + ghostsize

         k     = nghostplus
         SELECT CASE (ppm_dim)
         CASE (2)
            DO i=1,nghostplus
               !----------------------------------------------------------------
               !  first those at the south boundary
               !----------------------------------------------------------------
               IF (xt(2,i).GE.yminf.AND.xt(2,i).LT.ymini) THEN
                  k               = k + 1
                  ighost(k)       = ighost(i)
                  xt(1,k)         = xt(1,i)
                  xt(2,k)         = xt(2,i) + len_phys(2)
                  xt_offset(1,k)  = xt_offset(1,i)
                  xt_offset(2,k)  = len_phys(2)
                  xt_off_fac(1,k) = xt_off_fac(1,i)
                  xt_off_fac(2,k) = xt_off_fac(2,i)
               ENDIF
            ENDDO
         CASE (3)
            DO i=1,nghostplus
               !----------------------------------------------------------------
               !  first those at the south boundary
               !----------------------------------------------------------------
               IF (xt(2,i).GE.yminf.AND.xt(2,i).LT.ymini) THEN
                  k               = k + 1
                  ighost(k)       = ighost(i)
                  xt(1,k)         = xt(1,i)
                  xt(2,k)         = xt(2,i) + len_phys(2)
                  xt(3,k)         = xt(3,i)
                  xt_offset(1,k)  = xt_offset(1,i)
                  xt_offset(2,k)  = len_phys(2)
                  xt_off_fac(1,k) = xt_off_fac(1,i)
                  xt_off_fac(2,k) = xt_off_fac(2,i)
                  xt_offset(3,k)  = xt_offset(3,i)
                  xt_off_fac(3,k) = xt_off_fac(3,i)
               ENDIF
            ENDDO
         END SELECT

         IF (isymm.EQ.0) THEN
            !----------------------------------------------------------------
            !  then the north bc, but only if we are not using symmetry
            !----------------------------------------------------------------
            ymaxf = max_phys(2)
            ymaxi = max_phys(2) - ghostsize

            SELECT CASE (ppm_dim)
            CASE (2)
               DO i=1,nghostplus
                  IF (xt(2,i).GT.ymaxi.AND.xt(2,i).LT.ymaxf) THEN
                     k               = k + 1
                     ighost(k)       = ighost(i)
                     xt(1,k)         = xt(1,i)
                     xt(2,k)         = xt(2,i) - len_phys(2)
                     xt_offset(1,k)  = xt_offset(1,i)
                     xt_offset(2,k)  = -len_phys(2)
                     xt_off_fac(1,k) = xt_off_fac(1,i)
                     xt_off_fac(2,k) = xt_off_fac(2,i)
                  ENDIF
               ENDDO
            CASE (3)
               DO i=1,nghostplus
                  IF (xt(2,i).GT.ymaxi.AND.xt(2,i).LT.ymaxf) THEN
                     k               = k + 1
                     ighost(k)       = ighost(i)
                     xt(1,k)         = xt(1,i)
                     xt(2,k)         = xt(2,i) - len_phys(2)
                     xt(3,k)         = xt(3,i)
                     xt_offset(1,k)  = xt_offset(1,i)
                     xt_offset(2,k)  = -len_phys(2)
                     xt_off_fac(1,k) = xt_off_fac(1,i)
                     xt_off_fac(2,k) = xt_off_fac(2,i)
                     xt_offset(3,k)  = xt_offset(3,i)
                     xt_off_fac(3,k) = xt_off_fac(3,i)
                  ENDIF
               ENDDO
            END SELECT
         ENDIF

         !-------------------------------------------------------------------
         !  update the ghost counter
         !-------------------------------------------------------------------
         nghostplus = k
