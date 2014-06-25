         !----------------------------------------------------------
         ! (Re)allocate memory for the symmetry ghosts
         !----------------------------------------------------------
       
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
         FORALL (i=1:ppm_dim,j=nghostplus+1:2*nghostplus) xt_offset(i,j) = 0.0_MK
       
         !----------------------------------------------------------
         ! symmetry ghosts in the y-direction.
         !----------------------------------------------------------
       
         k = nghostplus
       
         !----------------------------------------------------------
         ! those at north boundary,
         ! even using symmetry inter-comunication.
         ! No '=' in the inequality checking.
         !----------------------------------------------------------
       
         ymaxf = max_phys(2)
         ymaxi = max_phys(2) - ghostsize
       
         SELECT CASE (ppm_dim)
         CASE (2)
           DO i = 1, nghostplus
             IF (xt(2,i).GT.ymaxi.AND.xt(2,i).LE.ymaxf) THEN
               k         = k + 1
               ighost(k) = ighost(i)
               xt(1,k)   = xt(1,i)
               xt(2,k)   = 2.0_MK*ymaxf - xt(2,i)
               xt_offset(1,k) = xt_offset(1,i)
               xt_offset(2,k) = 2.0_MK*ymaxf !2.0_MK*(ymaxf-xt(2,i))
               xt_off_fac(1,k) = xt_off_fac(1,i)
               xt_off_fac(2,k) = -1.0_MK
             END IF
           END DO
         CASE (3)
           DO i = 1, nghostplus
             IF (xt(2,i).GT.ymaxi.AND.xt(2,i).LE.ymaxf ) THEN
               k         = k + 1
               ighost(k) = ighost(i)
               xt(1,k)   = xt(1,i)
               xt(2,k)   = 2.0_MK*ymaxf - xt(2,i)
               xt(3,k)   = xt(3,i)
               xt_offset(1,k) = xt_offset(1,i)
               xt_offset(2,k) = 2.0_MK*ymaxf !2.0_MK*(ymaxf-xt(2,i))
               xt_offset(3,k) = xt_offset(3,i)
               xt_off_fac(1,k) = xt_off_fac(1,i)
               xt_off_fac(2,k) = -1.0_MK
               xt_off_fac(3,k) = xt_off_fac(3,i)
             END IF
           END DO
         END SELECT
       
         !----------------------------------------------------------
         ! update the ghost counter.
         !----------------------------------------------------------
       
         nghostplus = k
       
       
