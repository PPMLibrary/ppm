         !----------------------------------------------------------
         ! (Re)allocate memory for the symmetry ghosts.
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
       
         !----------------------------------------------------------
         ! symmetry ghosts in the x-direction.
         !----------------------------------------------------------
       
         k = nghostplus
       
         !----------------------------------------------------------
         ! the east bc, even we are using symmetry.
         ! No '=' in the inequality checking.
         !----------------------------------------------------------
       
         xmaxf = max_phys(1)
         xmaxi = max_phys(1) - ghostsize
       
         SELECT CASE (ppm_dim)
         CASE (2)
           DO i = 1, nghostplus
             IF (xt(1,i).GT.xmaxi.AND.xt(1,i).LE.xmaxf) THEN
               k         = k + 1
               ighost(k) = ighost(i)
               xt(1,k)   = 2.0_MK*xmaxf - xt(1,i)
               xt(2,k)   = xt(2,i)
               xt_offset(1,k) = 2.0_MK*xmaxf !2.0_MK*(xmaxf-xt(1,i))
               xt_offset(2,k) = xt_offset(2,i)
               xt_off_fac(1,k) = -1.0_MK
               xt_off_fac(2,k) = xt_off_fac(2,i)
             END IF
           END DO
         CASE (3)
           DO i = 1, nghostplus
             IF (xt(1,i).GT.xmaxi.AND.xt(1,i).LE.xmaxf) THEN
               k         = k + 1
               ighost(k) = ighost(i)
               xt(1,k)   = 2.0_MK*xmaxf - xt(1,i)
               xt(2,k)   = xt(2,i)
               xt(3,k)   = xt(3,i)
               xt_offset(1,k) = 2.0_MK*xmaxf !2.0_MK*(xmaxf-xt(1,i))
               xt_offset(2,k) = xt_offset(2,i)
               xt_offset(3,k) = xt_offset(3,i)
               xt_off_fac(1,k) = -1.0_MK
               xt_off_fac(2,k) = xt_off_fac(2,i)
               xt_off_fac(3,k) = xt_off_fac(3,i)
             END IF
           END DO
         END SELECT
       
         !----------------------------------------------------------
         ! update the ghost counter
         !----------------------------------------------------------
       
         nghostplus = k
       
       
