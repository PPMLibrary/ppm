  !----------------------------------------------------------
  ! (Re)allocate memory for the symmetry ghosts
  !----------------------------------------------------------
  
  iopt   = ppm_param_alloc_grow_preserve
  ldu(1) = ppm_dim
  ldu(2) = 2*nghostplus
  CALL ppm_alloc(xt,ldu,iopt,info) 
  CALL ppm_alloc(xt_offset,ldu,iopt,info) 
  IF (info.NE.0) THEN
     info = ppm_error_fatal
     CALL ppm_error(ppm_err_alloc,'ppm_map_part_ghost_get', &
          'xt',__LINE__,info)
     GOTO 9999
  ENDIF
  
  ldu(1) = ldu(2)
  CALL ppm_alloc(ighost,ldu,iopt,info) 
  IF (info.NE.0) THEN
     info = ppm_error_fatal
     CALL ppm_error(ppm_err_alloc,'ppm_map_part_ghost_get', &
          'ighost',__LINE__,info)
     GOTO 9999
  ENDIF
  
  !-------------------------------------------------------------------
  !  clear out the new segment of xt_offset to be able to copy the
  !  correct values from the previous loop through xt_offset
  !-------------------------------------------------------------------
  FORALL(i=1:ppm_dim,j=nghostplus+1:2*nghostplus) &
  &      xt_offset(i,j) = 0.0_MK
  
  !----------------------------------------------------------
  ! symmetry ghosts in the y-direction.
  !----------------------------------------------------------
  
  k = nghostplus

  !----------------------------------------------------------
  ! those at the south boundary.
  ! No '=' in the inequality checking.
  !----------------------------------------------------------
  
  yminf = min_phys(2)
#if    __MODE == __HOM
  ymini = min_phys(2) + ghostsize
#elif  __MODE == __INHOM
  ymini = min_phys(2) + ghostsize(2)
#endif 
  
  
  DO i = 1, nghostplus
     IF ( xt(2,i) > yminf .AND. xt(2,i) < ymini ) THEN
        k         = k + 1
        ighost(k) = ighost(i)
        xt(1,k)   = xt(1,i)
        xt(2,k)   = 2.0_MK*yminf - xt(2,i)
        xt_offset(1,k) = xt_offset(1,i)
        xt_offset(2,k) = -2.0_MK*(xt(2,i)-yminf)
        IF ( ppm_dim == 3 ) THEN
           xt(3,k)   = xt(3,i)
           xt_offset(3,k) = xt_offset(3,i)
        END IF
     END IF
  END DO
  
  !----------------------------------------------------------
  ! update the ghost counter.
  !----------------------------------------------------------
  
  nghostplus = k
  
  
