  !----------------------------------------------------------
  ! (Re)allocate memory for the symmetry ghosts.
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
  ! symmetry ghosts in the z-direction.
  !----------------------------------------------------------
  
  k = nghostplus
  
  !----------------------------------------------------------
  ! those at down boundary.
  ! No '=' in the inequality checking.
  !----------------------------------------------------------
  
  zminf = min_phys(3)
#if    __MODE == __HOM 
  zmini = min_phys(3) + ghostsize
#elif  __MODE == __INHOM
   !----------------------------------------------------------------
   !  For inhomogenous ghostlayers, we use the entire max
   !----------------------------------------------------------------
#if    __KIND == __SINGLE_PRECISION  
  zmini = min_phys(3) + topo%ghostsizes(3)
#else
  zmini = min_phys(3) + topo%ghostsized(3)
#endif
#endif
  
  DO i = 1, nghostplus
     
     IF ( xt(3,i) > zminf .AND. xt(3,i) < zmini ) THEN
        k         = k + 1
        ighost(k) = ighost(i)
        xt(1,k)   = xt(1,i)
        xt(2,k)   = xt(2,i)
        xt(3,k)   = 2.0_MK*zminf - xt(3,i)
        xt_offset(1,k) = xt_offset(1,i)
        xt_offset(2,k) = xt_offset(2,i)
        xt_offset(3,k) = -2.0_MK*(xt(3,i)-zminf)
     END IF
  END DO
  
  !----------------------------------------------------------
  ! update the ghost counter.
  !----------------------------------------------------------
  
  nghostplus = k
  
