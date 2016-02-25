#if __DIM == 2
      SUBROUTINE DTYPE(particles_initialize2d)(Pc,Npart_global,info, &
      &          distrib,topoid,minphys,maxphys,cutoff,name)
#elif __DIM == 3
      SUBROUTINE DTYPE(particles_initialize3d)(Pc,Npart_global,info, &
      &          distrib,topoid,minphys,maxphys,cutoff,name)
#endif
          !!!Set initial particle positions

          DEFINE_MK()
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(DTYPE(ppm_t_particles))                        :: Pc
          !!! Data structure containing the Particle set
          INTEGER,                               INTENT(INOUT) :: Npart_global
          !!! total number of particles that will be initialized
          INTEGER,                               INTENT(  OUT) :: info
          !!! Return status, on success 0.
          !-------------------------------------------------------------------------
          !  Optional arguments
          !-------------------------------------------------------------------------
          INTEGER,                     OPTIONAL, INTENT(IN   ) :: distrib
          !!! type of initial distribution. One of
          !!! ppm_param_part_init_cartesian (default)
          !!! ppm_param_part_init_random
          INTEGER,                     OPTIONAL, INTENT(IN   ) :: topoid
          !!! topology id (used only to get the extent of the physical domain)
          REAL(MK), DIMENSION(ppm_dim),OPTIONAL, INTENT(IN   ) :: minphys
          !!! extent of the physical domain. Only if topoid is not present.
          REAL(MK), DIMENSION(ppm_dim),OPTIONAL, INTENT(IN   ) :: maxphys
          !!! extent of the physical domain. Only if topoid is not present.
          REAL(MK),                    OPTIONAL, INTENT(IN   ) :: cutoff
          !!! cutoff of the particles
          CHARACTER(LEN=*),            OPTIONAL, INTENT(IN   ) :: name
          !!! name for this set of particles
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          TYPE(ppm_t_topo), POINTER :: topo

          REAL(MK), DIMENSION(ppm_dim)      :: min_phys,max_phys,len_phys
          REAL(MK), DIMENSION(:,:), POINTER :: xp
          REAL(MK), DIMENSION(:  ), POINTER :: randnb
          REAL(MK)                          :: y,z,h
          REAL(MK)                          :: hx,hy,hz
          REAL(MK)                          :: shift

          INTEGER :: ip,i,j,k,Npart,iopt
          INTEGER :: nijk(ppm_dim),nijk_global(ppm_dim)
          INTEGER :: remaining_rows
          INTEGER :: distribution

          CHARACTER(LEN = ppm_char) :: filename

          LOGICAL :: lperiodic

          start_subroutine("particles_initialize")

          distribution=MERGE(distrib,ppm_param_part_init_cartesian,PRESENT(distrib))

          !Get boundaries of computational domain
          IF (PRESENT(topoid).AND.(PRESENT(minphys).OR.PRESENT(maxphys))) THEN
             fail("probable conflict of optional arguments. Use topoid OR minphys")
          ENDIF

          lperiodic=.FALSE.

          IF (PRESENT(topoid)) THEN
             topo => ppm_topo(topoid)%t

             SELECT CASE (MK)
             CASE (ppm_kind_single)
                min_phys = topo%min_physs
                max_phys = topo%max_physs

             CASE (ppm_kind_double)
                min_phys = topo%min_physd
                max_phys = topo%max_physd

             END SELECT

             IF (ANY(topo%bcdef(1:2*ppm_dim).EQ.ppm_param_bcdef_periodic)) THEN
                lperiodic=.TRUE.
             ENDIF
          ELSE IF (PRESENT(minphys).AND.PRESENT(maxphys)) THEN
             min_phys = minphys
             max_phys = maxphys
          ELSE
             fail("optional arguments needed to define the domain boundaries",ppm_error=ppm_error_fatal)
          ENDIF

          len_phys=max_phys-min_phys

          check_true(<#MINVAL(len_phys(1:ppm_dim)).GT.0#>,&
          & "Domain length is <= 0 along one dimension. Check input parameters")


          h = (PRODUCT(len_phys)/REAL(Npart_global,MK))**(1./REAL(ppm_dim))
          nijk_global = FLOOR(len_phys/h)
          Npart_global = PRODUCT(nijk_global)
          remaining_rows = MOD(nijk_global(ppm_dim),ppm_nproc)

          !number of particles along x and z
          nijk(1:ppm_dim) = nijk_global(1:ppm_dim)
          !number of particles along y
          nijk(2) = nijk_global(ppm_dim)/ppm_nproc

          !number of particles on this processor
          Npart = PRODUCT(nijk)

          !last proc takes care of the additional rows (remainder)
          IF (ppm_rank.EQ.ppm_nproc-1) THEN
#if   __DIM == 2
             Npart = Npart + remaining_rows * nijk(1)
#elif __DIM == 3
             Npart = Npart + remaining_rows * nijk(1)*nijk(3)
#endif
          ENDIF

          CALL Pc%create(Npart,info,name=name)
          or_fail("Failed to create particle set")

          check_associated(<#Pc%xp#>)

          !use a shortcut, for convenience
          xp => Pc%xp

          !-----------------------------------------------------------------------
          ! set particles
          !-----------------------------------------------------------------------
          ip = 0
          shift = 0._MK
          !shifts positions of the particles. Set to 0.5 to place
          !particles in the middle of each cell, set to 0 to place them in the lower
          !left corner
          h=EPSILON(1._MK)

          SELECT CASE (distribution)
          CASE (ppm_param_part_init_cartesian)
             IF (lperiodic) THEN
#if __DIM == 3
                IF (ALL(topo%bcdef(1:2*ppm_dim:2).EQ.ppm_param_bcdef_periodic)) THEN
                   hz=len_phys(3)/REAL(nijk(3),MK)
                   hy=len_phys(2)/REAL(nijk_global(2),MK)
                   hx=len_phys(1)/REAL(nijk(1),MK)

                   DO k = 1,nijk(3)
                      z = min_phys(3) + hz*(k-1) + shift*hz
                   DO j = 1,nijk(2)
                      y = min_phys(2) + hy*(j-1 + ppm_rank*nijk(2)) + shift*hy
                      DO i = 1,nijk(1)
                         ip = ip + 1
                         xp(1,ip) = min_phys(1) + hx*(i-1) + shift*hx
                         xp(2,ip) = y
                         xp(3,ip) = z

                         ! impose periodic boundaries on X & Y
                         IF (xp(1,ip) .GE. max_phys(1)) xp(1,ip) = xp(1,ip) - len_phys(1)
                         IF (xp(2,ip) .GE. max_phys(2)) xp(2,ip) = xp(2,ip) - len_phys(2)
                         IF (xp(1,ip) .LT. min_phys(1)) xp(1,ip) = xp(1,ip) + len_phys(1)
                         IF (xp(2,ip) .LT. min_phys(2)) xp(2,ip) = xp(2,ip) + len_phys(2)
                         IF (xp(3,ip) .GE. max_phys(3)) xp(3,ip) = xp(3,ip) - len_phys(3)
                         IF (xp(3,ip) .LT. min_phys(3)) xp(3,ip) = xp(3,ip) + len_phys(3)
                      ENDDO
                   ENDDO
                   ENDDO

                   IF (ppm_rank.EQ.(ppm_nproc-1)) THEN
                      DO k = 1,nijk(3)
                         z = min_phys(3) + hz*(k-1) + shift*hz
                      DO j = 1,remaining_rows
                         y = min_phys(2) + hy*(j-1 + ppm_nproc*nijk(2)) + shift*hy
                         DO i = 1,nijk(1)
                            ip = ip + 1
                            xp(1,ip) = min_phys(1) + hx*(i-1) + shift*hx
                            xp(2,ip) = y
                            xp(3,ip) = z

                            ! impose periodic boundaries:
                            IF (xp(1,ip) .GE. max_phys(1)) xp(1,ip) = xp(1,ip) - len_phys(1)
                            IF (xp(2,ip) .GE. max_phys(2)) xp(2,ip) = xp(2,ip) - len_phys(2)
                            IF (xp(1,ip) .LT. min_phys(1)) xp(1,ip) = xp(1,ip) + len_phys(1)
                            IF (xp(2,ip) .LT. min_phys(2)) xp(2,ip) = xp(2,ip) + len_phys(2)
                            IF (xp(3,ip) .GE. max_phys(3)) xp(3,ip) = xp(3,ip) - len_phys(3)
                            IF (xp(3,ip) .LT. min_phys(3)) xp(3,ip) = xp(3,ip) + len_phys(3)
                         ENDDO
                      ENDDO
                      ENDDO
                   ENDIF
                ELSE IF (topo%bcdef(1).EQ.ppm_param_bcdef_periodic.AND. &
                &        topo%bcdef(3).EQ.ppm_param_bcdef_periodic) THEN
                   hz=(len_phys(3)-h)/MERGE(REAL(nijk(3),MK),REAL(nijk(3)-1,MK),shift.GT.0._MK)
                   hy=len_phys(2)/REAL(nijk_global(2),MK)
                   hx=len_phys(1)/REAL(nijk(1),MK)

                   DO k = 1,nijk(3)
                      z = min_phys(3) + hz*(k-1) + shift*hz
                   DO j = 1,nijk(2)
                      y = min_phys(2) + hy*(j-1 + ppm_rank*nijk(2)) + shift*hy
                      DO i = 1,nijk(1)
                         ip = ip + 1
                         xp(1,ip) = min_phys(1) + hx*(i-1) + shift*hx
                         xp(2,ip) = y
                         xp(3,ip) = z

                         ! impose periodic boundaries on X & Y
                         IF (xp(1,ip) .GE. max_phys(1)) xp(1,ip) = xp(1,ip) - len_phys(1)
                         IF (xp(2,ip) .GE. max_phys(2)) xp(2,ip) = xp(2,ip) - len_phys(2)
                         IF (xp(1,ip) .LT. min_phys(1)) xp(1,ip) = xp(1,ip) + len_phys(1)
                         IF (xp(2,ip) .LT. min_phys(2)) xp(2,ip) = xp(2,ip) + len_phys(2)
                      ENDDO
                   ENDDO
                   ENDDO

                   IF (ppm_rank.EQ.(ppm_nproc-1)) THEN
                      DO k = 1,nijk(3)
                         z = min_phys(3) + hz*(k-1) + shift*hz
                      DO j = 1,remaining_rows
                         y = min_phys(2) + hy*(j-1 + ppm_nproc*nijk(2)) + shift*hy
                         DO i = 1,nijk(1)
                            ip = ip + 1
                            xp(1,ip) = min_phys(1) + hx*(i-1) + shift*hx
                            xp(2,ip) = y
                            xp(3,ip) = z

                            ! impose periodic boundaries:
                            IF (xp(1,ip) .GE. max_phys(1)) xp(1,ip) = xp(1,ip) - len_phys(1)
                            IF (xp(2,ip) .GE. max_phys(2)) xp(2,ip) = xp(2,ip) - len_phys(2)
                            IF (xp(1,ip) .LT. min_phys(1)) xp(1,ip) = xp(1,ip) + len_phys(1)
                            IF (xp(2,ip) .LT. min_phys(2)) xp(2,ip) = xp(2,ip) + len_phys(2)
                         ENDDO
                      ENDDO
                      ENDDO
                   ENDIF
                ELSE IF (topo%bcdef(1).EQ.ppm_param_bcdef_periodic.AND. &
                &        topo%bcdef(5).EQ.ppm_param_bcdef_periodic) THEN
                   hz=len_phys(3)/REAL(nijk(3),MK)
                   hy=(len_phys(2)-h)/MERGE(REAL(nijk_global(2),MK),REAL(nijk_global(2)-1,MK),shift.GT.0._MK)
                   hx=len_phys(1)/REAL(nijk(1),MK)

                   DO k = 1,nijk(3)
                      z = min_phys(3) + hz*(k-1) + shift*hz
                   DO j = 1,nijk(2)
                      y = min_phys(2) + hy*(j-1 + ppm_rank*nijk(2)) + shift*hy
                      DO i = 1,nijk(1)
                         ip = ip + 1
                         xp(1,ip) = min_phys(1) + hx*(i-1) + shift*hx
                         xp(2,ip) = y
                         xp(3,ip) = z

                         ! impose periodic boundaries on X & Z
                         IF (xp(1,ip) .GE. max_phys(1)) xp(1,ip) = xp(1,ip) - len_phys(1)
                         IF (xp(1,ip) .LT. min_phys(1)) xp(1,ip) = xp(1,ip) + len_phys(1)
                         IF (xp(3,ip) .GE. max_phys(3)) xp(3,ip) = xp(3,ip) - len_phys(3)
                         IF (xp(3,ip) .LT. min_phys(3)) xp(3,ip) = xp(3,ip) + len_phys(3)
                      ENDDO
                   ENDDO
                   ENDDO

                   IF (ppm_rank.EQ.(ppm_nproc-1)) THEN
                      DO k = 1,nijk(3)
                         z = min_phys(3) + hz*(k-1) + shift*hz
                      DO j = 1,remaining_rows
                         y = min_phys(2) + hy*(j-1 + ppm_nproc*nijk(2)) + shift*hy
                         DO i = 1,nijk(1)
                            ip = ip + 1
                            xp(1,ip) = min_phys(1) + hx*(i-1) + shift*hx
                            xp(2,ip) = y
                            xp(3,ip) = z

                            ! impose periodic boundaries on X & Z
                            IF (xp(1,ip) .GE. max_phys(1)) xp(1,ip) = xp(1,ip) - len_phys(1)
                            IF (xp(1,ip) .LT. min_phys(1)) xp(1,ip) = xp(1,ip) + len_phys(1)
                            IF (xp(3,ip) .GE. max_phys(3)) xp(3,ip) = xp(3,ip) - len_phys(3)
                            IF (xp(3,ip) .LT. min_phys(3)) xp(3,ip) = xp(3,ip) + len_phys(3)
                         ENDDO
                      ENDDO
                      ENDDO
                   ENDIF
                ELSE IF (topo%bcdef(3).EQ.ppm_param_bcdef_periodic.AND. &
                &        topo%bcdef(5).EQ.ppm_param_bcdef_periodic) THEN
                   hz=len_phys(3)/REAL(nijk(3),MK)
                   hy=len_phys(2)/REAL(nijk_global(2),MK)
                   hx=(len_phys(1)-h)/MERGE(REAL(nijk(1),MK),REAL(nijk(1)-1,MK),shift.GT.0._MK)

                   DO k = 1,nijk(3)
                      z = min_phys(3) + hz*(k-1) + shift*hz
                   DO j = 1,nijk(2)
                      y = min_phys(2) + hy*(j-1 + ppm_rank*nijk(2)) + shift*hy
                      DO i = 1,nijk(1)
                         ip = ip + 1
                         xp(1,ip) = min_phys(1) + hx*(i-1) + shift*hx
                         xp(2,ip) = y
                         xp(3,ip) = z

                         ! impose periodic boundaries on Y & Z
                         IF (xp(2,ip) .GE. max_phys(2)) xp(2,ip) = xp(2,ip) - len_phys(2)
                         IF (xp(2,ip) .LT. min_phys(2)) xp(2,ip) = xp(2,ip) + len_phys(2)
                         IF (xp(3,ip) .GE. max_phys(3)) xp(3,ip) = xp(3,ip) - len_phys(3)
                         IF (xp(3,ip) .LT. min_phys(3)) xp(3,ip) = xp(3,ip) + len_phys(3)
                      ENDDO
                   ENDDO
                   ENDDO

                   IF (ppm_rank.EQ.(ppm_nproc-1)) THEN
                      DO k = 1,nijk(3)
                         z = min_phys(3) + hz*(k-1) + shift*hz
                      DO j = 1,remaining_rows
                         y = min_phys(2) + hy*(j-1 + ppm_nproc*nijk(2)) + shift*hy
                         DO i = 1,nijk(1)
                            ip = ip + 1
                            xp(1,ip) = min_phys(1) + hx*(i-1) + shift*hx
                            xp(2,ip) = y
                            xp(3,ip) = z

                            ! impose periodic boundaries:
                            IF (xp(2,ip) .GE. max_phys(2)) xp(2,ip) = xp(2,ip) - len_phys(2)
                            IF (xp(2,ip) .LT. min_phys(2)) xp(2,ip) = xp(2,ip) + len_phys(2)
                            IF (xp(3,ip) .GE. max_phys(3)) xp(3,ip) = xp(3,ip) - len_phys(3)
                            IF (xp(3,ip) .LT. min_phys(3)) xp(3,ip) = xp(3,ip) + len_phys(3)
                         ENDDO
                      ENDDO
                      ENDDO
                   ENDIF
                ELSE IF (topo%bcdef(1).EQ.ppm_param_bcdef_periodic) THEN
                   hz=(len_phys(3)-h)/MERGE(REAL(nijk(3),MK),REAL(nijk(3)-1,MK),shift.GT.0._MK)
                   hy=(len_phys(2)-h)/MERGE(REAL(nijk_global(2),MK),REAL(nijk_global(2)-1,MK),shift.GT.0._MK)
                   hx=len_phys(1)/REAL(nijk(1),MK)

                   DO k = 1,nijk(3)
                      z = min_phys(3) + hz*(k-1) + shift*hz
                   DO j = 1,nijk(2)
                      y = min_phys(2) + hy*(j-1 + ppm_rank*nijk(2)) + shift*hy
                      DO i = 1,nijk(1)
                         ip = ip + 1
                         xp(1,ip) = min_phys(1) + hx*(i-1) + shift*hx
                         xp(2,ip) = y
                         xp(3,ip) = z

                         ! impose periodic boundaries on X
                         IF (xp(1,ip) .GE. max_phys(1)) xp(1,ip) = xp(1,ip) - len_phys(1)
                         IF (xp(1,ip) .LT. min_phys(1)) xp(1,ip) = xp(1,ip) + len_phys(1)
                      ENDDO
                   ENDDO
                   ENDDO

                   IF (ppm_rank.EQ.(ppm_nproc-1)) THEN
                      DO k = 1,nijk(3)
                         z = min_phys(3) + hz*(k-1) + shift*hz
                      DO j = 1,remaining_rows
                         y = min_phys(2) + hy*(j-1 + ppm_nproc*nijk(2)) + shift*hy
                         DO i = 1,nijk(1)
                            ip = ip + 1
                            xp(1,ip) = min_phys(1) + hx*(i-1) + shift*hx
                            xp(2,ip) = y
                            xp(3,ip) = z

                            ! impose periodic boundaries:
                            IF (xp(1,ip) .GE. max_phys(1)) xp(1,ip) = xp(1,ip) - len_phys(1)
                            IF (xp(1,ip) .LT. min_phys(1)) xp(1,ip) = xp(1,ip) + len_phys(1)
                         ENDDO
                      ENDDO
                      ENDDO
                   ENDIF
                ELSE IF (topo%bcdef(3).EQ.ppm_param_bcdef_periodic) THEN
                   hz=(len_phys(3)-h)/MERGE(REAL(nijk(3),MK),REAL(nijk(3)-1,MK),shift.GT.0._MK)
                   hy=len_phys(2)/REAL(nijk_global(2),MK)
                   hx=(len_phys(1)-h)/MERGE(REAL(nijk(1),MK),REAL(nijk(1)-1,MK),shift.GT.0._MK)

                   DO k = 1,nijk(3)
                      z = min_phys(3) + hz*(k-1) + shift*hz
                   DO j = 1,nijk(2)
                      y = min_phys(2) + hy*(j-1 + ppm_rank*nijk(2)) + shift*hy
                      DO i = 1,nijk(1)
                         ip = ip + 1
                         xp(1,ip) = min_phys(1) + hx*(i-1) + shift*hx
                         xp(2,ip) = y
                         xp(3,ip) = z

                         ! impose periodic boundaries on Y
                         IF (xp(2,ip) .GE. max_phys(2)) xp(2,ip) = xp(2,ip) - len_phys(2)
                         IF (xp(2,ip) .LT. min_phys(2)) xp(2,ip) = xp(2,ip) + len_phys(2)
                      ENDDO
                   ENDDO
                   ENDDO

                   IF (ppm_rank.EQ.(ppm_nproc-1)) THEN
                      DO k = 1,nijk(3)
                         z = min_phys(3) + hz*(k-1) + shift*hz
                      DO j = 1,remaining_rows
                         y = min_phys(2) + hy*(j-1 + ppm_nproc*nijk(2)) + shift*hy
                         DO i = 1,nijk(1)
                            ip = ip + 1
                            xp(1,ip) = min_phys(1) + hx*(i-1) + shift*hx
                            xp(2,ip) = y
                            xp(3,ip) = z

                            ! impose periodic boundaries:
                            IF (xp(2,ip) .GE. max_phys(2)) xp(2,ip) = xp(2,ip) - len_phys(2)
                            IF (xp(2,ip) .LT. min_phys(2)) xp(2,ip) = xp(2,ip) + len_phys(2)
                         ENDDO
                      ENDDO
                      ENDDO
                   ENDIF
                ELSE IF (topo%bcdef(5).EQ.ppm_param_bcdef_periodic) THEN
                   hz=len_phys(3)/REAL(nijk(3),MK)
                   hy=(len_phys(2)-h)/MERGE(REAL(nijk_global(2),MK),REAL(nijk_global(2)-1,MK),shift.GT.0._MK)
                   hx=(len_phys(1)-h)/MERGE(REAL(nijk(1),MK),REAL(nijk(1)-1,MK),shift.GT.0._MK)

                   DO k = 1,nijk(3)
                      z = min_phys(3) + hz*(k-1) + shift*hz
                   DO j = 1,nijk(2)
                      y = min_phys(2) + hy*(j-1 + ppm_rank*nijk(2)) + shift*hy
                      DO i = 1,nijk(1)
                         ip = ip + 1
                         xp(1,ip) = min_phys(1) + hx*(i-1) + shift*hx
                         xp(2,ip) = y
                         xp(3,ip) = z

                         ! impose periodic boundaries on Z
                         IF (xp(3,ip) .GE. max_phys(3)) xp(3,ip) = xp(3,ip) - len_phys(3)
                         IF (xp(3,ip) .LT. min_phys(3)) xp(3,ip) = xp(3,ip) + len_phys(3)
                      ENDDO
                   ENDDO
                   ENDDO

                   IF (ppm_rank.EQ.(ppm_nproc-1)) THEN
                      DO k = 1,nijk(3)
                         z = min_phys(3) + hz*(k-1) + shift*hz
                      DO j = 1,remaining_rows
                         y = min_phys(2) + hy*(j-1 + ppm_nproc*nijk(2)) + shift*hy
                         DO i = 1,nijk(1)
                            ip = ip + 1
                            xp(1,ip) = min_phys(1) + hx*(i-1) + shift*hx
                            xp(2,ip) = y
                            xp(3,ip) = z

                            ! impose periodic boundaries:
                            IF (xp(3,ip) .GE. max_phys(3)) xp(3,ip) = xp(3,ip) - len_phys(3)
                            IF (xp(3,ip) .LT. min_phys(3)) xp(3,ip) = xp(3,ip) + len_phys(3)
                         ENDDO
                      ENDDO
                      ENDDO
                   ENDIF
                ENDIF
#else
                IF (ALL(topo%bcdef(1:2*ppm_dim:2).EQ.ppm_param_bcdef_periodic)) THEN
                   hy=len_phys(2)/REAL(nijk_global(2),MK)
                   hx=len_phys(1)/REAL(nijk(1),MK)

                   DO j = 1,nijk(2)
                      y = min_phys(2) + hy*(j-1 + ppm_rank*nijk(2)) + shift*hy
                      DO i = 1,nijk(1)
                         ip = ip + 1
                         xp(1,ip) = min_phys(1) + hx*(i-1) + shift*hx
                         xp(2,ip) = y

                         ! impose periodic boundaries:
                         IF (xp(1,ip) .GE. max_phys(1)) xp(1,ip) = xp(1,ip) - len_phys(1)
                         IF (xp(1,ip) .LT. min_phys(1)) xp(1,ip) = xp(1,ip) + len_phys(1)
                         IF (xp(2,ip) .GE. max_phys(2)) xp(2,ip) = xp(2,ip) - len_phys(2)
                         IF (xp(2,ip) .LT. min_phys(2)) xp(2,ip) = xp(2,ip) + len_phys(2)
                      ENDDO
                   ENDDO

                   IF (ppm_rank.EQ.(ppm_nproc-1)) THEN
                      DO j = 1,remaining_rows
                         y = min_phys(2) + hy*(j-1 + ppm_nproc*nijk(2)) + shift*hy
                         DO i = 1,nijk(1)
                            ip = ip + 1
                            xp(1,ip) = min_phys(1) + hx*(i-1) + shift*hx
                            xp(2,ip) = y

                            ! impose periodic boundaries:
                            IF (xp(1,ip) .GE. max_phys(1)) xp(1,ip) = xp(1,ip) - len_phys(1)
                            IF (xp(1,ip) .LT. min_phys(1)) xp(1,ip) = xp(1,ip) + len_phys(1)
                            IF (xp(2,ip) .GE. max_phys(2)) xp(2,ip) = xp(2,ip) - len_phys(2)
                            IF (xp(2,ip) .LT. min_phys(2)) xp(2,ip) = xp(2,ip) + len_phys(2)
                         ENDDO
                      ENDDO
                   ENDIF
                ELSE IF (topo%bcdef(1).EQ.ppm_param_bcdef_periodic) THEN
                   hy=(len_phys(2)-h)/MERGE(REAL(nijk_global(2),MK),REAL(nijk_global(2)-1,MK),shift.GT.0._MK)
                   hx=len_phys(1)/REAL(nijk(1),MK)

                   DO j = 1,nijk(2)
                      y = min_phys(2) + hy*(j-1 + ppm_rank*nijk(2)) + shift*hy
                      DO i = 1,nijk(1)
                         ip = ip + 1
                         xp(1,ip) = min_phys(1) + hx*(i-1) + shift*hx
                         xp(2,ip) = y

                         ! impose periodic boundaries:
                         IF (xp(1,ip) .GE. max_phys(1)) xp(1,ip) = xp(1,ip) - len_phys(1)
                         IF (xp(1,ip) .LT. min_phys(1)) xp(1,ip) = xp(1,ip) + len_phys(1)
                      ENDDO
                   ENDDO

                   IF (ppm_rank.EQ.(ppm_nproc-1)) THEN
                      DO j = 1,remaining_rows
                         y = min_phys(2) + hy*(j-1 + ppm_nproc*nijk(2)) + shift*hy
                         DO i = 1,nijk(1)
                            ip = ip + 1
                            xp(1,ip) = min_phys(1) + hx*(i-1) + shift*hx
                            xp(2,ip) = y

                            ! impose periodic boundaries:
                            IF (xp(1,ip) .GE. max_phys(1)) xp(1,ip) = xp(1,ip) - len_phys(1)
                            IF (xp(1,ip) .LT. min_phys(1)) xp(1,ip) = xp(1,ip) + len_phys(1)
                         ENDDO
                      ENDDO
                   ENDIF
                ELSE
                   hy=len_phys(2)/REAL(nijk_global(2),MK)
                   hx=(len_phys(1)-h)/MERGE(REAL(nijk(1),MK),REAL(nijk(1)-1,MK),shift.GT.0._MK)

                   DO j = 1,nijk(2)
                      y = min_phys(2) + hy*(j-1 + ppm_rank*nijk(2)) + shift*hy
                      DO i = 1,nijk(1)
                         ip = ip + 1
                         xp(1,ip) = min_phys(1) + hx*(i-1) + shift*hx
                         xp(2,ip) = y

                         ! impose periodic boundaries:
                         IF (xp(2,ip) .GE. max_phys(2)) xp(2,ip) = xp(2,ip) - len_phys(2)
                         IF (xp(2,ip) .LT. min_phys(2)) xp(2,ip) = xp(2,ip) + len_phys(2)
                      ENDDO
                   ENDDO

                   IF (ppm_rank.EQ.(ppm_nproc-1)) THEN
                      DO j = 1,remaining_rows
                         y = min_phys(2) + hy*(j-1 + ppm_nproc*nijk(2)) + shift*hy
                         DO i = 1,nijk(1)
                            ip = ip + 1
                            xp(1,ip) = min_phys(1) + hx*(i-1) + shift*hx
                            xp(2,ip) = y

                            ! impose periodic boundaries:
                            IF (xp(2,ip) .GE. max_phys(2)) xp(2,ip) = xp(2,ip) - len_phys(2)
                            IF (xp(2,ip) .LT. min_phys(2)) xp(2,ip) = xp(2,ip) + len_phys(2)
                         ENDDO
                      ENDDO
                   ENDIF
                ENDIF
#endif
             ELSE
             ! There is no periodic boundary condition
#if __DIM == 3
                hz=(len_phys(3)-h)/MERGE(REAL(nijk(3),MK),REAL(nijk(3)-1,MK),shift.GT.0._MK)
#endif
                hy=(len_phys(2)-h)/MERGE(REAL(nijk_global(2),MK),REAL(nijk_global(2)-1,MK),shift.GT.0._MK)
                hx=(len_phys(1)-h)/MERGE(REAL(nijk(1),MK),REAL(nijk(1)-1,MK),shift.GT.0._MK)

#if __DIM == 3
                DO k = 1,nijk(3)
                   z = min_phys(3) + hz*(k-1) + shift*hz
#endif
                DO j = 1,nijk(2)
                   y = min_phys(2) + hy*(j-1 + ppm_rank*nijk(2)) + shift*hy
                   DO i = 1,nijk(1)
                      ip = ip + 1
                      xp(1,ip) = min_phys(1) + hx*(i-1) + shift*hx
                      xp(2,ip) = y
#if __DIM == 3
                      xp(3,ip) = z
#endif
                   ENDDO
                ENDDO
#if __DIM == 3
                ENDDO
#endif

                IF (ppm_rank.EQ.(ppm_nproc-1)) THEN
#if __DIM == 3
                   DO k = 1,nijk(3)
                      z = min_phys(3) + hz*(k-1) + shift*hz
#endif
                   DO j = 1,remaining_rows
                      y = min_phys(2) + hy*(j-1 + ppm_nproc*nijk(2)) + shift*hy
                      DO i = 1,nijk(1)
                         ip = ip + 1
                         xp(1,ip) = min_phys(1) + hx*(i-1) + shift*hx
                         xp(2,ip) = y
#if __DIM == 3
                         xp(3,ip) = z
#endif
                      ENDDO
                   ENDDO
#if __DIM == 3
                   ENDDO
#endif
                ENDIF
             ENDIF !lperiodic
             Pc%flags(ppm_part_cartesian) = .TRUE.

          CASE (ppm_param_part_init_random)
             NULLIFY(randnb)
             iopt = ppm_param_alloc_fit
#ifdef same_random_sequence_nproc
             ldc(1) = ppm_dim*Npart_global
#else
             ldc(1) = ppm_dim*Npart
#endif
             CALL ppm_alloc(randnb,ldc(1:1),iopt,info)
             or_fail_alloc("randnb")

             IF (.NOT.ASSOCIATED(ppm_particles_seed)) THEN
                CALL RANDOM_SEED(SIZE=ppm_particles_seedsize)

                ldc(1) = ppm_particles_seedsize
                CALL ppm_alloc(ppm_particles_seed,ldc(1:1),iopt,info)
                or_fail_alloc("ppm_particles_seed")

                DO i=1,ppm_particles_seedsize
                   ppm_particles_seed(i)=i*i*i*i
                ENDDO

                CALL RANDOM_SEED(PUT=ppm_particles_seed)
             ENDIF
             CALL RANDOM_NUMBER(randnb)

#if __DIM == 3
             hz=len_phys(3)/REAL(nijk(3),MK)
#endif
             hy=len_phys(2)/REAL(nijk_global(2),MK)
             hx=len_phys(1)/REAL(nijk(1),MK)

#if __DIM == 3
             DO k = 1,nijk(3)
                z = min_phys(3) + hz*(k-1) + shift*hz
#endif
             DO j = 1,nijk(2)
                y = min_phys(2) + hy*(j-1 + ppm_rank*nijk(2))  + shift* hy

                DO i = 1,nijk(1)
                   ip = ip + 1
                   ! uniformly random in cells
#ifdef same_random_sequence_nproc
                   xp(1,ip) = min_phys(1) + hx*(i-1) + shift + &
                   &          randnb(ppm_dim*ppm_rank*PRODUCT(nijk)+ppm_dim*ip - 1)*hx
                   xp(2,ip) = y + randnb(ppm_dim*ppm_rank*PRODUCT(nijk)+ppm_dim*ip    )*hx
#if __DIM == 3
                   xp(3,ip) = z + randnb(ppm_dim*ppm_rank*PRODUCT(nijk)+ppm_dim*ip - 2)*hx
#endif
#else
                   xp(1,ip) = min_phys(1) + hx*(i-1) + shift + randnb(ppm_dim*ip - 1)*hx
                   xp(2,ip) = y + randnb(ppm_dim*ip    )*hx
#if __DIM == 3
                   xp(3,ip) = z + randnb(ppm_dim*ip - 2)*hx
#endif
#endif
                   ! impose periodic boundaries:
                   IF (xp(1,ip) .GE. max_phys(1)) xp(1,ip) = xp(1,ip) - len_phys(1)
                   IF (xp(2,ip) .GE. max_phys(2)) xp(2,ip) = xp(2,ip) - len_phys(2)
                   IF (xp(1,ip) .LT. min_phys(1)) xp(1,ip) = xp(1,ip) + len_phys(1)
                   IF (xp(2,ip) .LT. min_phys(2)) xp(2,ip) = xp(2,ip) + len_phys(2)
#if __DIM == 3
                   IF (xp(3,ip) .GE. max_phys(3)) xp(3,ip) = xp(3,ip) - len_phys(3)
                   IF (xp(3,ip) .LT. min_phys(3)) xp(3,ip) = xp(3,ip) + len_phys(3)
#endif
                ENDDO
             ENDDO
#if __DIM == 3
             ENDDO
#endif
             IF (ppm_rank.EQ.(ppm_nproc-1)) THEN
#if __DIM == 3
                DO k = 1,nijk(3)
                   z = min_phys(3) + hz*(k-1) + shift*hz
#endif
                   DO j = 1,remaining_rows
                      y = min_phys(2) + hy*(j-1 + ppm_nproc*nijk(2)) + shift*hy
                      DO i = 1,nijk(1)
                         ip = ip + 1
                         ! uniformly random in cells
#ifdef same_random_sequence_nproc
                         xp(1,ip) = min_phys(1) + hx*(i-1) + shift + &
                         &          randnb(ppm_dim*(ppm_nproc-1)*PRODUCT(nijk)+ppm_dim*ip - 1)*hx
                         xp(2,ip) = y + randnb(ppm_dim*(ppm_nproc-1)*PRODUCT(nijk)+ppm_dim*ip    )*hx
#if __DIM == 3
                         xp(3,ip) = z + randnb(ppm_dim*(ppm_nproc-1)*PRODUCT(nijk)+ppm_dim*ip - 2)*hx
#endif
#else
                         xp(1,ip) = min_phys(1) + hx*(i-1) + shift + randnb(ppm_dim*ip - 1)*hx
                         xp(2,ip) = y + randnb(ppm_dim*ip    )*hx
#if __DIM == 3
                         xp(3,ip) = z + randnb(ppm_dim*ip - 2)*hx
#endif
#endif
                         ! impose periodic boundaries:
                         IF (xp(1,ip) .GE. max_phys(1)) xp(1,ip) = xp(1,ip) - len_phys(1)
                         IF (xp(2,ip) .GE. max_phys(2)) xp(2,ip) = xp(2,ip) - len_phys(2)
                         IF (xp(1,ip) .LT. min_phys(1)) xp(1,ip) = xp(1,ip) + len_phys(1)
                         IF (xp(2,ip) .LT. min_phys(2)) xp(2,ip) = xp(2,ip) + len_phys(2)
#if __DIM == 3
                         IF (xp(3,ip) .GE. max_phys(3)) xp(3,ip) = xp(3,ip) - len_phys(3)
                         IF (xp(3,ip) .LT. min_phys(3)) xp(3,ip) = xp(3,ip) + len_phys(3)
#endif
                      ENDDO
                   ENDDO
#if __DIM == 3
                ENDDO
#endif
             ENDIF
             Pc%flags(ppm_part_cartesian) = .FALSE.

             iopt=ppm_param_dealloc
             CALL ppm_alloc(randnb,ldc,iopt,info)
             or_fail_dealloc("Could not deallocate randnb")

          CASE DEFAULT
             fail("Unknown distribution type (or not yet supported)")

          END SELECT

          xp => NULL()

          ! (global) average interparticle spacing
          Pc%h_avg = (PRODUCT(len_phys)/REAL(Npart_global))**(1./REAL(ppm_dim))
          ! min interparticle spacing (not needed now)
          Pc%h_min = -1._MK

          Pc%flags(ppm_part_areinside) = .TRUE.

          ! set cutoff to a default value
          Pc%ghostlayer=MERGE(cutoff,2.1_MK*Pc%h_avg,PRESENT(cutoff))

          IF (PRESENT(topoid)) THEN
             Pc%active_topoid = topoid
          ENDIF

          ! Even though convenient, this code should not be here
          ! Neighborlists have lots of options and are costly. The user should
          ! have the chance to decide himself where and how to create them.
          !CALL Pc%create_neighlist(Pc,info,name='self',&
          !    skin=0._MK,symmetry=.FALSE.,cutoff=cutoff)
          !    or_fail("failed to create neighbour list")

          !-----------------------------------------------------------------------
          ! Finalize
          !-----------------------------------------------------------------------
          end_subroutine()

#if __DIM == 2
      END SUBROUTINE DTYPE(particles_initialize2d)
#elif __DIM == 3
      END SUBROUTINE DTYPE(particles_initialize3d)
#endif

#undef __DIM
