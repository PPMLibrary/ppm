      !------------------------------------------------------------------------!
      !     Subroutine   :                 ppm_interp_to_mesh_all
      !------------------------------------------------------------------------!
      ! Copyright (c) 2012 CSE Lab (ETH Zurich), MOSAIC Group (ETH Zurich),
      !                    Center for Fluid Dynamics (DTU)
      !
      !
      ! This file is part of the Parallel Particle Mesh Library (PPM).
      !
      ! PPM is free software: you can redistribute it and/or modify
      ! it under the terms of the GNU Lesser General Public License
      ! as published by the Free Software Foundation, either
      ! version 3 of the License, or (at your option) any later
      ! version.
      !
      ! PPM is distributed in the hope that it will be useful,
      ! but WITHOUT ANY WARRANTY; without even the implied warranty of
      ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
      ! GNU General Public License for more details.
      !
      ! You should have received a copy of the GNU General Public License
      ! and the GNU Lesser General Public License along with PPM. If not,
      ! see <http://www.gnu.org/licenses/>.
      !
      ! Parallel Particle Mesh Library (PPM)
      ! ETH Zurich
      ! CH-8092 Zurich, Switzerland
      !------------------------------------------------------------------------!

      SUBROUTINE DTYPE(part_interp_to_mesh_all)(this,Mesh,kernel,info,p2m_bcdef)
      !!! This subroutine interpolate all the fields discretized on the particle set to
      !!! a mesh.  It calls this%interp to carry out particle to mesh interpolation.
      !!!
      !!! Currently 2 interpolation schemes are supported:
      !!!
      !!! * ppm_param_rmsh_kernel_bsp2
      !!! * ppm_param_rmsh_kernel_mp4
      !!!
      !!!
      !!! [TIP]
      !!! There is no need to perform a `ghost_get` before calling this routine
      !!! as the routine calls itself a `ghost_put` to the field after
      !!! interpolating from particles to the field.
      !------------------------------------------------------------------------!
      !  INCLUDES
      !------------------------------------------------------------------------!

      !------------------------------------------------------------------------!
      !  Modules
      !------------------------------------------------------------------------!
      IMPLICIT NONE

      DEFINE_MK()
      !-------------------------------------------------------------------------!
      ! Arguments
      !-------------------------------------------------------------------------!
      CLASS(DTYPE(ppm_t_particles))                   :: this
      CLASS(ppm_t_equi_mesh_)                         :: Mesh
      INTEGER                     ,     INTENT(IN   ) :: kernel
      !!! Choice of the kernel used to compute the weights.
      INTEGER                     ,     INTENT(  OUT) :: info
      !!! Returns status, 0 upon success
      INTEGER, DIMENSION(:  )     , POINTER, OPTIONAL :: p2m_bcdef
     !!! Boundary conditions used for the interpolation routines, they may be
     !!! different than the boundary conditions set in the topology.
     !!! The values in this array can be one of:
     !!!
     !!! - ppm_param_bcdef_symmetry
     !!! - ppm_param_bcdef_antisymmetry
     !-------------------------------------------------------------------------!
     ! Local variables
     !-------------------------------------------------------------------------!
!       REAL(MK) , DIMENSION(:,:)     , POINTER  :: xp => NULL()
!       REAL(MK) , DIMENSION(:)       , POINTER  :: wp_s => NULL()
!       REAL(MK) , DIMENSION(:,:)     , POINTER  :: wp_v => NULL()
      INTEGER                                  :: kernel_support
      INTEGER,  DIMENSION(ppm_dim+2)           :: ldu
      INTEGER                                  :: ip,ndim
      INTEGER                                  :: nb_part
      ! aliases
!       CLASS(ppm_t_subpatch_),POINTER           :: p    => NULL()
      CLASS(ppm_t_main_abstr),POINTER          :: abstr
      CLASS(ppm_t_field_),   POINTER           :: field
!       CLASS(DTYPE(ppm_t_part_prop)_),POINTER   :: prop => NULL()

      start_subroutine("ppm_remesh")

     !-------------------------------------------------------------------------!
     !  This is a hack! Somehow having trouble with constructors in the
     !  ppm_module_data_rmsh module
     !-------------------------------------------------------------------------!
      ppm_rmsh_kernelsize = (/1,2,2,4/)

      ndim = ppm_dim

     !-------------------------------------------------------------------------!
     !  Check arguments
     !-------------------------------------------------------------------------!
      IF (ppm_debug .GT. 0) THEN
        CALL check
        IF (info .NE. 0) GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------!
      !  If there is nothing to do, do nothing
      !-------------------------------------------------------------------------!
      IF(this%Npart.EQ.0) GOTO 9999

      !-------------------------------------------------------------------------
      !  Loop through all the fields that are discretized on the
      !  particle set and interpolate them on the mesh
      !-------------------------------------------------------------------------
      IF (.NOT. ASSOCIATED(this%field_ptr)) THEN
          stdout("No fields are currently discretized on this particle set,remeshing the postions only")
      ELSE
        abstr => this%field_ptr%begin()
        DO WHILE (ASSOCIATED(abstr))
            SELECT TYPE(field => abstr)
            CLASS IS (ppm_t_field_)
            CALL this%interp_to_mesh(Mesh,field,ppm_param_rmsh_kernel_mp4,info)
                or_fail("interp_to_mesh")
            END SELECT
            abstr => this%field_ptr%next()
        ENDDO
      ENDIF

      end_subroutine()
      RETURN
      CONTAINS
      SUBROUTINE check
        IF (.NOT. ppm_initialized) THEN
            fail("Please call ppm_init first!",&
            ppm_err_ppm_noinit,exit_point=8888)
        ENDIF
        IF ((kernel.LT.1).OR.(kernel.GT.4)) THEN
            fail("Wrong kernel definition",&
                ppm_err_ppm_noinit,exit_point=8888)
        END IF
        kernel_support = ppm_rmsh_kernelsize(kernel)*2
        IF(.NOT.((kernel_support.EQ.2).OR.(kernel_support.EQ.4) &
     &               .OR.(kernel_support.EQ.6))) THEN
            fail("Wrong kernel support",ppm_err_argument,exit_point=8888)
        END IF
        IF (this%Npart .GT. 0) THEN
           IF (SIZE(this%xp,2) .LT. this%Npart) THEN
            fail("this%Npart is wrong. Corrupted data structure?",exit_point=8888)
           ENDIF
        ENDIF
 8888   CONTINUE
      END SUBROUTINE check

      END SUBROUTINE DTYPE(part_interp_to_mesh_all)

