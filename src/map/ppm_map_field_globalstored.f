      !-------------------------------------------------------------------------
      !  Subroutine   :                 ppm_map_field_globalstored
      !-------------------------------------------------------------------------
      ! Copyright (c) 2010 CSE Lab (ETH Zurich), MOSAIC Group (ETH Zurich), 
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
      !-------------------------------------------------------------------------

      SUBROUTINE ppm_map_field_globalstored(topoid,target_topoid,            &
     &                                meshid,target_meshid,info)
      !!! This routine maps field data between two topologies
      !!! (which however need compatible meshes defined on
      !!! them) using a global mapping (i.e. every processor
      !!! communicates with every other). Source mesh must be
      !!! on the current field topology. Global lists with
      !!! all mesh blocks that have to be send and/or
      !!! received are built in this routine. Push, pop and
      !!! send will use these lists.
      !!!
      !!! [WARNING]
      !!! This routine assumes that ppm_map_field_init was called before
      !!!
      !!! [WARNING]
      !!! This routine has not been tested, reviewed, or checked. Comments
      !!! and documentation are wrong. This routine might kill your cat or
      !!! worse.
      !!!
      !!! [NOTE]
      !!! The first part of the send/recv lists contains the
      !!! on-processor data.
      !!!
      !!! [CAUTION]
      !!! Side effect: this routine uses the same global
      !!! send/recv buffers, pointers and lists as the
      !!! particle mapping routines (reason: these buffers
      !!! can be large => memory issues). Field and particle
      !!! mappings can therefore never overlap, but one must
      !!! be completed before the other starts.
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_data_mesh
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
      USE ppm_module_write
      USE ppm_module_check_id
      USE ppm_module_mesh_block_intersect
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      INTEGER                 , INTENT(IN   ) :: topoid
      !!! Topology ID of source 
      !!!
      !!! CAUTION: used to be target
      INTEGER                 , INTENT(IN   ) :: target_topoid
      !!! Topology ID of target
      INTEGER                 , INTENT(IN   ) :: meshid
      !!! Mesh ID of source
      INTEGER                 , INTENT(IN   ) :: target_meshid
      !!! Mesh ID of target
      INTEGER                 , INTENT(  OUT) :: info
      !!! Return status, 0 on success
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      INTEGER, DIMENSION(3)            :: ldu
      INTEGER                          :: i,j,idom
      INTEGER                          :: iopt,pdim
      CHARACTER(ppm_char)              :: mesg
      REAL(ppm_kind_double)            :: t0
      LOGICAL                          :: valid
      TYPE(ppm_t_topo),      POINTER   :: topo        => NULL()
      TYPE(ppm_t_topo),      POINTER   :: target_topo => NULL()
      TYPE(ppm_t_equi_mesh), POINTER   :: mesh        => NULL()
      TYPE(ppm_t_equi_mesh), POINTER   :: target_mesh => NULL()
      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_map_field_globalstored',t0,info)
      pdim = ppm_dim

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
        CALL check
        IF (info .NE. 0) GOTO 9999
      ENDIF

      topo => ppm_topo(topoid)%t
      target_topo => ppm_topo(target_topoid)%t

      SELECT TYPE (t => ppm_mesh%vec(meshid))
      TYPE IS (ppm_t_equi_mesh)
          mesh => t
      END SELECT

      SELECT TYPE (t => ppm_mesh%vec(target_meshid))
      TYPE IS (ppm_t_equi_mesh)
          target_mesh => t
      END SELECT



      IF (ppm_buffer_set .GT. 0) THEN
        info = ppm_error_warning
        CALL ppm_error(ppm_err_map_incomp,'ppm_map_field_global_symm',  &
     &      'Buffer was not empty. Possible loss of data!',__LINE__,info)
      ENDIF

      !-------------------------------------------------------------------------
      !  Save the map type for the subsequent calls (used for checks!)
      !-------------------------------------------------------------------------
      ppm_map_type = ppm_param_map_global

      !-------------------------------------------------------------------------
      !  Check if origin and target meshes are compatible (i.e. have the
      !  same number of grid points in the whole comput. domain)
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
        DO i=1,pdim
          IF (mesh%Nm(i) .NE. target_mesh%Nm(i)) THEN
              info = ppm_error_notice
              CALL ppm_error(ppm_err_bad_mesh,'ppm_map_field_globalstored',  &
     &            'source and destination meshes are incompatible',__LINE__,info)
          ENDIF
        ENDDO
      ENDIF

      !-------------------------------------------------------------------------
      !  Point the global mesh sendlists to the precomputed ones
      !-------------------------------------------------------------------------
      iopt   = ppm_param_dealloc
      
      ldu(1) = 1
      IF (ASSOCIATED(ppm_mesh_isendfromsub)) THEN
         CALL ppm_alloc(ppm_mesh_isendfromsub,ldu,iopt,info)
         IF (info .NE. 0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_dealloc,'ppm_map_field_globalstored',     &
     &        'source send sub list PPM_MESH_ISENDFROMSUB',__LINE__,info)
            GOTO 9999
         ENDIF
      END IF
      ppm_mesh_isendfromsub => mesh%mapping%isendfromsub
      
      ldu(2) = 1
      IF (ASSOCIATED(ppm_mesh_isendblkstart)) THEN
         CALL ppm_alloc(ppm_mesh_isendblkstart,ldu,iopt,info)
         IF (info .NE. 0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_dealloc,'ppm_map_field_globalstored',     &
     &        'source send sub list PPM_MESH_ISENDBLKSTART',__LINE__,info)
            GOTO 9999
         ENDIF
      END IF
      ppm_mesh_isendblkstart => mesh%mapping%isendblkstart
      
      IF (ASSOCIATED(ppm_mesh_isendblksize)) THEN
         CALL ppm_alloc(ppm_mesh_isendblksize,ldu,iopt,info)
         IF (info .NE. 0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_dealloc,'ppm_map_field_globalstored',     &
     &        'source send sub list PPM_MESH_ISENDBLKSIZE',__LINE__,info)
            GOTO 9999
         ENDIF
      END IF
      ppm_mesh_isendblksize => mesh%mapping%isendblksize
      
      
      !-------------------------------------------------------------------------
      !  Point the global mesh recvlists to the precomputed ones
      !-------------------------------------------------------------------------
      iopt   = ppm_param_dealloc
      
      ldu(1) = 1
      IF (ASSOCIATED(ppm_mesh_irecvtosub)) THEN
         CALL ppm_alloc(ppm_mesh_irecvtosub,ldu,iopt,info)
         IF (info .NE. 0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_dealloc,'ppm_map_field_globalstored',     &
     &        'source send sub list PPM_MESH_IRECVTOSUB',__LINE__,info)
            GOTO 9999
         ENDIF
      END IF
      ppm_mesh_irecvtosub => mesh%mapping%irecvtosub
      
      ldu(2) = 1
      IF (ASSOCIATED(ppm_mesh_irecvblkstart)) THEN
         CALL ppm_alloc(ppm_mesh_irecvblkstart,ldu,iopt,info)
         IF (info .NE. 0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_dealloc,'ppm_map_field_globalstored',     &
     &        'source send sub list PPM_MESH_IRECVBLKSTART',__LINE__,info)
            GOTO 9999
         ENDIF
      END IF
      ppm_mesh_irecvblkstart => mesh%mapping%irecvblkstart
      
      IF (ASSOCIATED(ppm_mesh_irecvblksize)) THEN
         CALL ppm_alloc(ppm_mesh_irecvblksize,ldu,iopt,info)
         IF (info .NE. 0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_dealloc,'ppm_map_field_globalstored',     &
     &        'source send sub list PPM_MESH_IRECVBLKSIZE',__LINE__,info)
            GOTO 9999
         ENDIF
      END IF
      ppm_mesh_irecvblksize => mesh%mapping%irecvblksize

      !-------------------------------------------------------------------------
      !  Point the global send/recv lists to the precomputed ones
      !-------------------------------------------------------------------------
      iopt   = ppm_param_dealloc
      ldu(1) = 1
      IF (ASSOCIATED(ppm_isendlist)) THEN
         CALL ppm_alloc(ppm_isendlist,ldu,iopt,info)
         IF (info .NE. 0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_dealloc,'ppm_map_field_globalstored',     &
     &        'source send sub list PPM_ISENDLIST',__LINE__,info)
            GOTO 9999
         ENDIF
      END IF
      ppm_isendlist => mesh%mapping%isendlist
      
      IF (ASSOCIATED(ppm_irecvlist)) THEN
         CALL ppm_alloc(ppm_irecvlist,ldu,iopt,info)
         IF (info .NE. 0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_dealloc,'ppm_map_field_globalstored',     &
     &        'source send sub list PPM_IRECVLIST',__LINE__,info)
            GOTO 9999
         ENDIF
      END IF
      ppm_irecvlist => mesh%mapping%irecvlist
      
      IF (ASSOCIATED(ppm_psendbuffer)) THEN
         CALL ppm_alloc(ppm_psendbuffer,ldu,iopt,info)
         IF (info .NE. 0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_dealloc,'ppm_map_field_globalstored',     &
     &        'source send sub list PPM_PSENDBUFFER',__LINE__,info)
            GOTO 9999
         ENDIF
      END IF
      ppm_psendbuffer => mesh%mapping%psendbuffer
      
      IF (ASSOCIATED(ppm_precvbuffer)) THEN
         CALL ppm_alloc(ppm_precvbuffer,ldu,iopt,info)
         IF (info .NE. 0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_dealloc,'ppm_map_field_globalstored',     &
     &        'source send sub list PPM_PRECVBUFFER',__LINE__,info)
            GOTO 9999
         ENDIF
      END IF
      ppm_precvbuffer => mesh%mapping%precvbuffer
      
      ppm_nsendlist = mesh%mapping%nsendlist
      ppm_nrecvlist = mesh%mapping%nrecvlist
      
      !-------------------------------------------------------------------------
      !  Reset the number of buffer entries
      !-------------------------------------------------------------------------
      ppm_buffer_set = 0
      
      !-------------------------------------------------------------------------
      !  Reset the sizes of buffer
      !-------------------------------------------------------------------------
      ppm_nsendbuffer    = 0
      ppm_nrecvbuffer    = 0


      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_map_field_globalstored',t0,info)
      RETURN

      CONTAINS

      SUBROUTINE check
          CALL ppm_check_topoid(target_topoid,valid,info)
          IF (.NOT. valid) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_map_field_globalstored',  &
     &            'target topoid not valid',__LINE__,info)
              GOTO 8888
          ENDIF
          CALL ppm_check_meshid(topoid,meshid,valid,info)
          IF (.NOT. valid) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_map_field_globalstored',  &
     &            'source meshid not valid',__LINE__,info)
              GOTO 8888
          ENDIF
          CALL ppm_check_meshid(target_topoid,target_meshid,valid,info)
          IF (.NOT. valid) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_map_field_globalstored',  &
     &            'destination meshid not valid',__LINE__,info)
              GOTO 8888
          ENDIF
 8888     CONTINUE
      END SUBROUTINE check

      END SUBROUTINE ppm_map_field_globalstored
