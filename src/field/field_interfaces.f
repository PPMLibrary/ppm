!-------------------------------------------------------------------------
!  Subroutine   :                  field_interfaces
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
!CREATE
SUBROUTINE field_create_(this,lda,name,info)
    IMPORT ppm_t_field_
    CLASS(ppm_t_field_)                      :: this
    INTEGER,                     INTENT(IN) :: lda
    CHARACTER(LEN=*),            INTENT(IN) :: name
    INTEGER,                    INTENT(OUT) :: info
END SUBROUTINE
!DESTROY
SUBROUTINE field_destroy_(this,info)
    IMPORT ppm_t_field_
    CLASS(ppm_t_field_)                 :: this
    INTEGER,               INTENT(OUT) :: info
END SUBROUTINE
!CREATE
SUBROUTINE mesh_discr_info_create_(this,meshID,lda,p_idx,flags,info)
    IMPORT ppm_t_mesh_discr_info_,ppm_mdata_lflags
    CLASS(ppm_t_mesh_discr_info_)           :: this
    INTEGER,                     INTENT(IN) :: meshID
    INTEGER,                     INTENT(IN) :: lda
    INTEGER,                     INTENT(IN) :: p_idx
    LOGICAL,DIMENSION(ppm_mdata_lflags)     :: flags
    INTEGER,                    INTENT(OUT) :: info
END SUBROUTINE
!DESTROY
SUBROUTINE mesh_discr_info_destroy_(this,info)
    IMPORT ppm_t_mesh_discr_info_
    CLASS(ppm_t_mesh_discr_info_)             :: this
    INTEGER,               INTENT(OUT) :: info
END SUBROUTINE

SUBROUTINE field_discretize_on_(this,mesh,info,datatype)
    IMPORT ppm_t_field_,ppm_t_equi_mesh_
    !!! Allocate field on a mesh
    !!! If the field has a procedure for initialization (e.g. an
    !!! initial condition), then the field is also initialized.
    CLASS(ppm_t_field_)                :: this
    CLASS(ppm_t_equi_mesh_)            :: mesh
    !!! mesh onto which this field is to be discretized
    INTEGER,               INTENT(OUT) :: info
    INTEGER, OPTIONAL                  :: datatype
    !!! By default, the type is assumed to be real, double-precision.
END SUBROUTINE 

