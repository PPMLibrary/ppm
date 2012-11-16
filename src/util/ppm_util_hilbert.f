      !--*- f90 -*--------------------------------------------------------------
      !  Subroutine   :                 ppm_hsort
      !-------------------------------------------------------------------------
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
      !-------------------------------------------------------------------------
#if   __KIND == __SINGLE_PRECISION
SUBROUTINE hsort_s(pset,info)
#elif __KIND == __DOUBLE_PRECISION
SUBROUTINE hsort_d(pset,info)
#endif
      !!! This routine provides Hilbert-curve sorting for particle-sets.
      !!!
      !!! This only a stub implementation. Many important features are in this
      !!! version NOT implemented. Only particles that carry no properties and
      !!! have no ghosts can be sorted for now.
      !!!
      !!! Also, no checks are performed on the particles data structure! It is
      !!! your responsibility to check that your particles are in a clean state
      !!! before you reorder them!
      !!! And yes, this routine might actually eat your cat.
      
      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------

      USE ppm_module_particles_typedef
      USE ppm_module_topo_typedef
      USE ppm_module_interfaces
      USE ppm_module_data
      USE ppm_module_topo
      USE ppm_module_alloc
      USE ppm_module_error
      USE iso_c_binding
      
      implicit none
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#ifdef __MPI
      INCLUDE 'mpif.h'
#endif

#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
    
      ! arguments
      
#if   __KIND == __SINGLE_PRECISION
      TYPE(ppm_t_particles_s)         :: pset
#elif __KIND == __DOUBLE_PRECISION
      TYPE(ppm_t_particles_d)         :: pset
#endif
      INTEGER,            INTENT(OUT)   :: info
      interface
        subroutine hilbert_sort_s(xp,np,mind,maxd,dim,info) bind(c)
        use iso_c_binding
        implicit none
        real(c_float)      :: xp(*)
        integer(c_int), value, intent(in) :: np
        real(c_float), intent(in) :: mind(2)
        real(c_float), intent(in) :: maxd(2)
        integer(c_int), value, intent(in) :: dim
        integer(c_int), intent(out) :: info
        end subroutine hilbert_sort_s
        subroutine hilbert_sort_d(xp,np,mind,maxd,dim,info) bind(c)
        use iso_c_binding
        implicit none
        real(c_double)      :: xp(*)
        integer(c_int), value, intent(in) :: np
        real(c_double), intent(in) :: mind(2)
        real(c_double), intent(in) :: maxd(2)
        integer(c_int), value, intent(in) :: dim
        integer(c_int), intent(out) :: info
        end subroutine hilbert_sort_d
      end interface
      TYPE(ppm_t_topo), POINTER         :: topo => NULL()
      INTEGER                           :: i
      start_subroutine("hilbert_sort")

      topo => ppm_topo(pset%active_topoid)%t


#if   __KIND == __SINGLE_PRECISION
      call hilbert_sort_s(pset%xp(:,1:pset%Npart),pset%Npart,&
      &                 topo%min_physs,topo%max_physs,ppm_dim,info)
#elif __KIND == __DOUBLE_PRECISION
      call hilbert_sort_d(pset%xp(:,1:pset%Npart),pset%Npart,&
      &                 topo%min_physd,topo%max_physd,ppm_dim,info)
#endif

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
      end_subroutine()

#if   __KIND == __SINGLE_PRECISION
END SUBROUTINE hsort_s
#elif __KIND == __DOUBLE_PRECISION
END SUBROUTINE hsort_d
#endif
