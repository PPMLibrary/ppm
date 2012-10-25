      !--*- f90 -*--------------------------------------------------------------
      !  Module       :              ppm_module_data_loadbal
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
     
      MODULE ppm_module_data_loadbal
      !!! This module holds data used by the load balancing routines.
      !!!
      !!! [NOTE]
      !!! The variables declared in this modules should not be directly accessed
      !!! by the user. They are managed interally by the library.
         !----------------------------------------------------------------------
         !  Modules
         !----------------------------------------------------------------------
         USE ppm_module_data, ONLY: ppm_kind_single,ppm_kind_double
         PRIVATE :: ppm_kind_double

         !----------------------------------------------------------------------
         !  Timing and load statistics
         !----------------------------------------------------------------------

         REAL(ppm_kind_double) :: ppm_loadbal_decomp_cost = 0.0_ppm_kind_double
         !!! estimated cost of redefining the topology

         INTEGER               :: ppm_loadbal_dcn = 0
         !!! counter of how many times decomp cost was measured (N)

         REAL(ppm_kind_double) :: ppm_loadbal_runsum = 0.0_ppm_kind_double
         !!! running sum of (maxtime - avgtime)

         REAL(ppm_kind_double) :: ppm_loadbal_old_sar = -1.0_ppm_kind_double
         !!! old value of the SAR function

         REAL(ppm_kind_double) :: ppm_loadbal_deltaold = 0.0_ppm_kind_double
         !!! old value of load imbalance (maxtime - avgtime)
         REAL(ppm_kind_double) :: ppm_loadbal_slpavg = 0.0_ppm_kind_double
         !!! slp average SAR needs
         INTEGER               :: ppm_loadbal_nold = 0
         !!! number of
         INTEGER               :: ppm_loadbal_nslp = 0

         INTEGER               :: ppm_loadbal_sendrank   = -1
         !!! The process to which this process will send some subdomains
         INTEGER               :: ppm_loadbal_recvrank   = -1
         !!! The process which will receive some subdomains
         INTEGER               :: ppm_loadbal_n_sendrank = -1
         !!! Neighboring process' choice to send the subdomain to
         INTEGER               :: ppm_loadbal_n_recvrank = -1
         !!! Neighboring process' choice to receive the subdomain from
         INTEGER               :: ppm_loadbal_npart_old  = 1
         !!! Current number of Npart on this processor
         INTEGER               :: ppm_loadbal_npart_new  = 1
         !!! The netto number of particles after partial mapping. Say, if
         !!! 10 particles leave the processor and 5 come in, the value becomes
         !!! -5. This value is later used for determining an estimate elapsed
         !!! time for the next time step.
         REAL(ppm_kind_single) :: ppm_loadbal_oldtime_s = 0.0_ppm_kind_single
         REAL(ppm_kind_double) :: ppm_loadbal_oldtime_d = 0.0_ppm_kind_double
         INTEGER               :: ppm_loadbal_subnpart = 0
         !!! Number of particles in the subdomain that will be sent/received
         REAL(ppm_kind_single), DIMENSION(:), POINTER :: ppm_loadbal_xps => NULL()
         REAL(ppm_kind_double), DIMENSION(:), POINTER :: ppm_loadbal_xpd => NULL()

         REAL(ppm_kind_single) :: ppm_loadbal_proccosts = -1._ppm_kind_single
         REAL(ppm_kind_double) :: ppm_loadbal_proccostd = -1._ppm_kind_double
         !!! DLB cost of this processor
         REAL(ppm_kind_single) :: ppm_loadbal_totcommcosts = -1._ppm_kind_single
         REAL(ppm_kind_double) :: ppm_loadbal_totcommcostd = -1._ppm_kind_double
         !!! DLB total communication cost given by the time model
         REAL(ppm_kind_single) :: ppm_loadbal_totcompcosts = -1._ppm_kind_single
         REAL(ppm_kind_double) :: ppm_loadbal_totcompcostd = -1._ppm_kind_double
         !!! DLB total computation cost given by the time model
         REAL(ppm_kind_single), DIMENSION(:), POINTER :: ppm_loadbal_subcosts => NULL()
         REAL(ppm_kind_double), DIMENSION(:), POINTER :: ppm_loadbal_subcostd => NULL()
         !!! DLB total cost of all the subs on that processor
         !!! It is given by the time model
         REAL(ppm_kind_single), DIMENSION(:), POINTER :: ppm_loadbal_subcommcosts => NULL()
         REAL(ppm_kind_double), DIMENSION(:), POINTER :: ppm_loadbal_subcommcostd => NULL()
         !!! DLB total communication cost of all the subs on that processor
         !!! It is given by the time model
         REAL(ppm_kind_single), DIMENSION(:), POINTER :: ppm_loadbal_subcompcosts => NULL()
         REAL(ppm_kind_double), DIMENSION(:), POINTER :: ppm_loadbal_subcompcostd => NULL()
         !!! DLB total computation COST of all the subs on that processor
         !!! It is given by the time model
         INTEGER                :: ppm_loadbal_comm_part_num = 0
         !!! The number of particles being communicated
!         REAL(ppm_kind_single)  :: ppm_loadbal_comm_part_costs
         REAL(ppm_kind_double)  :: ppm_loadbal_comm_part_cost = -1._ppm_kind_double
         !!! The average cost of a communicated particle
         REAL(ppm_kind_double)  :: ppm_loadbal_comp_part_cost = -1._ppm_kind_double
         !!! The average COST of a real particle
         !!! This is always 1. (normalized)
         !!! 'ppm_loadbal_comm_part_cost' is a multiple of this.
         REAL(ppm_kind_double)  :: ppm_loadbal_comp_ratio = -1._ppm_kind_double
         !!! Normalization ratio of comp time
         REAL(ppm_kind_double)  :: ppm_loadbal_comp_part_time = -1._ppm_kind_double
         !!! The average compute TIME of a real particle
         REAL(ppm_kind_double)  :: ppm_loadbal_comm_part_time = -1._ppm_kind_double
         !!! The average communication TIME of a ghost particle
         INTEGER                :: ppm_loadbal_comm_field_num
         !!! The number of mesh nodes being communicated
         REAL(ppm_kind_single)  :: ppm_loadbal_comm_field_costs
         REAL(ppm_kind_double)  :: ppm_loadbal_comm_field_costd
         !!! The average cost of a communicated mesh node
         REAL(ppm_kind_double)  :: ppm_loadbal_comm_time = 0._ppm_kind_double
         !!! Total communication TIME
         REAL(ppm_kind_double)  :: ppm_loadbal_partmap_time_old = 0._ppm_kind_double
         !!! Old timing for partial mapping
         REAL(ppm_kind_double)  :: ppm_loadbal_partmap_time_new = 0._ppm_kind_double
         !!! New timing for partial mapping
         REAL(ppm_kind_double)  :: ppm_loadbal_ghostmap_time = 0._ppm_kind_double
         !!! Timing for ghost mapping
         INTEGER,DIMENSION(:),POINTER :: ppm_loadbal_isendlist,ppm_loadbal_irecvlist
         !!! Arrays to store the original ppm_isendlist & ppm_irecvlist
         INTEGER                      :: ppm_loadbal_nsendlist
      END MODULE ppm_module_data_loadbal
