      !--*- f90 -*--------------------------------------------------------------
      !  Module       :                 ppm_module_data_rmsh
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  ETH Zurich
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      MODULE ppm_module_data_rmsh
      !!! This module holds data needed by the (re)meshing routines.
      !!!
      !!! [NOTE]
      !!! The variables declared in this modules should not be directly accessed
      !!! by the user. They are managed interally by the library.
        !-----------------------------------------------------------------------
        !  Includes
        !-----------------------------------------------------------------------
        USE ppm_module_data,ONLY:ppm_kind_single,ppm_kind_double
        !PRIVATE :: ppm_kind_single,ppm_kind_double
        IMPLICIT NONE
        ! time
        REAL(KIND(1.0D0)) :: t0

        !  number of currently implemented kernels
        INTEGER, PARAMETER :: max_defkernels = 4

        !  kernel sizes
        INTEGER, DIMENSION(4) :: ppm_rmsh_kernelsize 
        DATA ppm_rmsh_kernelsize /1,2,2,3/

        !  internal weights
        REAL(ppm_kind_single),   DIMENSION(:,:,:)      , POINTER :: wx1_s
        REAL(ppm_kind_double),   DIMENSION(:,:,:)      , POINTER :: wx1_d
        REAL(ppm_kind_single),   DIMENSION(:,:,:)      , POINTER :: wx2_s
        REAL(ppm_kind_double),   DIMENSION(:,:,:)      , POINTER :: wx2_d
        REAL(ppm_kind_single),   DIMENSION(:,:,:)      , POINTER :: wx3_s
        REAL(ppm_kind_double),   DIMENSION(:,:,:)      , POINTER :: wx3_d    

        !  internal fields
        REAL(ppm_kind_single),   DIMENSION(:,:,:,:)    , POINTER :: tuc_2ds
        REAL(ppm_kind_double),   DIMENSION(:,:,:,:)    , POINTER :: tuc_2dd
        REAL(ppm_kind_single),   DIMENSION(:,:,:,:,:)  , POINTER :: tuc_3ds
        REAL(ppm_kind_double),   DIMENSION(:,:,:,:,:)  , POINTER :: tuc_3dd

        !  internal particle lists
        INTEGER          ,       DIMENSION(:,:)        , POINTER :: list_sub
        INTEGER          ,       DIMENSION(:  )        , POINTER :: store_info


        !  fill data
        INTEGER, DIMENSION(5) :: imasf
        DATA imasf /1,2,3,4,5/


      END MODULE ppm_module_data_rmsh


      
