      !-------------------------------------------------------------------------
      !  Module       :                 ppm_module_data_rmsh
      !-------------------------------------------------------------------------
      !
      !  Purpose      : rmsh module
      !
      !  Remarks      : 
      !
      !  References   : 
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_data_rmsh.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:57  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.3  2006/09/04 18:34:52  pchatela
      !  Fixes and cleanups to make it compile for
      !  release-1-0
      !
      !  Revision 1.2  2004/08/12 10:44:22  michaebe
      !  changed some variable names
      !
      !  Revision 1.1  2004/08/09 16:19:14  michaebe
      !  initial impl.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      MODULE ppm_module_data_rmsh

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


      
