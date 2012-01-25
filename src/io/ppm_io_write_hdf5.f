      !-------------------------------------------------------------------------
      !  Subroutine   :                 ppm_io_write_hdf5
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
#if   __DIM == 0
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE io_write_hdf5_0s(adata,filename,datasetname,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE io_write_hdf5_0d(adata,filename,datasetname,info)
#elif __KIND == __INTEGER
      SUBROUTINE io_write_hdf5_0i(adata,filename,datasetname,info)
#endif
#elif __DIM == 1
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE io_write_hdf5_1s(adata,filename,datasetname,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE io_write_hdf5_1d(adata,filename,datasetname,info)
#elif __KIND == __INTEGER
      SUBROUTINE io_write_hdf5_1i(adata,filename,datasetname,info)
#endif
#elif __DIM == 2
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE io_write_hdf5_2s(adata,filename,datasetname,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE io_write_hdf5_2d(adata,filename,datasetname,info)
#elif __KIND == __INTEGER
      SUBROUTINE io_write_hdf5_2i(adata,filename,datasetname,info)
#endif
#elif __DIM == 3
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE io_write_hdf5_3s(adata,filename,datasetname,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE io_write_hdf5_3d(adata,filename,datasetname,info)
#elif __KIND == __INTEGER
      SUBROUTINE io_write_hdf5_3i(adata,filename,datasetname,info)
#endif
#elif __DIM == 4
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE io_write_hdf5_4s(adata,filename,datasetname,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE io_write_hdf5_4d(adata,filename,datasetname,info)
#elif __KIND == __INTEGER
      SUBROUTINE io_write_hdf5_4i(adata,filename,datasetname,info)
#endif
#elif __DIM == 5
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE io_write_hdf5_5s(adata,filename,datasetname,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE io_write_hdf5_5d(adata,filename,datasetname,info)
#elif __KIND == __INTEGER
      SUBROUTINE io_write_hdf5_5i(adata,filename,datasetname,info)
#endif
#endif


      !!! This routine registers a given variable to a dataset defined
      !!! in a parallel HDF5 file for checkpointing. If such a file
      !!! not exist, it will create one and add a new dataset filled
      !!! with this variable values.

      !!! NOTE: Writing hyperslabs by column is faster in FORTRAN.


      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      IMPLICIT NONE
#if    __KIND == __SINGLE_PRECISION | __KIND == __SINGLE_PRECISION_COMPLEX
      INTEGER, PARAMETER :: MK = ppm_kind_single
#else
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#ifdef __MPI
      INCLUDE 'mpif.h'
#endif
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
#if   __DIM == 0
#if   __KIND == __SINGLE_PRECISION | __KIND == __DOUBLE_PRECISION
      REAL(MK)                           , INTENT(IN) :: adata
#elif __KIND == __INTEGER
      INTEGER                            , INTENT(IN) :: adata
#endif
#elif __DIM == 1
#if   __KIND == __SINGLE_PRECISION | __KIND == __DOUBLE_PRECISION
      REAL(MK)   , DIMENSION(:)          , INTENT(IN) :: adata
#elif __KIND == __INTEGER
      INTEGER    , DIMENSION(:)          , INTENT(IN) :: adata
#endif
#elif __DIM == 2
#if   __KIND == __SINGLE_PRECISION | __KIND == __DOUBLE_PRECISION
      REAL(MK)   , DIMENSION(:,:)        , INTENT(IN) :: adata
#elif __KIND == __INTEGER
      INTEGER    , DIMENSION(:,:)        , INTENT(IN) :: adata
#endif
#elif __DIM == 3
#if   __KIND == __SINGLE_PRECISION | __KIND == __DOUBLE_PRECISION
      REAL(MK)   , DIMENSION(:,:,:)      , INTENT(IN) :: adata
#elif __KIND == __INTEGER
      INTEGER    , DIMENSION(:,:,:)      , INTENT(IN) :: adata
#endif
#elif __DIM == 4
#if   __KIND == __SINGLE_PRECISION | __KIND == __DOUBLE_PRECISION
      REAL(MK)   , DIMENSION(:,:,:,:)    , INTENT(IN) :: adata
#elif __KIND == __INTEGER
      INTEGER    , DIMENSION(:,:,:,:)    , INTENT(IN) :: adata
#endif
#elif __DIM == 5
#if   __KIND == __SINGLE_PRECISION | __KIND == __DOUBLE_PRECISION
      REAL(MK)   , DIMENSION(:,:,:,:,:)  , INTENT(IN) :: adata
#elif __KIND == __INTEGER
      INTEGER    , DIMENSION(:,:,:,:,:)  , INTENT(IN) :: adata
#endif
#endif

      !!! Data vector to be written. 1d array of either single, double,
      !!! single complex, double complex, integer or logical elements.
      !!! Needs to be allocated to proper size before calling this
      !!! routine. Its size determines the number of elements read from the
      !!! file.
      CHARACTER(LEN=*)         , INTENT(IN   ) :: filename
      !!! Name of the to-be-created HDF5 file
      CHARACTER(LEN=*)         , INTENT(IN   ) :: datasetname
      !!! Name of the dataset to-be-created HDF5 file
      INTEGER                         , INTENT(  OUT) :: info
      !!! Return status, 0 on success

      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      REAL(MK)                         :: t0
      LOGICAL                          :: lexist
      !!! TRUE, if an HDF5 exists
      INTEGER(HSIZE_T)                 :: num_of_val
      !!! Number of the elements in the 1D adata array.
      INTEGER(HID_T)                   :: filespace
      !!! Dataspace identifier in file
      INTEGER(HID_T)                   :: file_id
      !!! File identifier
      INTEGER(HID_T)                   :: dset_id
      !!! Dataset identifier
      INTEGER(HID_T)                   :: memspace
      !!! Dataspace identifier in memory
      INTEGER(HID_T)                   :: plist_id,plist_id2
      !!! Property list identifier
      INTEGER(HSIZE_T), DIMENSION(2)   :: count
      INTEGER(HSSIZE_T), DIMENSION(ppm_nproc)  :: offset

      INTEGER(HSIZE_T), DIMENSION(2)   :: dimsfi
      !!! Dataspace identifier in memory
      INTEGER(HID_T)                   :: dataspace_dim
      !!! Dataset dimensions
      CHARACTER(LEN=ppm_char)          :: dsetname,f_name,curr_path
      INTEGER                          :: adata_dim
      !!! number of dimensions in adata


      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_io_write_hdf5',t0,info)
#if   __DIM == 0   |  __DIM == 1
      adata_dim = 1
#elif __DIM == 2
      adata_dim = 2
#elif __DIM == 3
      adata_dim = 3
#elif __DIM == 4
      adata_dim = 4
#elif __DIM == 5
      adata_dim = 5
#endif

#if   __DIM == 0
      dimsfi(1) = 1
      dimsfi(2) = ppm_nproc
#else
      dimsfi(1) = MAXVAL(UBOUND(adata))
      dimsfi(2) = ppm_nproc !*adata_dim
#endif

      dataspace_dim = 2
      dsetname=datasetname
      f_name  =filename
!      WRITE(dsetname,*) trim(adjustl(datasetname))
!      WRITE(f_name,*) trim(adjustl(filename)),'.h5'
!      WRITE(*,*)'filename=',f_name
!      WRITE(*,*)'datasetname=',dsetname
      ! Initialize FORTRAN predefined datatypes for HDF5 writing
      CALL h5open_f(info)
      IF (info .NE. 0) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_ppm_noinit,'ppm_io_write_hdf5',    &
     &       'Cannot initialize FORTRAN predefined datatypes for HDF5',__LINE__,info)
         GOTO 9999
      ENDIF

      ! Setup file access property list with parallel I/O access.

      INQUIRE(file=f_name,exist=lexist)
      ! If the file exist then we do not have to create the file access prop
      IF (lexist) THEN
!          print*,'File exists'
	      CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, info)
	      IF (info .NE. 0) THEN
	         info = ppm_error_error
	         CALL ppm_error(ppm_err_open,'ppm_io_write_hdf5',   &
	 &                'FILE_ACCESS property failed.',__LINE__,info)
	         GOTO 9999
	      ENDIF
!	      print*,'h5pcreate_f::H5P_FILE_ACCESS_F--OK'

	      CALL h5pset_fapl_mpio_f(plist_id, ppm_comm, MPI_INFO_NULL, info)
          IF (info .NE. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_open,'ppm_io_write_hdf5',    &
     &       'h5pset_fapl_mpio_f error.',__LINE__,info)
            GOTO 9999
          ENDIF

!          print*,'h5pset_fapl_mpio_f :: H5P_FILE_ACCESS_F --OK'

	  ELSE ! We need to create both file access and file creation probs
!	      print*,'No file'
	      CALL h5pcreate_f(H5P_FILE_CREATE_F, plist_id2, info)
	      IF (info .NE. 0) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_open,'ppm_io_write_hdf5',   &
     &                'FILE_CREATE property failed',__LINE__,info)
             GOTO 9999
          ENDIF

!          print*,'h5pcreate_f::H5P_FILE_CREATE_F--OK'

          CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, info)
          IF (info .NE. 0) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_open,'ppm_io_write_hdf5',   &
     &                'FILE_ACCESS property failed',__LINE__,info)
             GOTO 9999
          ENDIF

!          print*,'h5pcreate_f::H5P_FILE_ACCESS_F--OK'

          CALL h5pset_fapl_mpio_f(plist_id, ppm_comm, MPI_INFO_NULL, info)
          IF (info .NE. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_open,'ppm_io_write_hdf5',    &
     &       'h5pset_fapl_mpio_f error.',__LINE__,info)
            GOTO 9999
          ENDIF

!          print*,'h5pset_fapl_mpio_f :: H5P_FILE_ACCESS_F --OK'

!          CALL h5pset_fapl_mpio_f(plist_id2, ppm_comm, MPI_INFO_NULL, info)
!          IF (info .NE. 0) THEN
!            info = ppm_error_error
!            CALL ppm_error(ppm_err_open,'ppm_io_write_hdf5',    &
!     &       'h5pset_fapl_mpio_f error.',__LINE__,info)
!            GOTO 9999
!          ENDIF
!
!          print*,'h5pset_fapl_mpio_f :: H5P_FILE_CREATE_F --OK'
      ENDIF

      ! H5Pset_fapl_mpio stores the user-supplied MPI IO parameters comm,
      ! for communicator,and info, for information, in the file access
      ! property list fapl_id. That property list can then be used to create
      ! and/or open a file. H5Pset_fapl_mpio is available only in the
      ! parallel HDF5 library and is not a collective function.
      !
      ! ppm_comm is the MPI communicator to be used for file open, as defined in
      ! MPI_FILE_OPEN of MPI-2. This function makes a duplicate of the
      ! communicator, so modifications to comm after this function call
      ! returns have no effect on the file access property list.
      ! mpi_info is the MPI Info object to be used for file open, as defined
      ! in MPI_FILE_OPEN of MPI-2. This function makes a duplicate copy of
      ! the Info object, so modifications to mpi_info object after this
      ! function call returns will have no effect on the file access
      ! property list.
      ! If the file access property list already contains previously-set
      ! communicator and mpi_info values, those values will be replaced
      ! and the old communicator and Info object will be freed.


      ! Create the file collectively.
      ! H5F_ACC_TRUNC :
      ! If the file already exists, the file is opened with read-write
      ! access, and new data will overwrite any existing data. If the file
      ! does not exist, it is created and opened with read-write access.
      ! H5F_ACC_EXCL  :
      ! If the file already exists, H5Fcreate fails. If the file does
      ! not exist, it is created and opened with read-write access. (Default)
!      print*,'came till here!!!'
      IF (lexist) THEN
        CALL h5fopen_f(f_name,H5F_ACC_RDWR_F,file_id,info, &
     &                   access_prp = plist_id)
      ELSE
        CALL h5fcreate_f(f_name,H5F_ACC_TRUNC_F,file_id,info, &
     &                   creation_prp = plist_id2, access_prp = plist_id)
      ENDIF
      IF (info .NE. 0) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_open,'ppm_io_write_hdf5',    &
     &                     'h5fcreate_f error.',__LINE__,info)
        GOTO 9999
      ENDIF

!      print*,'h5fcreate_f::--OK'

      ! Closes an existing property list class.
      ! Existing property lists of this class will continue to exist,
      ! but new ones are not able to be created.
      CALL h5pclose_f(plist_id, info)
      IF (.NOT.lexist) CALL h5pclose_f(plist_id2, info)
      IF (info .NE. 0) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_open,'ppm_io_write_hdf5',    &
     &       'h5pclose_f error.',__LINE__,info)
         GOTO 9999
      ENDIF

      ! Creates a new simple dataspace and opens it for access,
      ! returning a dataspace identifier (filespace).
      !
      ! rank         is the number of dimensions used in the dataspace.
      ! dimsfi       is a one-dimensional array of size rank specifying
      !                 the size of each dimension of the dataset.
      ! maximum_dims is an array of the same size specifying the upper
      !                 limit on the size of each dimension.
      !                 maximum_dims may be the null pointer, in which
      !                 case the upper limit is the same as current_dims.
      !                 Otherwise, no element of maximum_dims should be
      !                 smaller than the corresponding element of current_dims.
      !                 If an element of maximum_dims is H5S_UNLIMITED,
      !                 the maximum size of the corresponding dimension is unlimited.
      ! Any element of current_dims can be 0 (zero). Note that no data
      ! can be written to a dataset if the size of any dimension of its
      ! current dataspace is 0. This is sometimes a useful initial state
      ! for a dataset.
      ! Any dataset with an unlimited dimension must also be chunked;
      ! see H5Pset_chunk. Similarly, a dataset must be chunked if
      ! current_dims does not equal maximum_dims. The dataspace identifier
      ! returned from this function must be released with H5Sclose or leaks will occur.
!      print*,'h5fcreate_f is fine'
      CALL h5screate_simple_f(dataspace_dim, dimsfi, filespace, info)
      IF (info .NE. 0) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_open,'ppm_io_write_hdf5',    &
     &       'h5screate_simple_f error.',__LINE__,info)
         GOTO 9999
      ENDIF
!      print*,'h5screate_simple_f is fine'
      ! Create the dataset with default properties and links it to a location in the file.
#if   __KIND == __INTEGER
      CALL h5dcreate_f(file_id, dsetname,H5T_NATIVE_INTEGER,filespace, &
     &                 dset_id, info)
#elif __KIND == __SINGLE_PRECISION
      CALL h5dcreate_f(file_id, dsetname,H5T_NATIVE_REAL,filespace, &
     &                 dset_id, info)
#elif __KIND == __DOUBLE_PRECISION
      CALL h5dcreate_f(file_id, dsetname,H5T_NATIVE_DOUBLE,filespace, &
     &                 dset_id, info)
#endif
!      print*,'h5dcreate_f is fine'
      IF (info .NE. 0) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_open,'ppm_io_write_hdf5',    &
     &       'h5dcreate_f error.',__LINE__,info)
         GOTO 9999
      ENDIF
      ! Release and terminate access to a dataspace.
      CALL h5sclose_f(filespace, info)
      IF (info .NE. 0) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_open,'ppm_io_write_hdf5',    &
     &       'h5sclose_f error.',__LINE__,info)
         GOTO 9999
      ENDIF

      ! Each process defines dataset in memory and writes it to the hyperslab
      ! in the file.
      count(1) = dimsfi(1)
      count(2) = dimsfi(2)/ppm_nproc
      offset(1) = 0
      offset(2) = ppm_rank * count(2)
!      print*,'rank',ppm_rank,'dimsfi',dimsfi,'dataspace_dim',dataspace_dim
!      print*,'rank',ppm_rank, 'count',count,'offset',offset
      ! Process creates a new simple dataspace (memspace) and opens it for access,
      ! returning a dataspace identifier (memspace).
      CALL h5screate_simple_f(dataspace_dim, count, memspace, info)
      IF (info .NE. 0) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_open,'ppm_io_write_hdf5',    &
     &       'h5screate_simple_f (local dataspace creation) error.',__LINE__,info)
         GOTO 9999
      ENDIF
      ! Given a dataset id (dset_id), it returns a Dataspace identifier (filespace)
      CALL h5dget_space_f(dset_id, filespace, info)
      IF (info .NE. 0) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_open,'ppm_io_write_hdf5',    &
     &       'h5dget_space_f error.',__LINE__,info)
         GOTO 9999
      ENDIF
      ! Select hyperslab in the file.
      ! H5S_SELECT_SET_F  Replaces existing selection with parameters from this call.
      !                   Overlapping blocks are not supported with this operator.
      ! H5S_SELECT_OR_F   Adds new selection to the existing selection.(Binary OR)

      CALL h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset,count,info)
      IF (info .NE. 0) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_open,'ppm_io_write_hdf5',    &
     &       'h5sselect_hyperslab_f error.',__LINE__,info)
         GOTO 9999
      ENDIF

      ! Initialize data buffer with trivial data.
!      adata = ppm_rank

      ! Create property list for collective dataset write
      ! H5P_DATASET_XFER :: Properties for raw data transfer.
      CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, info)
      IF (info .NE. 0) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_open,'ppm_io_write_hdf5',    &
     &       'h5pcreate_f :: props for data transfer error.',__LINE__,info)
         GOTO 9999
      ENDIF
      ! Sets data transfer mode
      ! H5FD_MPIO_COLLECTIVE :: Use collective I/O access.
      CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, info)
      IF (info .NE. 0) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_open,'ppm_io_write_hdf5',    &
     &       'h5pset_dxpl_mpio_f :: collective I/O access error.',__LINE__,info)
         GOTO 9999
      ENDIF

      ! Write from the buffer to a dataset collectively
      ! dset_id            :: Identifier of the dataset to write to.
      ! H5T_NATIVE_INTEGER :: Memory datatype identifier
      ! adata              :: data buffer (scalar or array)
      ! dimsfi             :: Array to hold corresponding dimension sizes
      !                       of data buffer adata; dimsfi(k) has value
      !                       of the k-th dimension of buffer adata;
      !                       values are ignored if adata is a scalar
#if   __KIND == __INTEGER
      CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, adata, dimsfi,info, &
     &                file_space_id = filespace,mem_space_id=memspace, &
     &                xfer_prp = plist_id)
#elif __KIND == __SINGLE_PRECISION
      CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, adata, dimsfi,info, &
     &                file_space_id = filespace,mem_space_id=memspace, &
     &                xfer_prp = plist_id)
#elif __KIND == __DOUBLE_PRECISION
      CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, adata, dimsfi,info, &
     &                file_space_id = filespace,mem_space_id=memspace, &
     &                xfer_prp = plist_id)
#endif
      IF (info .NE. 0) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_open,'ppm_io_write_hdf5',    &
     &       'h5dwrite_vl_f :: collective file writing error.',__LINE__,info)
         GOTO 9999
      ENDIF

      ! Close dataspaces (local and global)
      CALL h5sclose_f(filespace, info)
      IF (info .NE. 0) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_open,'ppm_io_write_hdf5',    &
     &       'h5sclose_f :: filespace closing error.',__LINE__,info)
         GOTO 9999
      ENDIF
      CALL h5sclose_f(memspace, info)
      IF (info .NE. 0) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_open,'ppm_io_write_hdf5',    &
     &       'h5sclose_f :: memspace closing error.',__LINE__,info)
         GOTO 9999
      ENDIF

      ! Close the dataset and property list.
      CALL h5dclose_f(dset_id, info)
      IF (info .NE. 0) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_open,'ppm_io_write_hdf5',    &
     &       'h5dclose_f :: dataset closing error.',__LINE__,info)
         GOTO 9999
      ENDIF
      CALL h5pclose_f(plist_id, info)
      IF (info .NE. 0) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_open,'ppm_io_write_hdf5',    &
     &       'h5pclose_f :: proplist closing error.',__LINE__,info)
         GOTO 9999
      ENDIF

      ! Close the file
      CALL h5fclose_f(file_id, info)
      ! Close FORTRAN predefined datatypes.
      CALL h5close_f(info)
      IF (info .NE. 0) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_open,'ppm_io_write_hdf5',    &
     &       'h5close_f :: closing error.',__LINE__,info)
         GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_io_write_hdf5',t0,info)
      RETURN
#if   __DIM == 0
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE io_write_hdf5_0s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE io_write_hdf5_0d
#elif __KIND == __INTEGER
      END SUBROUTINE io_write_hdf5_0i
#endif
#elif __DIM == 1
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE io_write_hdf5_1s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE io_write_hdf5_1d
#elif __KIND == __INTEGER
      END SUBROUTINE io_write_hdf5_1i
#endif
#elif __DIM == 2
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE io_write_hdf5_2s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE io_write_hdf5_2d
#elif __KIND == __INTEGER
      END SUBROUTINE io_write_hdf5_2i
#endif
#elif __DIM == 3
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE io_write_hdf5_3s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE io_write_hdf5_3d
#elif __KIND == __INTEGER
      END SUBROUTINE io_write_hdf5_3i
#endif
#elif __DIM == 4
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE io_write_hdf5_4s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE io_write_hdf5_4d
#elif __KIND == __INTEGER
      END SUBROUTINE io_write_hdf5_4i
#endif
#elif __DIM == 5
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE io_write_hdf5_5s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE io_write_hdf5_5d
#elif __KIND == __INTEGER
      END SUBROUTINE io_write_hdf5_5i
#endif
#endif
