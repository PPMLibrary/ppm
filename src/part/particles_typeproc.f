minclude ppm_create_collection_procedures(DTYPE(particles),DTYPE(particles)_)

      SUBROUTINE DTYPE(part_create)(Pc,Npart,info,name)
          !!! create a set of particles
          !!! This allocates the particle positions.
          !-------------------------------------------------------------------------
          !  Modules
          !-------------------------------------------------------------------------
          IMPLICIT NONE

          DEFINE_MK()
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(DTYPE(ppm_t_particles))             :: Pc
          !!! Data structure containing the particles
          INTEGER,                    INTENT(IN   ) :: Npart
          !!! Number of particles
          INTEGER,                    INTENT(  OUT) :: info
          !!! Returns status, 0 upon success.
          !-------------------------------------------------------------------------
          !  Optional arguments
          !-------------------------------------------------------------------------
          CHARACTER(LEN=*), OPTIONAL, INTENT(IN   ) :: name
          !!! give a name to this Particle set
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          INTEGER :: iopt

          start_subroutine("part_create")

          !-----------------------------------------------------------------
          !  Destroy the DS if it already exists
          !-----------------------------------------------------------------
          IF (ASSOCIATED(Pc%xp)) THEN
             CALL Pc%destroy(info)
             or_fail("Pc%destroy")
          ENDIF

          ! Give a default name to this Particle set
          IF (PRESENT(name)) THEN
             Pc%name = ADJUSTL(TRIM(name))
          ELSE
             Pc%name = particles_dflt_partname()
          ENDIF

          !-----------------------------------------------------------------
          !  Allocate memory for the positions
          !-----------------------------------------------------------------
          iopt=ppm_param_alloc_fit
          ldc(1) = ppm_dim
          ldc(2) = Npart
          CALL ppm_alloc(Pc%xp,ldc(1:2),iopt,info)
          or_fail_alloc("Pc%xp")

          Pc%Npart = Npart
          Pc%Mpart = Npart

          Pc%flags(ppm_part_ghosts)       = .FALSE.
          Pc%flags(ppm_part_areinside)    = .FALSE.
          Pc%flags(ppm_part_partial)      = .FALSE.
          Pc%flags(ppm_part_reqput)       = .FALSE.
          Pc%flags(ppm_part_cartesian)    = .FALSE.
          Pc%flags(ppm_part_neighlists)   = .FALSE.
          Pc%flags(ppm_part_global_index) = .FALSE.

          Pc%active_topoid=-1
          ! No active topology yet

          Pc%ghostlayer=0.0_MK
          ! No ghost layer yet

          Pc%isymm=0
          ! NonSymmetric interactions per default

          IF (.NOT.ASSOCIATED(Pc%props)) THEN
             ALLOCATE(DTYPE(ppm_c_part_prop)::Pc%props,STAT=info)
             or_fail_alloc("could not allocate Pc%props")
          ELSE
             fail("property collection is already allocated. Call destroy() first?")
          ENDIF

          IF (.NOT.ASSOCIATED(Pc%neighs)) THEN
             ALLOCATE(DTYPE(ppm_c_neighlist)::Pc%neighs,STAT=info)
             or_fail_alloc("could not allocate Pc%neighs")
          ELSE
             fail("neighbour list collection is already allocated. Call destroy() first?")
          ENDIF

          Pc%NewNpart=0

          !TOCHECK
          !whether this is the correct allocation?
          !yaser
          IF (.NOT.ASSOCIATED(Pc%field_ptr)) THEN
             ALLOCATE(Pc%field_ptr,STAT=info)
             or_fail_alloc("Pc%field_ptr")
          ELSE
             fail("Pc%field_ptr was already associated. Use destroy() first")
          ENDIF

          IF (ASSOCIATED(Pc%ops)) THEN
             fail("collection of operator pointers is already associated. Call destroy() first?")
          ENDIF

          IF (.NOT.ASSOCIATED(Pc%stats)) THEN
             ALLOCATE(DTYPE(particles_stats)::Pc%stats,STAT=info)
             or_fail_alloc("could not allocate Pc%stats")
          ELSE
             fail("stats structure is already allocated. Call destroy() first?")
          ENDIF

          Pc%time =0.0_MK
          Pc%itime=0

          ! Particles have not been Initialized yet
          Pc%h_avg=-1.0_MK
          Pc%h_min=-1.0_MK

          CALL DTYPE(ppm_part)%vpush(Pc,info)
          or_fail("Failed to push the new Particle inside the Particle collection")

          !Creating a global ID for this particle
          Pc%ID = DTYPE(ppm_part)%nb

          end_subroutine()
      END SUBROUTINE DTYPE(part_create)

      SUBROUTINE DTYPE(part_destroy)(Pc,info)
          !!! Deallocate a ppm_t_particles data type
          !-------------------------------------------------------------------------
          !  Modules
          !-------------------------------------------------------------------------
          IMPLICIT NONE

          DEFINE_MK()
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(DTYPE(ppm_t_particles)) :: Pc
          !!! Data structure containing the Pc
          INTEGER,        INTENT(  OUT) :: info
          !!! Returns status, 0 upon success.
          !-------------------------------------------------------------------------
          ! local variables
          !-------------------------------------------------------------------------
          CLASS(ppm_t_main_abstr),  POINTER :: field
          CLASS(ppm_t_discr_info_), POINTER :: dinfo

          INTEGER :: i

          LOGICAL :: lalloc,ldealloc

          start_subroutine("part_destroy")

          ! first deallocate all content of Pc
          dealloc_pointer(Pc%xp)

          dealloc_pointer(Pc%pcost)

          dealloc_pointer(Pc%stats)

          !Deallocate neighbour lists
          destroy_collection_ptr(Pc%neighs)

          !Deallocate operators
          destroy_collection_ptr(Pc%ops)

          !Destroy the bookkeeping entries in the fields that are
          !discretized on this particle
          IF (ASSOCIATED(Pc%field_ptr)) THEN
             field => Pc%field_ptr%begin()
             field_loop: DO WHILE (ASSOCIATED(field))
                SELECT TYPE(field)
                CLASS IS (ppm_t_field_)
                   NULLIFY(dinfo)
                   IF (field%is_discretized_on(Pc,dinfo)) THEN
                      CALL field%discr_info%remove(info,dinfo)
                      or_fail("field%discr_info%remove")
                   ENDIF

                END SELECT
                field => Pc%field_ptr%next()
              ENDDO field_loop
          ENDIF

          !Deallocate pointers to fields
          destroy_collection_ptr(Pc%field_ptr)

          !Deallocate properties
          destroy_collection_ptr(Pc%props)

          CALL DTYPE(ppm_part)%vremove(info,Pc)
          or_fail("could not remove a detroyed object from a collection")

          Pc%ID = 0

          !-------------------------------------------------------------------------
          !  Finalize
          !-------------------------------------------------------------------------
          end_subroutine()
      END SUBROUTINE DTYPE(part_destroy)

#define __DIM 2
#include "part/ppm_particles_initialize.f"
#define __DIM 3
#include "part/ppm_particles_initialize.f"

      !!temporary hack to deal with both 2d and 3d
      SUBROUTINE DTYPE(part_initialize)(Pc,Npart_global,info,&
      &          distrib,topoid,minphys,maxphys,cutoff,name)
          !-----------------------------------------------------------------------
          ! Set initial particle positions
          !-----------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Modules
          !-------------------------------------------------------------------------
          IMPLICIT NONE

          DEFINE_MK()
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(DTYPE(ppm_t_particles))                         :: Pc
          !!! Data structure containing the particles
          INTEGER,                                INTENT(INOUT) :: Npart_global
          !!! total number of particles that will be initialized
          INTEGER,                                INTENT(  OUT) :: info
          !!! Return status, on success 0.
          !-------------------------------------------------------------------------
          !  Optional arguments
          !-------------------------------------------------------------------------
          INTEGER,                      OPTIONAL, INTENT(IN   ) :: distrib
          !!! type of initial distribution. One of
          !!! ppm_param_part_init_cartesian (default)
          !!! ppm_param_part_init_random
          INTEGER,                      OPTIONAL, INTENT(IN   ) :: topoid
          !!! topology id (used only to get the extent of the physical domain)
          REAL(MK), DIMENSION(ppm_dim), OPTIONAL, INTENT(IN   ) :: minphys
          !!! extent of the physical domain. Only if topoid is not present.
          REAL(MK), DIMENSION(ppm_dim), OPTIONAL, INTENT(IN   ) :: maxphys
          !!! extent of the physical domain. Only if topoid is not present.
          REAL(MK),                     OPTIONAL, INTENT(IN   ) :: cutoff
          !!! cutoff of the particles
          CHARACTER(LEN=*),             OPTIONAL, INTENT(IN   ) :: name
          !!! name for this set of particles
          !-------------------------------------------------------------------------
          ! local variables
          !-------------------------------------------------------------------------
          start_subroutine("part_initialize")

          SELECT CASE (ppm_dim)
          CASE (2)
             CALL DTYPE(particles_initialize2d)(Pc,Npart_global,info,&
             &    distrib,topoid,minphys,maxphys,cutoff,name=name)
          CASE DEFAULT
             CALL DTYPE(particles_initialize3d)(Pc,Npart_global,info,&
             &    distrib,topoid,minphys,maxphys,cutoff,name=name)
          END SELECT

          end_subroutine()
      END SUBROUTINE DTYPE(part_initialize)

      !TOCHECK
      !Yaser
      !The original implementation was totally wrong
      !as you remove particle ip the particles number will change
      !will cause deletion of the wrong particles, so I corrected
      !the implementaion by sorting the deleted particles and removing the
      !them in descending order from the last one.
      SUBROUTINE DTYPE(part_del_parts)(Pc,list_del_parts,nb_del,info)
          !!! remove some particles from a Particle set
          !!! WARNING: this implementation is NOT efficient
          !!! if the number of particles to delete is large.
          !-------------------------------------------------------------------------
          !  Modules
          !-------------------------------------------------------------------------
          USE ppm_module_util_qsort
          IMPLICIT NONE

          DEFINE_MK()
          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(DTYPE(ppm_t_particles))        :: Pc
          !!! Data structure containing the particles
          INTEGER, DIMENSION(:), POINTER       :: list_del_parts
          !!! list of particles to be deleted
          INTEGER,               INTENT(IN   ) :: nb_del
          !!! number of particles to be deleted
          INTEGER,               INTENT(  OUT) :: info
          !!! Returns status, 0 upon success.
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          CLASS(DTYPE(ppm_t_part_prop)_), POINTER :: prop

          INTEGER, DIMENSION(:), POINTER :: index_del_parts

          INTEGER :: i,ip,Npart,del_part,lda

          start_subroutine("part_del_parts")

          !-----------------------------------------------------------------
          !  check arguments
          !-----------------------------------------------------------------
          check_associated(<#Pc%xp#>, &
          & 'Pc structure had not been defined. Call allocate first')

          prop => Pc%props%begin()
          DO WHILE (ASSOCIATED(prop))
             IF (.NOT.prop%flags(ppm_ppt_partial)) THEN
                IF (prop%flags(ppm_ppt_map_parts)) THEN
                   fail('property not mapped, data will be lost')
                ENDIF
             ENDIF
             prop => Pc%props%next()
          ENDDO

          !Yaser
          !In order to make this implementation correct I have done
          !sorted the delted particles and will get rid of them in
          !descending order so the order of particles will not
          !deteriorate the deletion of the rest
          NULLIFY(index_del_parts)
          CALL ppm_util_qsort(list_del_parts,index_del_parts,info)
          or_fail("ppm_util_qsort")

          !-----------------------------------------------------------------
          !  Delete particles
          !-----------------------------------------------------------------
          Npart = Pc%Npart

          del_part = 0

          DO i=nb_del,1,-1
             ip = list_del_parts(index_del_parts(i))

             ! copying particles from the end of xp to the index that has
             ! to be removed
             Pc%xp(1:ppm_dim,ip) = Pc%xp(1:ppm_dim,Npart-i+1)

             del_part = del_part + 1

             prop => Pc%props%begin()
             DO WHILE (ASSOCIATED(prop))
                IF (prop%flags(ppm_ppt_partial)) THEN

                   lda = prop%lda

                   IF (lda.GE.2) THEN
                      SELECT CASE (prop%data_type)
                      CASE (ppm_type_int)
                         prop%data_2d_i(1:lda,ip) =prop%data_2d_i(1:lda,Npart-i+1)

                      CASE (ppm_type_longint)
                         prop%data_2d_li(1:lda,ip)=prop%data_2d_li(1:lda,Npart-i+1)

                      CASE (ppm_type_real,ppm_type_real_single)
                         prop%data_2d_r(1:lda,ip) =prop%data_2d_r(1:lda,Npart-i+1)

                      CASE (ppm_type_comp,ppm_type_comp_single)
                         prop%data_2d_c(1:lda,ip) =prop%data_2d_c(1:lda,Npart-i+1)

                      CASE (ppm_type_logical )
                         prop%data_2d_l(1:lda,ip) =prop%data_2d_l(1:lda,Npart-i+1)

                      END SELECT
                   ELSE
                      SELECT CASE (prop%data_type)
                      CASE (ppm_type_int)
                         prop%data_1d_i(ip) =prop%data_1d_i(Npart-i+1)

                      CASE (ppm_type_longint)
                         prop%data_1d_li(ip)=prop%data_1d_li(Npart-i+1)

                      CASE (ppm_type_real,ppm_type_real_single)
                         prop%data_1d_r(ip) =prop%data_1d_r(Npart-i+1)

                      CASE (ppm_type_comp,ppm_type_comp_single)
                         prop%data_1d_c(ip) =prop%data_1d_c(Npart-i+1)

                      CASE (ppm_type_logical )
                         prop%data_1d_l(ip) =prop%data_1d_l(Npart-i+1)

                      END SELECT
                   ENDIF
                   prop => Pc%props%next()
                ENDIF
             ENDDO !WHILE (ASSOCIATED(prop))
          ENDDO !i=nb_del,1,-1

          !New number of particles, after deleting some
          Pc%Npart = Npart - del_part

          CALL ppm_alloc(index_del_parts,ldc,ppm_param_dealloc,info)
          or_fail_dealloc("index_del_parts")

          end_subroutine()
      END SUBROUTINE DTYPE(part_del_parts)

minclude ppm_create_collection_procedures(DTYPE(part_prop),DTYPE(part_prop)_)

      SUBROUTINE DTYPE(prop_create)(prop,datatype,parts,npart,&
      &          lda,name,flags,info,field,zero)
          !!! Constructor for particle property data structure
          !-------------------------------------------------------------------------
          !  Modules
          !-------------------------------------------------------------------------
          IMPLICIT NONE

          DEFINE_MK()

          CLASS(DTYPE(ppm_t_part_prop))                               :: prop

          INTEGER,                                      INTENT(IN   ) :: datatype

          CLASS(ppm_t_discr_kind),     TARGET,          INTENT(IN   ) :: parts

          INTEGER,                                      INTENT(IN   ) :: npart
          INTEGER,                                      INTENT(IN   ) :: lda

          CHARACTER(LEN=*),                             INTENT(IN   ) :: name
          !!! name to this property

          LOGICAL, DIMENSION(ppm_param_length_pptflags),INTENT(IN   ) :: flags

          INTEGER,                                      INTENT(  OUT) :: info

          CLASS(ppm_t_field_),OPTIONAL,TARGET,          INTENT(IN   ) :: field

          LOGICAL,            OPTIONAL,                 INTENT(IN   ) :: zero
          !!! if true, then initialize the data to zero
          !-------------------------------------------------------------------------
          ! local variables
          !-------------------------------------------------------------------------
          INTEGER :: iopt

          LOGICAL :: zero_data

          start_subroutine("prop_create")

          prop%lda       = lda
          prop%data_type = datatype
          IF (PRESENT(field)) THEN
             prop%field_ptr => field
          ENDIF
          prop%name      = name
          prop%flags     = flags

          zero_data = MERGE(zero,.FALSE.,PRESENT(zero))

          prop%discr => parts

          iopt=ppm_param_alloc_grow

          IF (lda.GE.2) THEN
             ldc(1) = lda
             ldc(2) = npart

             SELECT CASE (datatype)
             CASE (ppm_type_int)
                CALL ppm_alloc(prop%data_2d_i,ldc,iopt,info)

             CASE (ppm_type_longint)
                CALL ppm_alloc(prop%data_2d_li,ldc,iopt,info)

             CASE (ppm_type_real,ppm_type_real_single)
                CALL ppm_alloc(prop%data_2d_r,ldc,iopt,info)

             CASE (ppm_type_comp,ppm_type_comp_single)
                CALL ppm_alloc(prop%data_2d_c,ldc,iopt,info)

             CASE (ppm_type_logical)
                CALL ppm_alloc(prop%data_2d_l,ldc,iopt,info)

             CASE DEFAULT
                fail("invalid type for particle property")

             END SELECT
             or_fail_alloc("allocating property failed",ppm_error=ppm_error_fatal)

             IF (zero_data) THEN
                SELECT CASE (datatype)
                CASE (ppm_type_int)
                   prop%data_2d_i(1:lda,1:npart) = 0

                CASE (ppm_type_longint)
                   prop%data_2d_li(1:lda,1:npart) = 0_ppm_kind_int64

                CASE (ppm_type_real,ppm_type_real_single)
                   prop%data_2d_r(1:lda,1:npart) = 0._MK

                CASE (ppm_type_comp,ppm_type_comp_single)
                   prop%data_2d_c(1:lda,1:npart) = 0._MK

                CASE (ppm_type_logical)
                   prop%data_2d_l(1:lda,1:npart) = .FALSE.

                END SELECT
             ENDIF
          ELSE
             ldc(1) = npart

             SELECT CASE (datatype)
             CASE (ppm_type_int)
                CALL ppm_alloc(prop%data_1d_i,ldc,iopt,info)

             CASE (ppm_type_longint)
                CALL ppm_alloc(prop%data_1d_li,ldc,iopt,info)

             CASE (ppm_type_real,ppm_type_real_single)
                CALL ppm_alloc(prop%data_1d_r,ldc,iopt,info)

             CASE (ppm_type_comp,ppm_type_comp_single)
                CALL ppm_alloc(prop%data_1d_c,ldc,iopt,info)

             CASE (ppm_type_logical)
                CALL ppm_alloc(prop%data_1d_l,ldc,iopt,info)

             CASE DEFAULT
                fail("invalid type for particle property")

             END SELECT
             or_fail_alloc("allocating property failed",ppm_error=ppm_error_fatal)

             IF (zero_data) THEN
                SELECT CASE (datatype)
                CASE (ppm_type_int)
                   prop%data_1d_i(1:npart) = 0

                CASE (ppm_type_longint)
                   prop%data_1d_li(1:npart) = 0_ppm_kind_int64

                CASE (ppm_type_real,ppm_type_real_single)
                   prop%data_1d_r(1:npart) = 0._MK

                CASE (ppm_type_comp,ppm_type_comp_single)
                   prop%data_1d_c(1:npart) = 0._MK

                CASE (ppm_type_logical)
                   prop%data_1d_l(1:npart) = .FALSE.

                END SELECT
             ENDIF
          ENDIF

          end_subroutine()
      END SUBROUTINE DTYPE(prop_create)
      !DESTROY ENTRY
      SUBROUTINE DTYPE(prop_destroy)(prop,info)
          IMPLICIT NONE
          CLASS(DTYPE(ppm_t_part_prop)) :: prop
          INTEGER,       INTENT(  OUT)  :: info
          !!! Returns status, 0 upon success.

          start_subroutine("prop_destroy")

          dealloc_pointer(prop%data_1d_i)
          dealloc_pointer(prop%data_2d_i)
          dealloc_pointer(prop%data_1d_li)
          dealloc_pointer(prop%data_2d_li)
          dealloc_pointer(prop%data_1d_r)
          dealloc_pointer(prop%data_2d_r)
          dealloc_pointer(prop%data_1d_c)
          dealloc_pointer(prop%data_2d_c)
          dealloc_pointer(prop%data_1d_l)
          dealloc_pointer(prop%data_2d_l)

          prop%data_type = ppm_type_none
          prop%lda = 0
          prop%flags = .FALSE.
          prop%field_ptr => NULL()
          prop%name = ""

          end_subroutine()
      END SUBROUTINE DTYPE(prop_destroy)

      SUBROUTINE DTYPE(prop_print_info)(prop,info,level,fileunit,propid)
          !-----------------------------------------------------------------------
          ! Print out summary information about this property
          !-----------------------------------------------------------------------
          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          DEFINE_MK()
          CLASS(DTYPE(ppm_t_part_prop))                          :: prop
          !!! Data structure containing the particles
          INTEGER,                            INTENT(  OUT)      :: info
          !!! Return status, on success 0.
          !-------------------------------------------------------------------------
          !  Optional arguments
          !-------------------------------------------------------------------------
          INTEGER, OPTIONAL,                  INTENT(IN   )      :: level
          !!! indentation level at which to printout the info. Default = 0
          INTEGER, OPTIONAL,                  INTENT(IN   )      :: fileunit
          !!! Already open file unit for printout. Default = stdout
          INTEGER, OPTIONAL,                  INTENT(IN   )      :: propid
          !!! id of this property in the parent struct
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          INTEGER :: lev,fileu,id

          CHARACTER(LEN = ppm_char) :: myformat

          start_subroutine("prop_print_info")

          fileu=MERGE(fileunit,    6,PRESENT(fileunit))
          lev  =MERGE(MAX(level,1),1,PRESENT(level))
          id   =MERGE(propid,      1,PRESENT(propid))

          WRITE(myformat,'(A,I0,A)') '(',4*lev,'X,A,I0,A,A,A,I0)'

          WRITE(fileu,myformat) 'Property ',id,': ', &
          & TRIM(color_print(prop%name,33)),         &
          & ' Type: ',prop%data_type

          lev = lev + 1

          WRITE(myformat,'(A,I0,A)') '(',4*lev,'X,A,I0)'
          WRITE(fileu,myformat) 'lda: ',prop%lda

          WRITE(myformat,'(A,I0,A,I0,A)') '(',4*lev,'X,A,', &
          & ppm_param_length_pptflags,'L)'
          WRITE(fileu,myformat) 'flags: ',prop%flags

          end_subroutine()
      END SUBROUTINE DTYPE(prop_print_info)

      SUBROUTINE DTYPE(part_prop_create)(this,info,field,part_prop, &
      &          discr_data,dtype,name,lda,zero,with_ghosts)
          !!! Adds a property to an existing particle set
          !-------------------------------------------------------------------------
          !  Modules
          !-------------------------------------------------------------------------
          IMPLICIT NONE

          CLASS(DTYPE(ppm_t_particles))                           :: this
          INTEGER,                                  INTENT(  OUT) :: info
          CLASS(ppm_t_field_),            OPTIONAL, INTENT(IN   ) :: field
          CLASS(DTYPE(ppm_t_part_prop)_), OPTIONAL, POINTER       :: part_prop
          !!! Pointer to the ppm_t_part_prop object for that property
          CLASS(ppm_t_discr_data),        OPTIONAL, POINTER       :: discr_data
          !!! Pointer to the ppm_t_discr_data object for that property
          INTEGER,                        OPTIONAL, INTENT(IN   ) :: dtype
          CHARACTER(LEN=*),               OPTIONAL, INTENT(IN   ) :: name
          INTEGER,                        OPTIONAL, INTENT(IN   ) :: lda
          !!! name to this property
          LOGICAL,                        OPTIONAL, INTENT(IN   ) :: zero
          !!! if true, then Initialize the data to zero
          LOGICAL,                        OPTIONAL, INTENT(IN   ) :: with_ghosts
          !!! if true, then allocate with Mpart instead of the default size of Npart
          !-------------------------------------------------------------------------
          ! local variables
          !-------------------------------------------------------------------------
          CLASS(DTYPE(ppm_t_part_prop)_), POINTER :: prop

          DEFINE_MK()
          INTEGER :: lda2,vec_size,npart,i,datatype

          CHARACTER(LEN=ppm_char) :: name2

          LOGICAL, DIMENSION(ppm_param_length_pptflags):: flags

          start_subroutine("particle_prop_create")

          npart = this%Npart
          flags = .FALSE.
          IF (PRESENT(with_ghosts)) THEN
             IF (with_ghosts) THEN
                IF (this%flags(ppm_part_ghosts)) THEN
                   npart = this%Mpart
                   flags(ppm_ppt_ghosts) = .TRUE.
                ELSE
                   fail("trying to init property for ghosts when ghosts are not computed")
                ENDIF
             ENDIF
          ENDIF
          flags(ppm_ppt_partial) = .TRUE.
          flags(ppm_ppt_map_ghosts) = .TRUE.
          flags(ppm_ppt_map_parts) = .TRUE.

          ALLOCATE(DTYPE(ppm_t_part_prop)::prop,STAT=info)
          or_fail_alloc("prop")

          IF (PRESENT(field)) THEN
             lda2 = field%lda
             name2 = field%name
             datatype = field%data_type

             SELECT CASE (PRESENT(zero))
             CASE (.TRUE.)
                ! Create the property
                CALL prop%create(datatype,this,npart,lda2,name2,flags,info,field,zero)
                or_fail("creating property array failed")
             CASE DEFAULT
                ! Create the property
                CALL prop%create(datatype,this,npart,lda2,name2,flags,info,field)
                or_fail("creating property array failed")
             END SELECT
          ELSE
             lda2 = 1
             IF (PRESENT(lda)) THEN
                IF (lda.GE.2) THEN
                   lda2 = lda
                ENDIF
             ENDIF
             IF (PRESENT(name)) THEN
                name2=TRIM(ADJUSTL(name))
             ELSE
                name2="default_ppt_name"
             ENDIF
             IF (PRESENT(dtype)) THEN
                datatype=dtype
             ELSE
                datatype=MERGE(ppm_type_real,ppm_type_real_single,MK.EQ.ppm_kind_double)
             ENDIF

             SELECT CASE (PRESENT(zero))
             CASE (.TRUE.)
                ! Create the property
                CALL prop%create(datatype,this,npart,lda2,name2,flags,info,zero=zero)
                or_fail("creating property array failed")
             CASE DEFAULT
                ! Create the property
                CALL prop%create(datatype,this,npart,lda2,name2,flags,info)
                or_fail("creating property array failed")
             END SELECT
          ENDIF

          IF (PRESENT(part_prop)) THEN
             part_prop => prop
          ENDIF
          IF (PRESENT(discr_data)) THEN
             discr_data => prop
          ENDIF

          CALL this%props%push(prop,info)
          or_fail("pushing new property into collection failed")

          end_subroutine()
      END SUBROUTINE DTYPE(part_prop_create)

      SUBROUTINE DTYPE(part_prop_destroy)(this,prop,info)
          !!! Destroy a property from an existing particle set
          !-------------------------------------------------------------------------
          !  Modules
          !-------------------------------------------------------------------------
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(DTYPE(ppm_t_particles))          :: this
          CLASS(ppm_t_discr_data), INTENT(INOUT) :: prop

          INTEGER,                 INTENT(  OUT) :: info
          !-------------------------------------------------------------------------
          ! local variables
          !-------------------------------------------------------------------------
          start_subroutine("part_prop_destroy")

          check_associated(<#this%props#>,&
          & "Particle set does not contain any discretized data")

          SELECT TYPE(prop)
          CLASS IS (DTYPE(ppm_t_part_prop))
              CALL this%props%remove(info,prop)
              or_fail("could not remove property from its container")

          CLASS DEFAULT
              fail("discretization data has to be of class ppm_t_part_prop")
          END SELECT

          end_subroutine()
      END SUBROUTINE DTYPE(part_prop_destroy)

      SUBROUTINE DTYPE(part_prop_realloc)(Pc,prop,info,with_ghosts,datatype,lda)
          !!! Reallocate the property array to the correct size
          !!! (e.g. if the number of particles has changed or if the type
          !!! of the data changes)
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(DTYPE(ppm_t_particles))                 :: Pc
          CLASS(DTYPE(ppm_t_part_prop)_), INTENT(INOUT) :: prop
          INTEGER,                        INTENT(  OUT) :: info
          LOGICAL, OPTIONAL,              INTENT(IN   ) :: with_ghosts
          !!! if true, then allocate with Mpart instead of the default size of Npart
          INTEGER, OPTIONAL,              INTENT(IN   ) :: datatype
          !!! deallocate the old data array and allocate a new one,
          !!! possibly of a different data type.
          INTEGER, OPTIONAL,              INTENT(IN   ) :: lda
          !!! deallocate the old data array and allocate a new one,
          !!! possibly of a different dimension
          !-------------------------------------------------------------------------
          ! local variables
          !-------------------------------------------------------------------------
          CLASS(ppm_t_field_), POINTER :: field

          INTEGER :: lda2,vec_size,npart,i,dtype

          LOGICAL, DIMENSION(ppm_param_length_pptflags) :: flags

          CHARACTER(LEN=ppm_char) :: name2

          start_subroutine("realloc_prop")

          flags = prop%flags
          name2 = prop%name
          npart = Pc%Npart
          IF (PRESENT(with_ghosts)) THEN
             IF (with_ghosts) THEN
                IF (Pc%flags(ppm_part_ghosts)) THEN
                   npart = Pc%Mpart
                   flags(ppm_ppt_ghosts) = .TRUE.
                ELSE
                   fail("trying to init property for ghosts when ghosts are not computed")
                ENDIF
             ENDIF
          ENDIF
          flags(ppm_ppt_partial) = .TRUE.

          dtype=MERGE(datatype,prop%data_type,PRESENT(datatype))
          lda2=MERGE(lda,prop%lda,PRESENT(lda))

          IF (lda2.NE.prop%lda.OR.dtype.NE.prop%data_type) THEN
             CALL prop%destroy(info)
          ENDIF

          IF (ASSOCIATED(prop%field_ptr)) THEN
             SELECT TYPE(f => prop%field_ptr)
             CLASS IS (ppm_t_field_)
                field => f
             END SELECT

             ! Create the property
             CALL prop%create(dtype,Pc,npart,lda2,name2,flags,info,field)
             or_fail("reallocating property array failed")
          ELSE
             ! Create the property
             CALL prop%create(dtype,Pc,npart,lda2,name2,flags,info)
             or_fail("reallocating property array failed")
          ENDIF

          end_subroutine()
      END SUBROUTINE DTYPE(part_prop_realloc)

      SUBROUTINE DTYPE(part_get_discr)(this,Field,prop,info)
          !!! Returns a pointer to the ppm_t_discr_data object that is
          !!! the discretization of that Field on this Particle set.
          !!! Fails with an error if the Field is not discretized here.
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(DTYPE(ppm_t_particles))          :: this
          CLASS(ppm_t_field_),     TARGET        :: Field
          CLASS(ppm_t_discr_data), POINTER       :: prop

          INTEGER,                 INTENT(  OUT) :: info
          !-------------------------------------------------------------------------
          ! local variables
          !-------------------------------------------------------------------------
          start_subroutine("part_get_discr")

          prop => this%props%begin()
          loop: DO WHILE(ASSOCIATED(prop))
              IF (ASSOCIATED(prop%field_ptr,Field)) THEN
                  EXIT loop
              ENDIF
              prop => this%props%next()
          ENDDO loop

          check_associated(prop,"Failed to find a discretization of this field for this particle set")

          end_subroutine()
      END SUBROUTINE

      SUBROUTINE DTYPE(part_prop_zero)(this,Field,info)
          !!! Reset values of a property to zero
          !!! (a bit unrolled)
          IMPLICIT NONE

          DEFINE_MK()
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(DTYPE(ppm_t_particles))      :: this
          CLASS(ppm_t_field_), TARGET        :: Field

          INTEGER,             INTENT(  OUT) :: info
          !-------------------------------------------------------------------------
          ! local variables
          !-------------------------------------------------------------------------
          INTEGER :: lda

          start_subroutine("part_prop_zero")

          lda = Field%lda
          SELECT CASE(lda)
          CASE (1)
              foreach p in particles(this) with sca_fields(w=Field) prec(DTYPE(prec))
                  w_p = 0._MK
              end foreach
          CASE (2)
              foreach p in particles(this) with vec_fields(w=Field) prec(DTYPE(prec))
                  w_p(1) = 0._MK
                  w_p(2) = 0._MK
              end foreach
          CASE (3)
              foreach p in particles(this) with vec_fields(w=Field) prec(DTYPE(prec))
                  w_p(1) = 0._MK
                  w_p(2) = 0._MK
                  w_p(3) = 0._MK
              end foreach
          CASE (4)
              foreach p in particles(this) with vec_fields(w=Field) prec(DTYPE(prec))
                  w_p(1) = 0._MK
                  w_p(2) = 0._MK
                  w_p(3) = 0._MK
                  w_p(4) = 0._MK
              end foreach
          CASE (5)
              foreach p in particles(this) with vec_fields(w=Field) prec(DTYPE(prec))
                  w_p(1) = 0._MK
                  w_p(2) = 0._MK
                  w_p(3) = 0._MK
                  w_p(4) = 0._MK
                  w_p(5) = 0._MK
              end foreach
          CASE DEFAULT
              foreach p in particles(this) with vec_fields(w=Field) prec(DTYPE(prec))
                  w_p(1:lda) = 0._MK
              end foreach
          END SELECT

          end_subroutine()
      END SUBROUTINE

      SUBROUTINE DTYPE(part_comp_global_index)(Pc,info)
          !!! Compute a global index for particles
          !!! (Uses MPI communications)
          !-------------------------------------------------------------------------
          !  Modules
          !-------------------------------------------------------------------------
          USE ppm_module_mpi
          IMPLICIT NONE

          DEFINE_MK()
          !-------------------------------------------------------------------------
          !  Include
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          ! Arguments
          !-------------------------------------------------------------------------
          CLASS(DTYPE(ppm_t_particles)) :: Pc

          INTEGER,       INTENT(   OUT) :: info
          !!! return status. On success, 0

          !-------------------------------------------------------------------------
          ! local variables
          !-------------------------------------------------------------------------
          INTEGER                        :: i,offset
          INTEGER, DIMENSION(:), POINTER :: wp => NULL()
#ifdef __MPI3
          INTEGER                        :: request
#endif
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          start_subroutine("part_comp_global_index")

          offset=0
#ifdef __MPI
#ifdef __MPI3
          CALL MPI_Iexscan(Pc%Npart,offset,1,MPI_INTEGER,MPI_SUM,ppm_comm,request,info)
          or_fail_MPI("MPI_Iexscan")
#else
          CALL MPI_Exscan(Pc%Npart,offset,1,MPI_INTEGER,MPI_SUM,ppm_comm,info)
          or_fail_MPI("MPI_Exscan")
#endif
#endif

          IF (.NOT.Pc%flags(ppm_part_global_index)) THEN
             CALL Pc%create_prop(info,part_prop=Pc%gi, &
             &    dtype=ppm_type_int,name="GlobalIndex")
             Pc%flags(ppm_part_global_index)=.TRUE.
          END IF

          CALL Pc%get(Pc%gi,wp,info)
#ifdef __MPI3
          CALL MPI_Wait(request,MPI_STATUS_IGNORE,info)
          or_fail_MPI("MPI_Wait")
#endif
          FORALL (i=1:Pc%Npart) wp(i) = offset + i !- 1 !uncomment if index from 0

          CALL Pc%set(Pc%gi,wp,info)

          end_subroutine()
      END SUBROUTINE DTYPE(part_comp_global_index)

      !!----------------------------------------------------------------
      !! Procedures for Particle Sets DS
      !!----------------------------------------------------------------
      SUBROUTINE DTYPE(get_xp)(this,xp,info,with_ghosts)
          IMPLICIT NONE

          DEFINE_MK()
          CLASS(DTYPE(ppm_t_particles))           :: this

          REAL(MK), DIMENSION(:,:), POINTER       :: xp

          INTEGER,                  INTENT(  OUT) :: info

          LOGICAL,        OPTIONAL, INTENT(IN   ) :: with_ghosts
          !-------------------------------------------------------------------------
          ! local variables
          !-------------------------------------------------------------------------
          start_subroutine("get_xp")

          check_associated(<#this%xp#>)

          IF (PRESENT(with_ghosts)) THEN
             IF (with_ghosts) THEN
                IF (this%flags(ppm_part_ghosts)) THEN
                   xp => this%xp(1:ppm_dim,1:this%Mpart)
                ELSE
                   stdout("WARNING: tried to get xp with ghosts when ghosts are not up-to-date")

                   xp => NULL()
                ENDIF
                RETURN
             ENDIF
          ENDIF

          xp => this%xp(1:ppm_dim,1:this%Npart)

          end_subroutine()
      END SUBROUTINE DTYPE(get_xp)

      SUBROUTINE DTYPE(set_xp)(this,xp,info,read_only,ghosts_ok)
          IMPLICIT NONE

          DEFINE_MK()
          CLASS(DTYPE(ppm_t_particles))                    :: this

          REAL(MK), DIMENSION(:,:), POINTER                :: xp

          INTEGER,                           INTENT(  OUT) :: info

          LOGICAL,                 OPTIONAL, INTENT(IN   ) :: read_only
          LOGICAL,                 OPTIONAL, INTENT(IN   ) :: ghosts_ok

          !-------------------------------------------------------------------------
          ! local variables
          !-------------------------------------------------------------------------
          CLASS(DTYPE(ppm_t_part_prop)_), POINTER :: prop
          CLASS(DTYPE(ppm_t_neighlist)_), POINTER :: Nl
          CLASS(ppm_t_operator_discr_),   POINTER :: op

          INTEGER :: i

          start_subroutine("set_xp")

          IF (PRESENT(ghosts_ok)) THEN
             IF (ghosts_ok) THEN
                xp => NULL()
                RETURN
             ENDIF
          ENDIF

          IF (PRESENT(read_only)) THEN
             IF (read_only) THEN
                xp => NULL()
                RETURN
             ENDIF
          ENDIF

          this%flags(ppm_part_areinside) = .FALSE.
          this%flags(ppm_part_partial)   = .FALSE.
          this%flags(ppm_part_ghosts)    = .FALSE.
          this%flags(ppm_part_cartesian) = .FALSE.

          prop => this%props%begin()
          DO WHILE (ASSOCIATED(prop))
             prop%flags(ppm_ppt_ghosts)  = .FALSE.
             prop%flags(ppm_ppt_partial) = .FALSE.

             prop => this%props%next()
          ENDDO

          Nl => this%neighs%begin()
          DO WHILE (ASSOCIATED(Nl))
             Nl%uptodate = .FALSE.

             Nl => this%neighs%next()
          ENDDO

          IF (ASSOCIATED(this%ops)) THEN
             op => this%ops%begin()
             DO WHILE (ASSOCIATED(op))
                op%flags(ppm_ops_iscomputed) = .FALSE.

                op => this%ops%next()
             ENDDO
          ENDIF

          xp => NULL()

          end_subroutine()
      END SUBROUTINE DTYPE(set_xp)

      SUBROUTINE DTYPE(part_move)(Pc,disp,info)
          !!!  Move all particles according to some displacement field
          !!!  The size of disp must match the size of xp
          !-------------------------------------------------------------------------
          !  Modules
          !-------------------------------------------------------------------------
          IMPLICIT NONE

          DEFINE_MK()
          !-------------------------------------------------------------------------
          ! Arguments
          !-------------------------------------------------------------------------
          CLASS(DTYPE(ppm_t_particles))           :: Pc
          !!! Data structure containing the particles

          REAL(MK), DIMENSION(:,:), TARGET       :: disp
          !!! Data structure containing the particles

          INTEGER,                  INTENT(  OUT) :: info
          !!! Return status, on success 0.
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          CLASS(DTYPE(ppm_t_part_prop)_), POINTER :: prop

          CLASS(ppm_t_operator_discr_),   POINTER :: op

          REAL(MK), DIMENSION(:,:), POINTER :: xp => NULL()

          INTEGER :: ip

          start_subroutine("particle_move")

          !-----------------------------------------------------------------
          !  checks
          !-----------------------------------------------------------------
          CALL Pc%get_xp(xp,info)
          or_fail("Particle positions cannot be accessed")

          !FORALL (ip=1:Pc%Npart) xp(1:ppm_dim,ip) = xp(1:ppm_dim,ip) + disp(1:ppm_dim,ip)
          DO ip=1,Pc%Npart
             xp(1:ppm_dim,ip) = xp(1:ppm_dim,ip) + disp(1:ppm_dim,ip)
          ENDDO

          CALL Pc%set_xp(xp,info)
          or_fail("set_xp")

          !-----------------------------------------------------------------
          !  update states
          !-----------------------------------------------------------------
          Pc%flags(ppm_part_ghosts) = .FALSE.

          prop => Pc%props%begin()
          DO WHILE (ASSOCIATED(prop))
             prop%flags(ppm_ppt_ghosts) = .FALSE.
             prop => Pc%props%next()
          ENDDO

          IF (ASSOCIATED(Pc%ops)) THEN
             op => Pc%ops%begin()
             DO WHILE (ASSOCIATED(op))
                op%flags(ppm_ops_iscomputed) = .FALSE.
                op => Pc%ops%next()
             ENDDO
          ENDIF

          Pc%flags(ppm_part_partial) = .FALSE.
          Pc%flags(ppm_part_cartesian) = .FALSE.

          end_subroutine()
      END SUBROUTINE DTYPE(part_move)

      SUBROUTINE DTYPE(part_apply_bc)(Pc,info)
          !!!  Apply boundary conditions for particles positions
          !!!  Assumptions:
          !!! * Particles positions need to have been mapped onto the topology
          !-------------------------------------------------------------------------
          !  Modules
          !-------------------------------------------------------------------------
          IMPLICIT NONE

          DEFINE_MK()
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(DTYPE(ppm_t_particles)) :: Pc
          !!! Data structure containing the particles

          INTEGER,        INTENT(  OUT) :: info
          !!! Return status, on success 0.
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          CLASS(DTYPE(ppm_t_part_prop)_), POINTER :: prop

          TYPE(ppm_t_topo), POINTER :: topo
          !!! pointer to topology

          REAL(MK), DIMENSION(:,:), POINTER :: xp
          !!! pointer to positions
          REAL(MK), DIMENSION(:), POINTER :: min_phys
          REAL(MK), DIMENSION(:), POINTER :: max_phys
          !!! computational domain corners
          REAL(MK), DIMENSION(ppm_dim)    :: len_phys
          !!! length of the computational domain
          REAL(MK)                        :: almostone

          INTEGER, DIMENSION(:), POINTER :: list_del_parts
          INTEGER                        :: di,ip
          INTEGER                        :: topoid
          INTEGER                        :: Npart,del_part

          start_subroutine("particles_apply_bc")

          !-----------------------------------------------------------------
          !  checks
          !-----------------------------------------------------------------
          IF (.NOT.ASSOCIATED(Pc%xp)) THEN
             fail('Particles structure had not been defined. Call allocate first')
          ENDIF

          topoid   =Pc%active_topoid
          topo     => ppm_topo(topoid)%t
          xp       => Pc%xp
          Npart    =Pc%Npart
          almostone=NEAREST(1._MK,-1._MK)

          !-----------------------------------------------------------------
          !  Move particles if needed
          !-----------------------------------------------------------------
#if   __KIND == __SINGLE_PRECISION
          min_phys => topo%min_physs
          max_phys => topo%max_physs
#elif __KIND == __DOUBLE_PRECISION
          min_phys => topo%min_physd
          max_phys => topo%max_physd
#endif
          len_phys(1:ppm_dim)=max_phys(1:ppm_dim)-min_phys(1:ppm_dim)

          del_part = 0
          DO di=1,ppm_dim
             SELECT CASE (topo%bcdef(di))
             CASE (ppm_param_bcdef_periodic)
                DO ip=1,Npart
                   IF (xp(di,ip) .EQ. max_phys(di)) &
                   &   xp(di,ip) = xp(di,ip) - len_phys(di)*almostone
                   IF (xp(di,ip) .GT. max_phys(di)) &
                   &   xp(di,ip) = xp(di,ip) - len_phys(di)
                   IF (xp(di,ip) .LT. min_phys(di)) &
                   &   xp(di,ip) = xp(di,ip) + len_phys(di)
                ENDDO

             CASE (ppm_param_bcdef_freespace)
                !delete particles that have crossed the boundary
                DO ip=Npart,1,-1
                   IF (xp(di,ip).GT.max_phys(di).OR.xp(di,ip).LT.min_phys(di)) THEN
                      del_part = del_part+1
                   ENDIF
                ENDDO

             CASE DEFAULT
                fail("this type of BC is not implemented/tested in this version")

             END SELECT
          ENDDO !di=1,ppm_dim

          IF (del_part.GT.0) THEN
             NULLIFY(list_del_parts)
             ldc(1) = del_part
             CALL ppm_alloc(list_del_parts,ldc,ppm_param_alloc_fit,info)
             or_fail_alloc("list_del_parts")

             del_part = 0
             DO di=1,ppm_dim
                IF (topo%bcdef(di).EQ.ppm_param_bcdef_freespace) THEN
                   DO ip=Npart,1,-1
                      !TOCHECK
                      !Yaser I changed this to GT instead of GE as for free space boundary conditions
                      !you might have some particles on the max_phys border of the domain
                      IF (xp(di,ip).GT.max_phys(di).OR.xp(di,ip).LT.min_phys(di)) THEN
                         del_part = del_part+1
                         list_del_parts(del_part)=ip
                      ENDIF
                   ENDDO
                ENDIF
             ENDDO

             CALL Pc%del_parts(list_del_parts,del_part,info)
             or_fail("Pc%del_parts: could not delete particles")

             CALL ppm_alloc(list_del_parts,ldc,ppm_param_dealloc,info)
             or_fail_dealloc("list_del_parts")
          ENDIF

          ! Update states
          Pc%Npart = Npart
          ! Particles are no longer on the right processors
          Pc%flags(ppm_part_partial) = .FALSE.
          ! But they are now all inside the computational domain
          Pc%flags(ppm_part_areinside) = .TRUE.
          ! Dangereous to use the ghosts
          Pc%flags(ppm_part_ghosts) = .FALSE.
          ! ghosts values for properties are also dangerous to use
          prop => Pc%props%begin()
          DO WHILE (ASSOCIATED(prop))
             prop%flags(ppm_ppt_ghosts) = .FALSE.
             prop => Pc%props%next()
          ENDDO

          !-----------------------------------------------------------------------
          ! Finalize
          !-----------------------------------------------------------------------

          end_subroutine()
      END SUBROUTINE DTYPE(part_apply_bc)

      SUBROUTINE DTYPE(part_print_info)(Pc,info,level,fileunit)
          !-----------------------------------------------------------------------
          ! Print out summary information about this Particle set
          ! (list of properties, operators, etc...)
          !-----------------------------------------------------------------------
          IMPLICIT NONE

          DEFINE_MK()
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(DTYPE(ppm_t_particles))    :: Pc
          !!! Data structure containing the particles
          INTEGER,           INTENT(  OUT) :: info
          !!! Return status, on success 0.
          !-------------------------------------------------------------------------
          !  Optional arguments
          !-------------------------------------------------------------------------
          INTEGER, OPTIONAL, INTENT(IN   ) :: level
          !!! indentation level at which to printout the info. Default = 0
          INTEGER, OPTIONAL, INTENT(IN   ) :: fileunit
          !!! Already open file unit for printout. Default = stdout
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          CLASS(DTYPE(ppm_t_part_prop)_), POINTER :: prop

          INTEGER :: lev,fileu

          CHARACTER(LEN = ppm_char) :: myformat

          start_subroutine("part_print_info")

          IF (PRESENT(fileunit)) THEN
             fileu = fileunit
          ELSE
             fileu = 6
          ENDIF
          IF (PRESENT(level)) THEN
             lev = MAX(level,1)
          ELSE
             lev = 1
          ENDIF

          WRITE(myformat,'(A,I0,A)') '(',4*lev,'X,A,A,2X,2(A,I0),A)'

          WRITE(fileu,myformat) 'Particle set: ',TRIM(color_print(Pc%name,31)),&
          & '(N = ',Pc%Npart,' M = ',Pc%Mpart,')'

          lev = lev + 1

          WRITE(myformat,'(A,I0,A,I0,A)') '(',4*lev,'X,A,',&
          & ppm_param_length_partflags,'L)'
          WRITE(fileu,myformat) 'flags: ',Pc%flags

          prop => Pc%props%begin()
          DO WHILE (ASSOCIATED(prop))
             CALL prop%print_info(info,lev,fileunit,Pc%props%iter_id)
             prop => Pc%props%next()
          ENDDO

          end_subroutine()
      END SUBROUTINE DTYPE(part_print_info)

      SUBROUTINE DTYPE(part_set_cutoff)(Pc,cutoff,info,Nlist)
          !!! Set a cutoff radius for a Particle set and update the
          !!! the ghostlayer sizes.
          !!! The cutoff radius concerns the default neighbor list, unless
          !!! specified otherwise.
          !!! If the cutoff is increased from its previous value, the neighbour
          !!! list is flagged as "not up-to-date" and will have to be recomputed
          !!! before it is used again.
          !!! If the ghostlayer sizes are increased from their previous values,
          !!! the ghosts are flagged as "not up-to-date" and will have to be
          !!! recomputed before they are used again.
          !-------------------------------------------------------------------------
          !  Modules
          !-------------------------------------------------------------------------
          IMPLICIT NONE

          DEFINE_MK()
          !-------------------------------------------------------------------------
          ! Arguments
          !-------------------------------------------------------------------------
          CLASS(DTYPE(ppm_t_particles))                           :: Pc
          REAL(MK),                                 INTENT(IN   ) :: cutoff
          !!! cutoff radius
          INTEGER,                                  INTENT(  OUT) :: info
          !!! return status. On success, 0
          CLASS(DTYPE(ppm_t_neighlist)_), OPTIONAL, INTENT(INOUT) :: NList
          !!! Neighbor list for which this cutoff radius
          !!! applies. By default, this is the "standard" Verlet list, with neighbours
          !!! sought within the particle set itself.

          !-------------------------------------------------------------------------
          ! local variables
          !-------------------------------------------------------------------------
          CLASS(DTYPE(ppm_t_neighlist)_), POINTER :: nl

          start_subroutine("part_set_cutoff")

          !-------------------------------------------------------------------------
          !  Set new cutoff
          !-------------------------------------------------------------------------
          IF (PRESENT(NList)) THEN
             check_true(<#Pc%neighs%has(NList)#>,&
             & "Neighbour list does not concern this particle set")
             IF (cutoff .LT. NList%cutoff) NList%uptodate = .FALSE.
             NList%cutoff = cutoff
          ELSE
             nl => Pc%get_neighlist()
             check_associated(nl,"Compute neighbour lists first")
             IF (cutoff .LT. nl%cutoff) nl%uptodate = .FALSE.
             nl%cutoff = cutoff
          ENDIF

          ! Compute ghostlayer sizes
          IF (cutoff.GT.Pc%ghostlayer) THEN
             !If the new cutoff is larger than the current ghostsize
             ! then the new ghostsize is the new cutoff
             Pc%ghostlayer = cutoff
             ! update states
             Pc%flags(ppm_part_ghosts) = .FALSE.
          ELSE IF (cutoff .LT. Pc%ghostlayer) THEN
             !Else, we find the new maximum cutoff radius amongst
             !all existing neighbor lists on this Particle set
             Pc%ghostlayer = 0._MK

             nl => Pc%neighs%begin()
             DO WHILE (ASSOCIATED(nl))
                IF (nl%cutoff .GT. Pc%ghostlayer) THEN
                   Pc%ghostlayer = nl%cutoff
                ENDIF
                nl => Pc%neighs%next()
             ENDDO
             !no need to update states: ghosts are still ok.
          ENDIF

          end_subroutine()
      END SUBROUTINE DTYPE(part_set_cutoff)

minclude ppm_create_collection_procedures(DTYPE(neighlist),DTYPE(neighlist)_)

      SUBROUTINE DTYPE(part_neigh_create)(this,Part_src,info,&
      &          name,skin,symmetry,cutoff,Nlist)
          !!! Create a data structure to store a neighbour list
          !-------------------------------------------------------------------------
          !  Modules
          !-------------------------------------------------------------------------
          IMPLICIT NONE

          DEFINE_MK()

          CLASS(DTYPE(ppm_t_particles))                           :: this
          CLASS(DTYPE(ppm_t_particles)_),           TARGET        :: Part_src
          !!! Particle set to which the neighbours belong (can be the same as this)
          INTEGER,                                  INTENT(  OUT) :: info
          CHARACTER(LEN=*),               OPTIONAL, INTENT(IN   ) :: name
          !!! name of this neighbour list
          REAL(MK),                       OPTIONAL, INTENT(IN   ) :: skin
          REAL(MK),                       OPTIONAL, INTENT(IN   ) :: cutoff
          LOGICAL,                        OPTIONAL, INTENT(IN   ) :: symmetry
          CLASS(DTYPE(ppm_t_neighlist)_), OPTIONAL, POINTER       :: Nlist
          !!! returns a pointer to the newly created verlet list
          !-------------------------------------------------------------------------
          ! local variables
          !-------------------------------------------------------------------------
          CLASS(DTYPE(ppm_t_neighlist)_), POINTER :: Nl

!           REAL(MK), DIMENSION(:), POINTER :: rcp => NULL()

          INTEGER :: vec_size,i

          start_subroutine("particle_neigh_create")

          ! Create the neighbour list
          ALLOCATE(DTYPE(ppm_t_neighlist)::Nl,STAT=info)
          or_fail_alloc("Nl")

          IF (PRESENT(name)) THEN
             Nl%name = name
          ELSE
             WRITE(Nl%name,*) 'Nl',TRIM(ADJUSTL(this%name)),'_',&
             & TRIM(ADJUSTL(Part_src%name))
          ENDIF

          check_associated(<#Part_src%xp#>,"Invalid particle set Part_src")

          Nl%Part => Part_src

          Nl%cutoff=MERGE(cutoff,this%ghostlayer,PRESENT(cutoff))
          Nl%skin  =MERGE(skin,0.0_MK,PRESENT(skin))

          IF (Part_src%isymm.EQ.1) THEN
             !Yaser
             !You can not have symmetric particles with NonSymmetric
             !neighborlist, (the particle symmetry is the dominant factor
             !in neighborlist symmetry)
             Nl%isymm = 1
          ELSE
             !You can have a NonSymmetric particles with symmetric neighborlist
             IF (PRESENT(symmetry)) THEN
                Nl%isymm = MERGE(1,0,symmetry)
             ELSE
                Nl%isymm = 0
             ENDIF
          ENDIF

          Nl%uptodate = .FALSE.
          Nl%nneighmin = 0
          Nl%nneighmax = 0

          !returning a pointer to the neighbour list, before it is pushed
          ! into the collection.
          IF (PRESENT(Nlist)) Nlist => Nl

          CALL this%neighs%push(Nl,info)
          or_fail("pushing new neighbour list into collection failed")

          end_subroutine()
      END SUBROUTINE DTYPE(part_neigh_create)

      SUBROUTINE DTYPE(part_neigh_destroy)(this,Nlist,info)
          !!! Destroy a property from an existing particle set
          !-------------------------------------------------------------------------
          !  Modules
          !-------------------------------------------------------------------------
          IMPLICIT NONE

          CLASS(DTYPE(ppm_t_particles))                 :: this
          CLASS(DTYPE(ppm_t_neighlist)_), POINTER       :: Nlist

          INTEGER,                        INTENT(  OUT) :: info
          !-------------------------------------------------------------------------
          ! local variables
          !-------------------------------------------------------------------------
          start_subroutine("part_neigh_destroy")

          CALL this%neighs%remove(info,Nlist)
          or_fail("could not remove Nlist from its collection")

          end_subroutine()
      END SUBROUTINE DTYPE(part_neigh_destroy)

      SUBROUTINE DTYPE(part_neighlist)(this,info,P_xset,name, &
      &          skin,symmetry,cutoff,lstore,incl_ghosts,knn)
          !!!  Neighbor lists for particles
          !!!  Compute the Verlet lists for the target particles, using neighbours
          !!!  from the particle set P_set (the default is that P_xset is the
          !!!  same set as the target particles)
          !!!-----------------------------------------------------------------
          !!!  Assumptions:
          !!! * Particles positions need to have been mapped onto the topology
          !!! * Ghost positions have been computed
          !-------------------------------------------------------------------------
          !  Modules
          !-------------------------------------------------------------------------
          USE ppm_module_mpi
          USE ppm_module_time, ONLY : ppm_time
          USE ppm_module_neighlist
#ifdef __WITH_CNL
          USE ppm_module_cnl
#endif
          USE ppm_module_inl_vlist
          USE ppm_module_inl_xset_vlist
          USE ppm_module_inl_k_vlist
          USE ppm_module_kdtree
          IMPLICIT NONE

          DEFINE_MK()
          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(DTYPE(ppm_t_particles)),            TARGET        :: this
          !!! Data structure containing the particles
          INTEGER,                                  INTENT(  OUT) :: info
          !!! Return status, on success 0.
          !-------------------------------------------------------------------------
          !  Optional arguments
          !-------------------------------------------------------------------------
          CLASS(DTYPE(ppm_t_particles)_), OPTIONAL, TARGET        :: P_xset
          !!! Particle set from which the neighbours are sought
          CHARACTER(LEN=*),               OPTIONAL, INTENT(IN   ) :: name
          !!! name of this neighbour list
          REAL(MK),                       OPTIONAL, INTENT(IN   ) :: skin
          !!! skin
          LOGICAL,                        OPTIONAL, INTENT(IN   ) :: symmetry
          !!! if using symmetry
          REAL(MK),                       OPTIONAL, INTENT(IN   ) :: cutoff
          !!! cutoff radius
          LOGICAL,                        OPTIONAL, INTENT(IN   ) :: lstore
          !!! store verlet lists
          LOGICAL,                        OPTIONAL, INTENT(IN   ) :: incl_ghosts
          !!! if true, then verlet lists are computed for all particles, incl. ghosts.
          !!! Default is false.
          INTEGER,                        OPTIONAL, INTENT(IN   ) :: knn
          !!! if present, neighbour lists are constructed such that each particle
          !!! has at least knn neighbours.
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          CLASS(ppm_t_operator_discr_),   POINTER :: op
          CLASS(DTYPE(ppm_t_particles)_), POINTER :: Part_src
          CLASS(DTYPE(ppm_t_neighlist)_), POINTER :: Nlist

          TYPE(ppm_t_topo), POINTER :: topo

          TYPE(DTYPE(kdtree)), POINTER :: tree

          TYPE(DTYPE(kdtree_result)), DIMENSION(:), ALLOCATABLE, TARGET :: results

          REAL(MK)                      :: lskin
          REAL(MK),DIMENSION(2*ppm_dim) :: ghostlayer
          REAL(ppm_kind_double)         :: t1,t2
          REAL(MK)                      :: cutoff_,skin_

          INTEGER :: topoid
          INTEGER :: nneighmin,nneighmax
          INTEGER :: op_id,np_target,i
          INTEGER :: ip,ineigh
          !!! index variable


          LOGICAL :: lknn,lsymm
          !!! uses a neighbour-finding algorithm that finds enough neighbours
          LOGICAL :: xset_neighlists

          LOGICAL :: symmetry_,lstore_

          CHARACTER(LEN=ppm_char) :: name_

          start_subroutine("part_comp_neighlist")

          !-----------------------------------------------------------------
          !  Checks
          !-----------------------------------------------------------------
          check_associated(<#this%xp#>,&
          & "Particles structure had not been defined. Call allocate first")

          check_true(<#this%flags(ppm_part_partial)#>,&
          & "Particles not mapped. Do a partial/global mapping")

          xset_neighlists = .FALSE.
          IF (PRESENT(P_xset)) THEN
             Part_src => P_xset

             check_associated(<#Part_src%xp#>,&
             & "Cross-Set particles have not been defined. Call allocate first")

             check_true(<#Part_src%flags(ppm_part_partial)#>,&
             & "Particles not mapped. Do a partial/global mapping")

             IF (.NOT.ASSOCIATED(Part_src,this)) THEN
                xset_neighlists = .TRUE.
             ENDIF
          ELSE
             Part_src => this
          ENDIF

          check_true(<#Part_src%flags(ppm_part_ghosts)#>,&
          & "Ghosts have not been updated. They are needed for neighlists")

          check_associated(<#this%neighs#>)

          !check whether the neighbour list already exists
          IF (this%has_neighlist(Part_src)) THEN
             Nlist => this%get_neighlist(Part_src)
             IF (PRESENT(skin).OR.PRESENT(symmetry).OR.PRESENT(cutoff)) THEN
                stdout("the optional arguments skin,",&
                & "symmetry or cutoff will not be used",&
                & " because the neighbour list already exists. We",&
                & " should perhaps change the API?  ")
                fail("Need to destroy/re-create this neighbour list first")
             ENDIF
          ELSE
             symmetry_=MERGE(symmetry,.FALSE.,PRESENT(symmetry))
             cutoff_  =MERGE(cutoff,this%ghostlayer,PRESENT(cutoff))
             skin_    =MERGE(skin,0.0_MK,PRESENT(skin))

             IF (PRESENT(name)) THEN
                name_ = name
             ELSE
                WRITE(name_,*) 'Nl',TRIM(ADJUSTL(this%name)),'_', &
                & TRIM(ADJUSTL(Part_src%name))
             ENDIF

             NULLIFY(Nlist)
             CALL this%create_neighlist(Part_src,info,name=TRIM(name_), &
             &    skin=skin_,symmetry=symmetry_,cutoff=cutoff_,Nlist=Nlist)
             or_fail("failed to create neighbour list")
          ENDIF

          check_associated(Nlist)

          !check that we have a cutoff radius
          check_true(<#Nlist%cutoff.GT.0#>,&
          & "cutoff is negative or zero - do we really want neighbour lists?")

          lsymm =Nlist%isymm.EQ.1
          lknn  =PRESENT(knn)
          lskin =Nlist%skin
          topoid=this%active_topoid

          do_something: IF (Nlist%uptodate .OR. this%Npart.EQ.0) THEN
             !neighbor lists are already up-to-date, or no particles on this proc
             !nothing to do
             IF (Nlist%uptodate) THEN
                fail('neighlists are already up-to-date, NOTHING to do', &
                & 999,exit_point=no,ppm_error=ppm_error_notice)
                info = 0
             ELSE
                Nlist%nneighmin = 0
                Nlist%nneighmax = 0
             ENDIF
          ELSE do_something
             !hack to build (potentially incomplete) neighbour lists even
             !for ghost particles
             np_target = this%Npart
             IF (PRESENT(incl_ghosts)) THEN
                IF (incl_ghosts) THEN
                   np_target = this%Mpart
                   topo => ppm_topo(topoid)%t
                   IF (MK.EQ.ppm_kind_single) THEN
                      topo%min_subs(:,:) = topo%min_subs(:,:) - topo%ghostsizes
                      topo%max_subs(:,:) = topo%max_subs(:,:) + topo%ghostsizes
                   ELSE IF (MK.EQ.ppm_kind_double) THEN
                      topo%min_subd(:,:) = topo%min_subd(:,:) - topo%ghostsized
                      topo%max_subd(:,:) = topo%max_subd(:,:) + topo%ghostsized
                   ENDIF
                ENDIF
             ENDIF

             IF (lknn) THEN
                this%stats%nb_kdtree = this%stats%nb_kdtree+1

                CALL ppm_time(t1,info)

                AllOCATE(tree,STAT=info)
                or_fail_alloc("tree")

                CALL tree%create(Part_src%xp(1:ppm_dim,1:Part_src%Mpart), &
                &    info,sort=.TRUE.,rearrange=.TRUE.)
                or_fail("tree%create")

                ALLOCATE(results(knn+1),STAT=info)
                or_fail_alloc("results")

                ldc(1) = knn
                ldc(2) = np_target
                CALL ppm_alloc(Nlist%vlist,ldc,ppm_param_alloc_grow,info)
                or_fail_alloc("Nlist%vlist")

                DO ip=1,np_target
                   CALL kdtree_n_nearest(tree,this%xp(1:ppm_dim,ip),&
                   &    knn+1,results,info)
                   or_fail("kdtree_n_nearest")

                   ! If the tree is not sorted we need to remove the
                   ! particle ip from the list of neighbors
                   Nlist%vlist(1:knn,ip)=results(2:knn+1)%idx
                ENDDO

                ldc(1) = np_target
                CALL ppm_alloc(Nlist%nvlist,ldc,ppm_param_alloc_grow,info)
                or_fail_alloc("Nlist%nvlist")

                Nlist%nvlist=knn

                CALL tree%destroy(info)
                or_fail("tree%destroy")

                DEALLOCATE(tree,results,STAT=info)
                or_fail_dealloc("tree & results")
                NULLIFY(tree)

                CALL ppm_time(t2,info)

                this%stats%t_kdtree = this%stats%t_kdtree+(t2-t1)

             ELSE
                lstore_=MERGE(lstore,.TRUE.,PRESENT(lstore))

                IF (xset_neighlists) THEN
                   this%stats%nb_xset_nl = this%stats%nb_xset_nl + 1

                   CALL ppm_time(t1,info)

                   CALL ppm_inl_xset_vlist(topoid,this%xp,this%Npart,         &
                   &    this%Mpart,Part_src%xp,Part_src%Npart,Part_src%Mpart, &
                   &    Nlist%cutoff,lskin,ghostlayer,info,Nlist%vlist,       &
                   &    Nlist%nvlist,lstore_)
                   or_fail("ppm_inl_xset_vlist failed")

                   CALL ppm_time(t2,info)

                   this%stats%t_xset_nl = this%stats%t_xset_nl + (t2 - t1)
                ELSE
                   this%stats%nb_nl = this%stats%nb_nl+1

                   CALL ppm_time(t1,info)

                   CALL ppm_neighlist_vlist(topoid,this%xp,this%Mpart,     &
                   &    Nlist%cutoff,lskin,lsymm,Nlist%vlist,Nlist%nvlist, &
                   &    info,lstore=lstore_)
                   or_fail("ppm_neighlist_vlist failed")

                   CALL ppm_time(t2,info)

                   this%stats%t_nl = this%stats%t_nl + (t2 - t1)
                ENDIF ! XSET
             ENDIF

             !restore subdomain sizes (revert hack)
             IF (PRESENT(incl_ghosts)) THEN
                IF (incl_ghosts) THEN
                   IF (MK.EQ.ppm_kind_single) THEN
                      topo%min_subs(:,:) = topo%min_subs(:,:) + topo%ghostsizes
                      topo%max_subs(:,:) = topo%max_subs(:,:) - topo%ghostsizes
                   ELSE IF (MK.EQ.ppm_kind_double) THEN
                      topo%min_subd(:,:) = topo%min_subd(:,:) + topo%ghostsized
                      topo%max_subd(:,:) = topo%max_subd(:,:) - topo%ghostsized
                   ENDIF
                   topo => NULL()
                ENDIF
             ENDIF

             !-----------------------------------------------------------------------
             !Update state
             !-----------------------------------------------------------------------
             Nlist%uptodate = .TRUE.

             Nlist%nneighmin = MINVAL(Nlist%nvlist(1:this%Npart))
             Nlist%nneighmax = MAXVAL(Nlist%nvlist(1:np_target))

             ! DC operators that do not use a xset neighbour list, if they exist,
             ! are no longer valid (they depend on the neighbour lists)
             IF (ASSOCIATED(this%ops)) THEN
                op => this%ops%begin()
                DO WHILE (ASSOCIATED(op))
                   IF (.NOT.op%flags(ppm_ops_interp)) THEN
                      op%flags(ppm_ops_iscomputed) = .FALSE.
                   ENDIF
                   op => this%ops%next()
                ENDDO
             ENDIF

             ! We want to distinguish between "self" neighbour lists
             ! and cross-set ones.
             IF (ASSOCIATED(Nlist%Part,this)) THEN
                this%flags(ppm_part_neighlists) = .TRUE.
             ENDIF

             Nlist => NULL()

          ENDIF do_something

          end_subroutine()
      END SUBROUTINE DTYPE(part_neighlist)

      FUNCTION DTYPE(has_ghosts)(this,Field) RESULT(res)
          !!! Returns true if the discretization of the field on this
          !!! particle set has its ghosts up-to-date.
          !!! If the Field argument is not present, then the function
          !!! returns whether the particles themselves have their
          !!! ghosts up-to-date.
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(DTYPE(ppm_t_particles)) :: this
          CLASS(ppm_t_field_), OPTIONAL :: Field
          LOGICAL                       :: res
          !-------------------------------------------------------------------------
          ! local variables
          !-------------------------------------------------------------------------
          CLASS(ppm_t_discr_data), POINTER :: prop

          start_function("has_ghosts")

          IF (PRESENT(Field)) THEN
             NULLIFY(prop)
             CALL this%get_discr(Field,prop,info)
             or_fail("this field is not discretized on this particle set")

             res = prop%flags(ppm_ppt_ghosts)
          ELSE
             res = this%flags(ppm_part_ghosts)
          ENDIF

          end_function()

      END FUNCTION DTYPE(has_ghosts)

      SUBROUTINE DTYPE(part_map_create)(Pc,id,source_topoid,target_topoid,info)
          !!! Adds a property to an existing particle set
          !-------------------------------------------------------------------------
          !  Modules
          !-------------------------------------------------------------------------
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(DTYPE(ppm_t_particles)) :: Pc

          INTEGER,        INTENT(  OUT) :: id
          INTEGER,        INTENT(IN   ) :: source_topoid
          INTEGER,        INTENT(IN   ) :: target_topoid
          INTEGER,        INTENT(  OUT) :: info

          !-------------------------------------------------------------------------
          ! local variables
          !-------------------------------------------------------------------------
          TYPE(DTYPE(ppm_t_ptr_part_mapping)), DIMENSION(:), POINTER:: vec_tmp

          CLASS(DTYPE(ppm_t_part_mapping)_), POINTER:: map

          INTEGER :: vec_size,npart,i

          LOGICAL, DIMENSION(ppm_param_length_pptflags):: flags

          start_subroutine("part_map_create")

          !Generate a new id (we should use templating here...)
          ASSOCIATE (cont => Pc%maps )
             id = 0
             IF (cont%nb.LT.cont%vec_size) THEN
                !there is at least one empty slot in the array
                ! of mapping pointers
                id = id + 1
                DO WHILE (ASSOCIATED(cont%vec(id)%t))
                   id = id + 1
                ENDDO
             ELSE
                IF (.NOT. ASSOCIATED(cont%vec)) THEN
                   !need to allocate the array of mapping pointers
                   vec_size=16
                   ALLOCATE(cont%vec(vec_size),STAT=info)
                   or_fail_alloc("cont%vec")

                   id = 1
                ELSE
                   !need to resize the array of mapping pointers
                   vec_size=MAX(2*cont%vec_size,16)

                   ALLOCATE(vec_tmp(vec_size),STAT=info)
                   or_fail_alloc("vec_tmp")

                   DO i=1,cont%vec_size
                      vec_tmp(i)%t => cont%vec(i)%t
                   ENDDO

                   DEALLOCATE(cont%vec)
                   cont%vec => vec_tmp
                ENDIF

                cont%vec_size = vec_size

                id = cont%nb + 1
             ENDIF
             cont%nb = cont%nb + 1

             IF (id .GT. cont%max_id) cont%max_id = id
             IF (id .LT. cont%min_id) cont%min_id = id

          END ASSOCIATE

          IF (.NOT. ASSOCIATED(Pc%maps%vec(id)%t)) THEN
             ALLOCATE(DTYPE(ppm_t_part_mapping)::Pc%maps%vec(id)%t,STAT=info)
             or_fail_alloc("Pc%maps%vec(id)%t")
          ENDIF

          map => Pc%maps%vec(id)%t

          ! Create the mapping
          CALL map%create(source_topoid,target_topoid,info)
          or_fail("map%create")

          end_subroutine()
      END SUBROUTINE DTYPE(part_map_create)

      SUBROUTINE DTYPE(part_map_destroy)(Pc,id,info)
          !!! Destroy a mapping from an existing particle set
          !-------------------------------------------------------------------------
          !  Modules
          !-------------------------------------------------------------------------
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(DTYPE(ppm_t_particles)) :: Pc

          INTEGER,        INTENT(INOUT) :: id
          INTEGER,        INTENT(  OUT) :: info
          !-------------------------------------------------------------------------
          ! local variables
          !-------------------------------------------------------------------------
          start_subroutine("part_map_destroy")

          ASSOCIATE (cont => Pc%maps)
             IF (id .LE. 0 .OR. id .GT. cont%vec_size) THEN
                fail("mapping id larger than size of mappings array")
             ENDIF

             CALL cont%vec(id)%t%destroy(info)
             or_fail("Pc%maps%vec(id)%t%destroy")

             DEALLOCATE(cont%vec(id)%t,STAT=info)
             or_fail_dealloc("Failed to deallocate property")
             NULLIFY(cont%vec(id)%t)

             cont%nb = cont%nb - 1
             IF (id .EQ. cont%max_id) THEN
                cont%max_id = cont%max_id - 1
                IF (cont%max_id .GT. 0) THEN
                   DO WHILE(.NOT.ASSOCIATED(cont%vec(cont%max_id)%t))
                      cont%max_id = cont%max_id - 1
                      IF (cont%max_id .EQ. 0) EXIT
                   ENDDO
                ENDIF
             ENDIF
             IF (cont%nb.EQ.0) THEN
                cont%min_id = ppm_big_i
             ELSE IF (id .EQ. cont%min_id) THEN
                cont%min_id = cont%min_id + 1
                IF (cont%min_id .LE. cont%vec_size) THEN
                   DO WHILE(.NOT.ASSOCIATED(cont%vec(cont%min_id)%t))
                      cont%min_id = cont%min_id + 1
                      IF (cont%min_id .GT. cont%vec_size) THEN
                         fail("fatal error in the data structure")
                      ENDIF
                  ENDDO
                ENDIF
             ENDIF
          END ASSOCIATE

          end_subroutine()
      END SUBROUTINE DTYPE(part_map_destroy)

      SUBROUTINE DTYPE(part_prop_push)(Pc,prop_id,info)

          !!! wrapper for ppm_map_part_push
          !-------------------------------------------------------------------------
          !  Modules
          !-------------------------------------------------------------------------
          USE ppm_module_map
          IMPLICIT NONE

          DEFINE_MK()
          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(DTYPE(ppm_t_particles)) :: Pc
          !!! Data structure containing the particles
          INTEGER,        INTENT(IN   ) :: prop_id
          !!! id of the property to be pushed
          INTEGER,        INTENT(  OUT) :: info
          !!! Returns status, 0 upon success.
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          CLASS(DTYPE(ppm_t_part_prop)_), POINTER :: prop

          INTEGER :: lda

          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          start_subroutine("map_part_push")

          !-----------------------------------------------------------------
          !  Call ppm_map_part_push
          !-----------------------------------------------------------------
          prop => Pc%props%vec(prop_id)%t

          lda = prop%lda

          IF (lda.GE.2) THEN
             SELECT CASE (prop%data_type)
             CASE (ppm_type_int)
                CALL ppm_map_part_push(prop%data_2d_i,lda,Pc%Npart,info)

             CASE (ppm_type_longint)
                fail('Type not supported for mappings.')

             CASE (ppm_type_real,ppm_type_real_single)
                CALL ppm_map_part_push(prop%data_2d_r,lda,Pc%Npart,info)

             CASE (ppm_type_comp,ppm_type_comp_single)
                CALL ppm_map_part_push(prop%data_2d_c,lda,Pc%Npart,info)

             CASE (ppm_type_logical)
                CALL ppm_map_part_push(prop%data_2d_l,lda,Pc%Npart,info)

             END SELECT
          ELSE
             SELECT CASE (prop%data_type)
             CASE (ppm_type_int)
                CALL ppm_map_part_push(prop%data_1d_i,Pc%Npart,info)

             CASE (ppm_type_longint)
                fail('Type not supported for mappings.')

             CASE (ppm_type_real,ppm_type_real_single)
                CALL ppm_map_part_push(prop%data_1d_r,Pc%Npart,info)

             CASE (ppm_type_comp,ppm_type_comp_single)
                CALL ppm_map_part_push(prop%data_1d_c,Pc%Npart,info)

             CASE (ppm_type_logical)
                CALL ppm_map_part_push(prop%data_1d_l,Pc%Npart,info)

             END SELECT
          ENDIF
          or_fail("ppm_map_part_push failed!")

          end_subroutine()
      END SUBROUTINE DTYPE(part_prop_push)

      SUBROUTINE DTYPE(part_prop_pop)(Pc,prop_id,Npart_new,info)
          !!! wrapper for ppm_map_part_pop

          !-------------------------------------------------------------------------
          !  Modules
          !-------------------------------------------------------------------------
          USE ppm_module_map
          IMPLICIT NONE

          DEFINE_MK()
          !-------------------------------------------------------------------------
          ! Include
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(DTYPE(ppm_t_particles)) :: Pc
          !!! Data structure containing the particles
          INTEGER,       INTENT(IN   )  :: prop_id
          !!! id of the property to be pooped
          INTEGER,       INTENT(IN   )  :: Npart_new
          !!! number of particles to pop from the buffer
          INTEGER,       INTENT(  OUT)  :: info
          !!! Returns status, 0 upon success.
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          CLASS(DTYPE(ppm_t_part_prop)_), POINTER :: prop

          INTEGER :: lda

          start_subroutine("map_part_pop")

          !-----------------------------------------------------------------
          !  Call ppm_map_part_pop
          !-----------------------------------------------------------------
          prop => Pc%props%vec(prop_id)%t

          lda = prop%lda

          IF (lda.GE.2) THEN
             SELECT CASE (prop%data_type)
             CASE (ppm_type_int)
                CALL ppm_map_part_pop(prop%data_2d_i,lda,Pc%Npart,Npart_new,info)

             CASE (ppm_type_longint)
                fail('Type not supported for mappings.')

             CASE (ppm_type_real,ppm_type_real_single)
                CALL ppm_map_part_pop(prop%data_2d_r,lda,Pc%Npart,Npart_new,info)

             CASE (ppm_type_comp,ppm_type_comp_single)
                CALL ppm_map_part_pop(prop%data_2d_c,lda,Pc%Npart,Npart_new,info)

             CASE (ppm_type_logical)
                CALL ppm_map_part_pop(prop%data_2d_l,lda,Pc%Npart,Npart_new,info)

             END SELECT
          ELSE
             SELECT CASE (prop%data_type)
             CASE (ppm_type_int)
                CALL ppm_map_part_pop(prop%data_1d_i,Pc%Npart,Npart_new,info)

             CASE (ppm_type_longint)
                fail('Type not supported for mappings.')

             CASE (ppm_type_real,ppm_type_real_single)
                CALL ppm_map_part_pop(prop%data_1d_r,Pc%Npart,Npart_new,info)

             CASE (ppm_type_comp,ppm_type_comp_single)
                CALL ppm_map_part_pop(prop%data_1d_c,Pc%Npart,Npart_new,info)

             CASE (ppm_type_logical)
                CALL ppm_map_part_pop(prop%data_1d_l,Pc%Npart,Npart_new,info)

             END SELECT

          ENDIF
          or_fail("ppm_map_part_pop")

          end_subroutine()
      END SUBROUTINE DTYPE(part_prop_pop)

      SUBROUTINE DTYPE(part_prop_ghost_pop)(Pc,prop_id,Mpart,info)

          !-------------------------------------------------------------------------
          !  Modules
          !-------------------------------------------------------------------------
          USE ppm_module_map
          IMPLICIT NONE

          DEFINE_MK()
          !-------------------------------------------------------------------------
          ! Include
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(DTYPE(ppm_t_particles)) :: Pc
          !!! Data structure containing the particles
          INTEGER,        INTENT(IN   ) :: prop_id
          !!! id of the property to be pooped
          INTEGER,        INTENT(IN   ) :: Mpart
          !!! The number of real + ghost particles (on the processor) to pop from the buffer
          INTEGER,        INTENT(  OUT) :: info
          !!! Returns status, 0 upon success.
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          CLASS(DTYPE(ppm_t_part_prop)_), POINTER :: prop

          INTEGER :: lda

          start_subroutine("map_part_ghost_pop")

          !-----------------------------------------------------------------
          !  Call ppm_map_part_ghost_pop
          !-----------------------------------------------------------------
          prop => Pc%props%vec(prop_id)%t

          lda = prop%lda

          IF (lda.GE.2) THEN
             SELECT CASE (prop%data_type)
             CASE (ppm_type_int)
                CALL ppm_map_part_ghost_pop(prop%data_2d_i,lda,Pc%Npart,Mpart,info)

             CASE (ppm_type_longint)
                fail('Type not supported for mappings.')

             CASE (ppm_type_real,ppm_type_real_single)
                CALL ppm_map_part_ghost_pop(prop%data_2d_r,lda,Pc%Npart,Mpart,info)

             CASE (ppm_type_comp,ppm_type_comp_single)
                CALL ppm_map_part_ghost_pop(prop%data_2d_c,lda,Pc%Npart,Mpart,info)

             CASE (ppm_type_logical)
                CALL ppm_map_part_ghost_pop(prop%data_2d_l,lda,Pc%Npart,Mpart,info)

             END SELECT
          ELSE
             SELECT CASE (prop%data_type)
             CASE (ppm_type_int)
                CALL ppm_map_part_ghost_pop(prop%data_1d_i,lda,Pc%Npart,Mpart,info)

             CASE (ppm_type_longint)
                fail('Type not supported for mappings.')

             CASE (ppm_type_real,ppm_type_real_single)
                CALL ppm_map_part_ghost_pop(prop%data_1d_r,lda,Pc%Npart,Mpart,info)

             CASE (ppm_type_comp,ppm_type_comp_single)
                CALL ppm_map_part_ghost_pop(prop%data_1d_c,lda,Pc%Npart,Mpart,info)

             CASE (ppm_type_logical)
                CALL ppm_map_part_ghost_pop(prop%data_1d_l,lda,Pc%Npart,Mpart,info)

             END SELECT

          ENDIF
          or_fail("ppm_map_part_ghost_pop")

          end_subroutine()
      END SUBROUTINE DTYPE(part_prop_ghost_pop)

      !ghost mapping
      !
      SUBROUTINE DTYPE(part_map_ghost_get)(Pc,info,ghostsize)
          !!! This routine maps/adds the ghost particles on the current topology.
          !!! This routine is similar to the partial mapping routine
          !!! (`ppm_map_part_partial`) in the sense that the ghost particles are
          !!! assumed to be located on neighbouring processors only, and thus only
          !!! require a nearest neighbour communication.
          !-------------------------------------------------------------------------
          !  Modules
          !-------------------------------------------------------------------------
          USE ppm_module_time, ONLY : ppm_time
          USE ppm_module_map
          IMPLICIT NONE

          DEFINE_MK()
          !-------------------------------------------------------------------------
          ! Include
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(DTYPE(ppm_t_particles))     :: Pc

          !!! Data structure containing the particles
          INTEGER,            INTENT(  OUT) :: info
          !!! Return status, on success 0.
          !-------------------------------------------------------------------------
          !  Optional Arguments
          !-------------------------------------------------------------------------
          REAL(MK), OPTIONAL, INTENT(IN   ) :: ghostsize
          !!! size of the ghost layers. Default is to use the particles cutoff
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          TYPE(ppm_t_topo), POINTER :: topo

          REAL(MK)              :: cutoff
          !!! cutoff radius
          REAL(ppm_kind_double) :: t1,t2

          INTEGER :: topoid
          !!! index variable

          start_subroutine("part_map_ghost_get")

          !-----------------------------------------------------------------
          !  Checks
          !-----------------------------------------------------------------
          !check that particles are allocated
          check_associated(<#Pc%xp#>,&
          & "Particles structure had not been defined. Call allocate first")
          !check that particles are mapped onto this topology
          check_true(<#Pc%flags(ppm_part_partial)#>,&
          & "Partial/global mapping required before doing a ghost mapping")
          !check that particles are inside the domain
          check_true(<#Pc%flags(ppm_part_areinside)#>,&
          & "Some particles may be outside the domain. Apply BC first")

          topoid = Pc%active_topoid

          topo => ppm_topo(topoid)%t

          cutoff = Pc%ghostlayer
          IF (PRESENT(ghostsize)) THEN
             IF (ghostsize .LT. cutoff) THEN
                fail("using ghostsize < cutoff+skin. Increase ghostsize.")
             ELSE
                cutoff = ghostsize
             ENDIF
          ENDIF

#if   __KIND == __SINGLE_PRECISION
          IF (cutoff.GT.topo%ghostsizes) THEN
             stdout("cutoff for ghost mapping       = ",cutoff)
             stdout("cutoff used to create topology = ",'topo%ghostsizes')
             fail("ghostsize of topology may be smaller than that of particles")
          ENDIF
#elif   __KIND == __DOUBLE_PRECISION
          IF (cutoff .GT. topo%ghostsized) THEN
             stdout("cutoff for ghost mapping       = ",cutoff)
             stdout("cutoff used to create topology = ",'topo%ghostsized')
             fail("ghostsize of topology may be smaller than that of particles")
          ENDIF
#endif
          IF (cutoff.GT.0._MK) THEN
             IF (Pc%flags(ppm_part_ghosts)) THEN
                IF (ppm_map_type_isactive(ppm_param_map_ghost_get).AND. &
                &   ppm_buffer_set.EQ.1) THEN
                   !TOCHECK
                   !Yaser : I have made this change as to avoid extra
                   ! expensive ghost_get when it is already available.
                   !skipping ghost get
                ELSE
                   Pc%stats%nb_ghost_get = Pc%stats%nb_ghost_get + 1

                   CALL ppm_time(t1,info)

                   CALL ppm_map_part_ghost_get(topoid,Pc%xp,ppm_dim, &
                   &    Pc%Npart,Pc%isymm,cutoff,info)
                   or_fail("ppm_map_part_ghost_get failed")

                   CALL ppm_time(t2,info)

                   Pc%stats%t_ghost_get = Pc%stats%t_ghost_get + (t2-t1)
                ENDIF
             ENDIF !Pc%flags(ppm_part_ghosts)
          ENDIF !cutoff.GT.0._MK

          end_subroutine()
      END SUBROUTINE DTYPE(part_map_ghost_get)

      SUBROUTINE DTYPE(part_map_ghost_put)(Pc,info)
          !!! Put back ghost particle values/properties to the
          !!! corresponding real particles. This is very useful in the case of
          !!! symmetric interactions as there the ghost particles are also updated.
          !!!
          !!! [IMPORTANT]
          !!! This routine can only be called after ghost particles have been
          !!! created using the `ppm_map_part_ghost_get` (+push/send/pop sequence)
          !!! routine.
          !-------------------------------------------------------------------------
          !  Modules
          !-------------------------------------------------------------------------
          USE ppm_module_map
          IMPLICIT NONE

          DEFINE_MK()
          !-------------------------------------------------------------------------
          ! Include
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(DTYPE(ppm_t_particles)) :: Pc

          !!! Data structure containing the particles
          INTEGER,        INTENT(  OUT) :: info
          !!! Return status, on success 0.
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          INTEGER :: topoid
          !!! index variable

          start_subroutine("part_map_ghost_put")

          !-----------------------------------------------------------------
          !  Checks
          !-----------------------------------------------------------------
          check_true(<#ppm_map_type_isactive(ppm_param_map_ghost_get)#>, &
          & "Need to call ghost_get (+push/send/pop sequence) before ghost_put")

          topoid = Pc%active_topoid

          IF (topoid.GT.0) THEN
             CALL ppm_map_part_ghost_put(topoid,info)
             or_fail("ppm_map_part_ghost_put failed")
          ENDIF

          end_subroutine()
      END SUBROUTINE DTYPE(part_map_ghost_put)

      SUBROUTINE DTYPE(part_map_ghost_push)(Pc,info,Field)
          !!! Push ghost particles properties into the send buffer.
          !!! If ghost positions are not up-to-date, they are re-computed
          !!! with a call to ghost_get. Otherwise ghost_get is skipped.
          !!! If the Field argument is present, only the discretization of that
          !!! Field on this particle set is pushed. Otherwise, all the properties
          !!! (that are not already up-to-date) are pushed, except those that
          !!! have the ppm_ppt_map_ghosts flag set to false.
          !!!  Assumptions:
          !!! * Particles positions need to have been mapped onto the topology
          !!!
          !-------------------------------------------------------------------------
          !  Modules
          !-------------------------------------------------------------------------
          USE ppm_module_time, ONLY : ppm_time
          USE ppm_module_map
          IMPLICIT NONE

          DEFINE_MK()
          !-------------------------------------------------------------------------
          ! Include
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(DTYPE(ppm_t_particles))                    :: Pc

          !!! Data structure containing the particles
          INTEGER,                           INTENT(  OUT) :: info
          !!! Return status, on success 0.
          !-------------------------------------------------------------------------
          !  Optional Arguments
          !-------------------------------------------------------------------------
          CLASS(ppm_t_main_abstr), OPTIONAL, TARGET        :: Field
          !!! Push only this field. Default is to push all of them.
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          CLASS(ppm_t_discr_data), POINTER :: prop

          !!! cutoff radius
          REAL(ppm_kind_double) :: t1,t2

          start_subroutine("part_map_ghost_push")

          !-----------------------------------------------------------------
          !  Checks
          !-----------------------------------------------------------------
          !check that particles are allocated
          check_associated(<#Pc%xp#>,&
          & "Particles structure had not been defined. Call allocate first")
          !check that particles are mapped onto this topology
          check_true(<#Pc%flags(ppm_part_partial)#>,&
          & "Do a partial/global mapping before doing a ghost mapping")
          !check that particles are inside the domain
          check_true(<#Pc%flags(ppm_part_areinside)#>,&
          & "some particles may be outside the domain. Apply BC first")

          IF (Pc%ghostlayer.GT.0._MK) THEN
             check_true(<#(ppm_map_type_isactive(ppm_param_map_ghost_get).OR.ppm_map_type_isactive(ppm_param_map_ghost_put))#>, &
             & "Need to call ghost_get before ghost_push")

             IF (PRESENT(Field)) THEN

                SELECT TYPE (Field)
                CLASS IS (ppm_t_field_)
                   NULLIFY(prop)
                   !Get discretization
                   !Points the iterator (props%iter_id) to the discretization
                   !of Field on the particle set Pc.
                   CALL Pc%get_discr(Field,prop,info)
                   or_fail("could not get discretization for this Field")

                   IF (.NOT.prop%flags(ppm_ppt_map_ghosts)) THEN
                      fail("Cannot map ghosts for a property for which the flag ppm_ppt_map_ghosts is not true")
                   ENDIF

                CLASS IS (DTYPE(ppm_t_part_prop)_)
                   Pc%props%iter_id=Pc%props%get_id(Field)

                   IF (.NOT.Field%flags(ppm_ppt_map_ghosts)) THEN
                      fail("Cannot map ghosts for a property for which the flag ppm_ppt_map_ghosts is not true")
                   ENDIF

                END SELECT

                Pc%stats%nb_ghost_push = Pc%stats%nb_ghost_push + 1

                CALL ppm_time(t1,info)

                CALL Pc%map_part_push_legacy(Pc%props%iter_id,info)
                or_fail("map_part_push")

                CALL ppm_time(t2,info)

                Pc%stats%t_ghost_push = Pc%stats%t_ghost_push + (t2-t1)

             ELSE
                !Update the ghost for the properties if
                ! 1) they have been mapped to this topology,
                ! 2) the ghosts have not yet been updated, and
                ! 3) the user wants them to be updated
                prop => Pc%props%begin()
                DO WHILE (ASSOCIATED(prop))
                   IF (prop%flags(ppm_ppt_map_ghosts)) THEN
                      IF (.NOT.prop%flags(ppm_ppt_ghosts)) THEN
                         IF (prop%flags(ppm_ppt_partial)) THEN

                            Pc%stats%nb_ghost_push = Pc%stats%nb_ghost_push + 1

                            CALL ppm_time(t1,info)

                            CALL Pc%map_part_push_legacy(Pc%props%iter_id,info)
                            or_fail("map_part_push")

                            CALL ppm_time(t2,info)

                            Pc%stats%t_ghost_push = Pc%stats%t_ghost_push + (t2-t1)

                         ELSE
                            stdout("pushing property ",'Pc%props%iter_id','TRIM(prop%name)')
                            fail("getting ghost for a property thats not mapped")
                         ENDIF !prop%flags(ppm_ppt_partial)
                      ENDIF !.NOT.prop%flags(ppm_ppt_ghosts)
                   ENDIF !prop%flags(ppm_ppt_map_ghosts)
                   prop => Pc%props%next()
                ENDDO ! ASSOCIATED(prop)
             ENDIF !PRESENT(Field)

          ELSE ! if cutoff .le. 0
             !stdout("cutoff = 0, nothing to do")
          ENDIF ! Pc%ghostlayer .GT. 0._MK

          end_subroutine()
      END SUBROUTINE DTYPE(part_map_ghost_push)

      SUBROUTINE DTYPE(part_map_ghost_send)(Pc,info)

          !!!  Send buffers during ghost mapping for particles
          !!!  Assumptions:
          !!! * Particles positions need to have been mapped onto the topology
          !!!
          !-------------------------------------------------------------------------
          !  Modules
          !-------------------------------------------------------------------------
          USE ppm_module_time, ONLY : ppm_time
          USE ppm_module_map
          IMPLICIT NONE

          DEFINE_MK()
          !-------------------------------------------------------------------------
          ! Include
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(DTYPE(ppm_t_particles)) :: Pc

          !!! Data structure containing the particles
          INTEGER,        INTENT(  OUT) :: info
          !!! Return status, on success 0.
          !-------------------------------------------------------------------------
          !  Optional Arguments
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          REAL(ppm_kind_double) :: t1,t2

          INTEGER :: Mpart

          start_subroutine("part_map_ghost_send")

          !-----------------------------------------------------------------
          !  Checks
          !-----------------------------------------------------------------
          !check that particles are allocated
          check_associated(<#Pc%xp#>, &
          & "Particles structure had not been defined. Call allocate first")
          !check that particles are mapped onto this topology
          check_true(<#Pc%flags(ppm_part_partial)#>, &
          & "Do a partial/global mapping before doing a ghost mapping")

          Pc%stats%nb_ghost_send = Pc%stats%nb_ghost_send + 1

          IF (ppm_map_type_isactive(ppm_param_map_ghost_get)) THEN
             CALL ppm_time(t1,info)

             !-----------------------------------------------------------------
             !  Send the buffer
             !-----------------------------------------------------------------
             CALL ppm_map_part_send(Pc%Npart,Pc%Mpart,info)
             or_fail("ppm_map_part_send")

             CALL ppm_time(t2,info)
          ELSE IF (ppm_map_type_isactive(ppm_param_map_ghost_put)) THEN
             CALL ppm_time(t1,info)

             !-----------------------------------------------------------------
             !  Send the buffer
             !-----------------------------------------------------------------
             CALL ppm_map_part_send(Pc%Npart,Mpart,info)
             or_fail("ppm_map_part_send")

             CALL ppm_time(t2,info)
          ENDIF

          Pc%stats%t_ghost_send = Pc%stats%t_ghost_send + (t2-t1)

          end_subroutine()
      END SUBROUTINE DTYPE(part_map_ghost_send)

      SUBROUTINE DTYPE(part_map_ghost_isend)(Pc,info,sendrecv)

          !!!  Send buffers during ghost mapping for particles
          !!!  Assumptions:
          !!! * Particles positions need to have been mapped onto the topology
          !!!
          !-------------------------------------------------------------------------
          !  Modules
          !-------------------------------------------------------------------------
          USE ppm_module_time, ONLY : ppm_time
          USE ppm_module_map
          IMPLICIT NONE

          DEFINE_MK()
          !-------------------------------------------------------------------------
          ! Include
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(DTYPE(ppm_t_particles))    :: Pc

          !!! Data structure containing the particles
          INTEGER,           INTENT(  OUT) :: info
          !!! Return status, on success 0.
          LOGICAL, OPTIONAL, INTENT(IN   ) :: sendrecv
          !-------------------------------------------------------------------------
          !  Optional Arguments
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          REAL(ppm_kind_double) :: t1,t2

          INTEGER :: Mpart

          start_subroutine("part_map_ghost_send")

          !-----------------------------------------------------------------
          !  Checks
          !-----------------------------------------------------------------
          IF (PRESENT(sendrecv)) THEN
             IF (sendrecv) THEN
                !check that particles are allocated
                check_associated(<#Pc%xp#>, &
                & "Particles structure had not been defined. Call allocate first")
                !check that particles are mapped onto this topology
                check_true(<#Pc%flags(ppm_part_partial)#>, &
                & "Do a partial/global mapping before doing a ghost mapping")

                Pc%stats%nb_ghost_send = Pc%stats%nb_ghost_send + 1

                IF (ppm_map_type_isactive(ppm_param_map_ghost_get)) THEN

                   CALL ppm_time(t1,info)
                   !-----------------------------------------------------------------
                   !  Send the buffer
                   !-----------------------------------------------------------------
                   CALL ppm_map_part_isend(Pc%Npart,Pc%Mpart,info,sendrecv=.TRUE.)
                   or_fail("ppm_map_part_isend")

                   CALL ppm_time(t2,info)
                ELSE IF (ppm_map_type_isactive(ppm_param_map_ghost_put)) THEN

                   CALL ppm_time(t1,info)
                   !-----------------------------------------------------------------
                   !  Send the buffer
                   !-----------------------------------------------------------------
                   CALL ppm_map_part_isend(Pc%Npart,Mpart,info,sendrecv=.TRUE.)
                   or_fail("ppm_map_part_isend")

                   CALL ppm_time(t2,info)
                ENDIF
             ELSE
                IF (ppm_map_type_isactive(ppm_param_map_ghost_get)) THEN

                   CALL ppm_time(t1,info)
                   !-----------------------------------------------------------------
                   !  Send the buffer
                   !-----------------------------------------------------------------
                   CALL ppm_map_part_isend(Pc%Npart,Pc%Mpart,info,sendrecv=.FALSE.)
                   or_fail("ppm_map_part_isend")

                   CALL ppm_time(t2,info)
                ELSE IF (ppm_map_type_isactive(ppm_param_map_ghost_put)) THEN

                   CALL ppm_time(t1,info)
                   !-----------------------------------------------------------------
                   !  Send the buffer
                   !-----------------------------------------------------------------
                   CALL ppm_map_part_isend(Pc%Npart,Mpart,info,sendrecv=.FALSE.)
                   or_fail("ppm_map_part_isend")

                   CALL ppm_time(t2,info)
                ENDIF
             ENDIF
          ELSE
             !check that particles are allocated
             check_associated(<#Pc%xp#>, &
             & "Particles structure had not been defined. Call allocate first")
             !check that particles are mapped onto this topology
             check_true(<#Pc%flags(ppm_part_partial)#>, &
             & "Do a partial/global mapping before doing a ghost mapping")

             Pc%stats%nb_ghost_send = Pc%stats%nb_ghost_send + 1

             IF (ppm_map_type_isactive(ppm_param_map_ghost_get)) THEN

                CALL ppm_time(t1,info)
                !-----------------------------------------------------------------
                !  Send the buffer
                !-----------------------------------------------------------------
                CALL ppm_map_part_isend(Pc%Npart,Pc%Mpart,info)
                or_fail("ppm_map_part_isend")

                CALL ppm_time(t2,info)
             ELSE IF (ppm_map_type_isactive(ppm_param_map_ghost_put)) THEN

                CALL ppm_time(t1,info)
                !-----------------------------------------------------------------
                !  Send the buffer
                !-----------------------------------------------------------------
                CALL ppm_map_part_isend(Pc%Npart,Mpart,info)
                or_fail("ppm_map_part_isend")

                CALL ppm_time(t2,info)
             ENDIF
          ENDIF

          Pc%stats%t_ghost_send = Pc%stats%t_ghost_send + (t2-t1)

          end_subroutine()
      END SUBROUTINE DTYPE(part_map_ghost_isend)

      SUBROUTINE DTYPE(part_map_ghost_pop)(Pc,info,Field)
          !!! This routine pops the contents of the receive buffer
          !!! If mapping is ghost_put when sending the ghosts back,
          !!! the value popped will be *added* to array passed to the routine.
          !-------------------------------------------------------------------------
          !  Modules
          !-------------------------------------------------------------------------
          USE ppm_module_mpi
          USE ppm_module_map
          IMPLICIT NONE

          DEFINE_MK()
          !-------------------------------------------------------------------------
          ! Include
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(DTYPE(ppm_t_particles))                    :: Pc

          !!! Data structure containing the particles
          INTEGER,                           INTENT(  OUT) :: info
          !!! Return status, on success 0.
          !-------------------------------------------------------------------------
          !  Optional Arguments
          !-------------------------------------------------------------------------
          CLASS(ppm_t_main_abstr), OPTIONAL, TARGET        :: Field
          !!! Pop only this field. Default is to pop all of them
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          CLASS(ppm_t_discr_data), POINTER :: prop

          start_subroutine("part_map_ghost_pop")

          !-----------------------------------------------------------------
          !  Checks
          !-----------------------------------------------------------------
          !check that particles are allocated
          check_associated(<#Pc%xp#>, &
          & "Particles structure had not been defined. Call allocate first")
          !check that particles are mapped onto this topology
          check_true(<#Pc%flags(ppm_part_partial)#>, &
          & "Do a partial/global mapping before doing a ghost mapping")

          IF (Pc%ghostlayer.GT.0._MK) THEN
             IF (ppm_map_type_isactive(ppm_param_map_ghost_put)) THEN
                IF (PRESENT(Field)) THEN
                   SELECT TYPE (Field)
                   CLASS IS (ppm_t_field_)
                      NULLIFY(prop)
                      !Get discretization
                      !Points the iterator (props%iter_id) to the discretization
                      !of Field on the particle set Pc.
                      CALL Pc%get_discr(Field,prop,info)
                      or_fail("could not get discretization for this Field")

                      IF (.NOT.prop%flags(ppm_ppt_map_ghosts)) THEN
                         fail("Cannot map ghosts for a property for which the flag ppm_ppt_map_ghosts is not true")
                      ENDIF

                   CLASS IS (DTYPE(ppm_t_part_prop)_)
                      Pc%props%iter_id=Pc%props%get_id(Field)

                      IF (.NOT.Field%flags(ppm_ppt_map_ghosts)) THEN
                         fail("Cannot map ghosts for a property for which the flag ppm_ppt_map_ghosts is not true")
                      ENDIF

                      prop => Pc%props%at(Pc%props%iter_id)

                   END SELECT

                   CALL Pc%map_part_ghost_pop_legacy(Pc%props%iter_id,Pc%Mpart,info)
                   or_fail("map_part_ghost_pop_legacy")

                   prop%flags(ppm_ppt_ghosts) = .FALSE.
                ELSE
                   prop => Pc%props%last()
                   DO WHILE (ASSOCIATED(prop))
                      IF (prop%flags(ppm_ppt_map_ghosts)) THEN
                         IF (.NOT.prop%flags(ppm_ppt_ghosts)) THEN
                            IF (prop%flags(ppm_ppt_partial)) THEN
                               CALL Pc%map_part_ghost_pop_legacy(Pc%props%iter_id,Pc%Mpart,info)
                               or_fail("map_part_ghost_pop_legacy")

                               prop%flags(ppm_ppt_ghosts) = .FALSE.
                            ENDIF !prop%flags(ppm_ppt_partial)
                         ENDIF !.NOT.prop%flags(ppm_ppt_ghosts)
                      ENDIF !prop%flags(ppm_ppt_map_ghosts)
                      prop => Pc%props%prev()
                   ENDDO !ASSOCIATED(prop)
                ENDIF !PRESENT(Field)
             ELSE IF (Pc%flags(ppm_part_ghosts).OR. &
             &   ppm_map_type_isactive(ppm_param_map_ghost_get)) THEN
                IF (PRESENT(Field)) THEN
                   SELECT TYPE (Field)
                   CLASS IS (ppm_t_field_)
                      NULLIFY(prop)
                      !Get discretization
                      !Points the iterator (props%iter_id) to the discretization
                      !of Field on the particle set Pc.
                      CALL Pc%get_discr(Field,prop,info)
                      or_fail("could not get discretization for this Field")

                      IF (.NOT.prop%flags(ppm_ppt_map_ghosts)) THEN
                         fail("Cannot map ghosts for a property for which the flag ppm_ppt_map_ghosts is not true")
                      ENDIF

                   CLASS IS (DTYPE(ppm_t_part_prop)_)
                      Pc%props%iter_id=Pc%props%get_id(Field)

                      IF (.NOT.Field%flags(ppm_ppt_map_ghosts)) THEN
                         fail("Cannot map ghosts for a property for which the flag ppm_ppt_map_ghosts is not true")
                      ENDIF

                      prop => Pc%props%at(Pc%props%iter_id)

                   END SELECT

                   CALL Pc%map_part_pop_legacy(Pc%props%iter_id,Pc%Mpart,info)
                   or_fail("map_part_pop")

                   prop%flags(ppm_ppt_ghosts) = .TRUE.
                ELSE
                   prop => Pc%props%last()
                   DO WHILE (ASSOCIATED(prop))
                      IF (prop%flags(ppm_ppt_map_ghosts)) THEN
                         IF (.NOT.prop%flags(ppm_ppt_ghosts)) THEN
                            IF (prop%flags(ppm_ppt_partial)) THEN
                               CALL Pc%map_part_pop_legacy(Pc%props%iter_id,Pc%Mpart,info)
                               or_fail("map_part_pop")

                               prop%flags(ppm_ppt_ghosts) = .TRUE.
                            ENDIF !prop%flags(ppm_ppt_partial)
                         ENDIF !.NOT.prop%flags(ppm_ppt_ghosts)
                      ENDIF !prop%flags(ppm_ppt_map_ghosts)
                      prop => Pc%props%prev()
                   ENDDO !ASSOCIATED(prop)
                ENDIF !PRESENT(Field)
             ELSE
                fail("Ghost buffer invalid. Correct sequence is ghost_get, ghost_push, ghost_send and ghost_pop.")
             ENDIF
          ELSE ! if cutoff .le. 0
             ! Update states
             prop => Pc%props%begin()
             DO WHILE (ASSOCIATED(prop))
                IF (prop%flags(ppm_ppt_map_ghosts)) THEN
                   prop%flags(ppm_ppt_ghosts) = .TRUE.
                ENDIF
                prop => Pc%props%next()
             ENDDO
          ENDIF !Pc%ghostlayer.GT.0._MK

          end_subroutine()
      END SUBROUTINE DTYPE(part_map_ghost_pop)

      SUBROUTINE DTYPE(part_map_ghost_push_pos)(Pc,info)
          !!! Push ghost particles positions from the send buffer.
          !!!
          !-------------------------------------------------------------------------
          !  Modules
          !-------------------------------------------------------------------------
          USE ppm_module_map
          IMPLICIT NONE

          DEFINE_MK()
          !-------------------------------------------------------------------------
          ! Include
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(DTYPE(ppm_t_particles)) :: Pc

          !!! Data structure containing the particles
          INTEGER,        INTENT(  OUT) :: info
          !!! Return status, on success 0.
          !-------------------------------------------------------------------------
          !  Optional Arguments
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------

          start_subroutine("part_map_ghost_push_pos")

          !-----------------------------------------------------------------
          !  Checks
          !-----------------------------------------------------------------
          !check that particles are allocated
          check_associated(<#Pc%xp#>, &
          & "Particles structure had not been defined. Call allocate first")
          !check that particles are mapped onto this topology
          check_true(<#Pc%flags(ppm_part_partial)#>, &
          & "Do a partial/global mapping before doing a ghost mapping")

          IF (Pc%ghostlayer .GT. 0._MK) THEN
             check_true(<#Pc%flags(ppm_part_ghosts)#>,&
             & "flags(ppm_part_ghosts) need to be set to .true.")

             check_true(<#ppm_map_type_isactive(ppm_param_map_ghost_get)#>, &
             & "Ghost buffer invalid. Correct sequence is ghost_get, ghost_push, ghost_send and ghost_pop.")

             CALL ppm_map_part_push(Pc%xp,ppm_dim,Pc%Npart,info,pushpp=.TRUE.)
             or_fail("map_part_push")

          ELSE ! if cutoff .le. 0
             !stdout("cutoff = 0, nothing to do")
             !stdout("setting all %has_ghost properties to true")
          ENDIF

          end_subroutine()
      END SUBROUTINE DTYPE(part_map_ghost_push_pos)

      SUBROUTINE DTYPE(part_map_ghost_pop_pos)(Pc,info)
          !!! Pop ghost particles positions from the send buffer.
          !!!  Assumptions:
          !!! * Needs to be called at the end of the
          !!!  ghost_get, ghost_push, send and ghost_pop call sequence.
          !!!
          !-------------------------------------------------------------------------
          !  Modules
          !-------------------------------------------------------------------------
          USE ppm_module_map
          IMPLICIT NONE

          DEFINE_MK()
          !-------------------------------------------------------------------------
          ! Include
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(DTYPE(ppm_t_particles)) :: Pc

          !!! Data structure containing the particles
          INTEGER,      INTENT(  OUT)   :: info
          !!! Return status, on success 0.
          !-------------------------------------------------------------------------
          !  Optional Arguments
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------

          start_subroutine("part_map_ghost_pop_pos")

          !-----------------------------------------------------------------
          !  Checks
          !-----------------------------------------------------------------
          !check that particles are allocated
          check_associated(<#Pc%xp#>, &
          & "Particles structure had not been defined. Call allocate first")

          !check that particles are mapped onto this topology
          check_true(<#Pc%flags(ppm_part_partial)#>, &
          & "Do a partial/global mapping before doing a ghost mapping")

          IF (Pc%ghostlayer.GT.0._MK) THEN
             check_true(<#Pc%flags(ppm_part_ghosts)#>,&
             & "flags(ppm_part_ghosts) need to be set to .true.")

             check_true(<#(ppm_map_type_isactive(ppm_param_map_ghost_get).AND.ppm_buffer_set.EQ.1)#>,&
             & "Ghost buffer invalid. Correct sequence is ghost_get, ghost_push, ghost_send and ghost_pop.")

             CALL ppm_map_part_pop(Pc%xp,ppm_dim,Pc%Npart,Pc%Mpart,info)
             or_fail("map_part_pop")
          ELSE ! if cutoff .le. 0
             !stdout("cutoff = 0, nothing to do")
             !stdout("setting all %has_ghost properties to true")
          ENDIF

          ! Update states
          !   ghosts have been computed
          Pc%flags(ppm_part_ghosts) = .TRUE.

          end_subroutine()
      END SUBROUTINE DTYPE(part_map_ghost_pop_pos)

      SUBROUTINE DTYPE(part_map_ghosts)(Pc,info,ghostsize)

          !!!  Ghost mapping for particles
          !!!  Assumptions:
          !!! * Particles positions need to have been mapped onto the topology
          !!!
          !-------------------------------------------------------------------------
          !  Modules
          !-------------------------------------------------------------------------
          USE ppm_module_time, ONLY : ppm_time
          USE ppm_module_map
          IMPLICIT NONE

          DEFINE_MK()
          !-------------------------------------------------------------------------
          ! Include
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(DTYPE(ppm_t_particles))      :: Pc

          !!! Data structure containing the particles
          INTEGER,             INTENT(  OUT) :: info
          !!! Return status, on success 0.
          !-------------------------------------------------------------------------
          !  Optional Arguments
          !-------------------------------------------------------------------------
          REAL(MK),  OPTIONAL, INTENT(IN   ) :: ghostsize
          !!! size of the ghost layers. Default is to use the particles cutoff
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          CLASS(ppm_t_discr_data), POINTER :: prop

          TYPE(ppm_t_topo), POINTER :: topo

          REAL(MK)              :: cutoff
          !!! cutoff radius
          REAL(ppm_kind_double) :: t1,t2

          INTEGER :: topoid
          !!! index variable

          LOGICAL :: skip_ghost_get
          LOGICAL :: skip_send

          start_subroutine("part_map_ghosts")

          skip_ghost_get = .FALSE.
          skip_send = .TRUE.
          !we must not call ppm_map_part_send unless ppm_map_part_push (or ghost_get)
          !has been called (in which case, skip_send is set to FALSE)

          !-----------------------------------------------------------------
          !  Checks
          !-----------------------------------------------------------------
          !check that particles are allocated
          check_associated(<#Pc%xp#>, &
          & "Particles structure had not been defined. Call allocate first")
          !check that particles are mapped onto this topology
          check_true(<#Pc%flags(ppm_part_partial)#>, &
          & "Do a partial/global mapping before doing a ghost mapping")
          !check that particles are inside the domain
          check_true(<#Pc%flags(ppm_part_areinside)#>, &
          & "some particles may be outside the domain. Apply BC first")

          topoid=Pc%active_topoid

          topo => ppm_topo(topoid)%t

          cutoff = Pc%ghostlayer
          IF (PRESENT(ghostsize)) THEN
             IF (ghostsize.LT.cutoff) THEN
                fail("using ghostsize < cutoff+skin. Increase ghostsize.")
             ELSE
                cutoff = ghostsize
             ENDIF
          ENDIF

#if   __KIND == __SINGLE_PRECISION
          IF (cutoff .GT. topo%ghostsizes) THEN
             stdout("cutoff for ghost mapping       = ",cutoff)
             stdout("cutoff used to create topology = ",'topo%ghostsizes')
             fail("ghostsize of topology may be smaller than that of particles")
          ENDIF
#elif   __KIND == __DOUBLE_PRECISION
          IF (cutoff .GT. topo%ghostsized) THEN
             stdout("cutoff for ghost mapping       = ",cutoff)
             stdout("cutoff used to create topology = ",'topo%ghostsized')
             fail("ghostsize of topology may be smaller than that of particles")
          ENDIF
#endif

          IF (cutoff .GT. 0._MK) THEN
             IF (Pc%flags(ppm_part_ghosts)) THEN
                IF (ppm_map_type_isactive(ppm_param_map_ghost_get)) THEN
                   skip_ghost_get = .TRUE.
                ENDIF
             ENDIF

             IF (skip_ghost_get) THEN
                !skip ghost get
             ELSE
                Pc%stats%nb_ghost_get = Pc%stats%nb_ghost_get + 1

                CALL ppm_time(t1,info)

                CALL ppm_map_part_ghost_get(topoid,Pc%xp,ppm_dim,&
                &    Pc%Npart,Pc%isymm,cutoff,info)
                or_fail("ppm_map_part_ghost_get failed")

                CALL ppm_time(t2,info)

                Pc%stats%t_ghost_get = Pc%stats%t_ghost_get + (t2-t1)

                skip_send = .FALSE.
             ENDIF

             !Update the ghost for the properties if
             ! 1) they have been mapped to this topology,
             ! 2) the ghosts have not yet been updated, and
             ! 3) the user wants them to be updated
             prop => Pc%props%begin()
             DO WHILE (ASSOCIATED(prop))
                IF (prop%flags(ppm_ppt_map_ghosts)) THEN
                   IF (.NOT.prop%flags(ppm_ppt_ghosts)) THEN
                      IF (prop%flags(ppm_ppt_partial)) THEN
                         Pc%stats%nb_ghost_push = Pc%stats%nb_ghost_push + 1

                         CALL ppm_time(t1,info)

                         CALL Pc%map_part_push_legacy(Pc%props%iter_id,info)
                         or_fail("map_part_push")

                         CALL ppm_time(t2,info)

                         Pc%stats%t_ghost_push=Pc%stats%t_ghost_push+(t2-t1)

                         skip_send = .FALSE.
                      ELSE
                         stdout("pushing property ",'Pc%props%iter_id','TRIM(prop%name)')
                         fail("getting ghosts for a property thats not mapped")
                      ENDIF
                   ENDIF
                ENDIF

                prop => Pc%props%next()
             ENDDO

             IF (.NOT.skip_send) THEN
                CALL ppm_map_part_isend(Pc%Npart,Pc%Mpart,info)
                or_fail("ppm_map_part_isend")

                prop => Pc%props%last()
                DO WHILE (ASSOCIATED(prop))

                   IF (prop%flags(ppm_ppt_map_ghosts)) THEN
                      IF (.NOT.prop%flags(ppm_ppt_ghosts)) THEN
                         IF (prop%flags(ppm_ppt_partial)) THEN
                            CALL Pc%map_part_pop_legacy(Pc%props%iter_id,Pc%Mpart,info)
                            or_fail("map_part_pop")

                            prop%flags(ppm_ppt_ghosts) = .TRUE.
                         ENDIF
                      ENDIF
                   ENDIF

                   prop => Pc%props%prev()
                ENDDO

                IF (.NOT.skip_ghost_get) THEN
                   CALL ppm_map_part_pop(Pc%xp,ppm_dim,Pc%Npart,Pc%Mpart,info)
                   or_fail("map_part_pop")
                ENDIF
             ENDIF !.NOT.skip_send

          ELSE ! if cutoff .le. 0
             !stdout("cutoff = 0, nothing to do")
             !stdout("setting all %has_ghost properties to true")

             prop => Pc%props%begin()
             DO WHILE (ASSOCIATED(prop))
                IF (prop%flags(ppm_ppt_map_ghosts)) THEN
                   prop%flags(ppm_ppt_ghosts) = .TRUE.
                ENDIF
                prop => Pc%props%next()
             ENDDO
          ENDIF
          ! Update states
          ! ghosts have been computed
          Pc%flags(ppm_part_ghosts) = .TRUE.
          ! the states for the properties have already been updated above

          end_subroutine()
      END SUBROUTINE DTYPE(part_map_ghosts)

      SUBROUTINE DTYPE(part_map)(Pc,info,global,topoid)
          !!!  Partial/Global mapping for particles
          !!!  Assumptions:
          !!! * All the particles have to be inside the domain
          !!!   (otherwise -> "unassigned particle error")
          !-------------------------------------------------------------------------
          !  Modules
          !-------------------------------------------------------------------------
          USE ppm_module_time, ONLY : ppm_time
          USE ppm_module_map
          IMPLICIT NONE

          DEFINE_MK()
          !-------------------------------------------------------------------------
          ! Include
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(DTYPE(ppm_t_particles))    :: Pc
          !!! Data structure containing the particles
          INTEGER,           INTENT(  OUT) :: info
          !!! Return status, on success 0.
          !-------------------------------------------------------------------------
          !  Optional Arguments
          !-------------------------------------------------------------------------
          LOGICAL, OPTIONAL, INTENT(IN   ) :: global
          !!! does a global mapping. Default is false (i.e. partial mapping)
          INTEGER, OPTIONAL, INTENT(IN   ) :: topoid
          !!! topology id
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          CLASS(DTYPE(ppm_t_part_prop)_), POINTER :: prop
          CLASS(DTYPE(ppm_t_neighlist)_), POINTER :: nl
          CLASS(ppm_t_operator_discr_),   POINTER :: op

          REAL(ppm_kind_double) :: t1,t2

          INTEGER :: Npart_new
          !!! new number of particles on this processor
          INTEGER :: ltopoid
          !!! index variable

          LOGICAL :: partial

          start_subroutine("particles_mapping")

          !-----------------------------------------------------------------
          !  Checks
          !-----------------------------------------------------------------
          check_associated(<#Pc%xp#>,&
          & "Particles structure had not been defined. Call allocate first")
          !check that particles are inside the domain
          check_true(<#Pc%flags(ppm_part_areinside)#>,&
          & "some particles may be outside the domain. Apply BC first")

          partial = .TRUE.
          IF (PRESENT(global)) THEN
             IF (global) THEN
                IF (.NOT.PRESENT(topoid)) THEN
                   fail("need the topoid parameter for global mapping")
                ENDIF
                partial = .FALSE.
             ENDIF
          ENDIF

          ltopoid=MERGE(Pc%active_topoid,topoid,partial)

          CALL ppm_time(t1,info)

          !-----------------------------------------------------------------------
          !  Map the particles onto the topology
          !-----------------------------------------------------------------------
          IF (partial.AND.Pc%flags(ppm_part_partial)) THEN
             !Particles have already been mapped onto this topology
             !nothing to do
          ELSE
             IF (partial) THEN
                CALL ppm_map_part_partial(ltopoid,Pc%xp,Pc%Npart,info)
                or_fail("ppm_map_part_partial")

                Pc%stats%nb_part_map = Pc%stats%nb_part_map + 1
             ELSE
                CALL ppm_map_part_global(ltopoid,Pc%xp,Pc%Npart,info)
                or_fail("ppm_map_part_global")

                Pc%stats%nb_global_map = Pc%stats%nb_global_map + 1
             ENDIF

             prop => Pc%props%begin()
             DO WHILE (ASSOCIATED(prop))
                IF (prop%flags(ppm_ppt_map_parts)) THEN
                   CALL Pc%map_part_push_legacy(Pc%props%iter_id,info)
                   or_fail("ppm_map_part_push_legacy")
                ENDIF
                prop => Pc%props%next()
             ENDDO

             CALL ppm_map_part_isend(Pc%Npart,Npart_new,info)
             or_fail("ppm_map_part_isend")

             prop => Pc%props%last()
             DO WHILE (ASSOCIATED(prop))
                IF (prop%flags(ppm_ppt_map_parts)) THEN
                   CALL Pc%map_part_pop_legacy(Pc%props%iter_id,Npart_new,info)
                   or_fail("ppm_map_part_pop_legacy")

                   prop%flags(ppm_ppt_partial) = .TRUE.
                ENDIF
                prop => Pc%props%prev()
             ENDDO

             CALL ppm_map_part_pop(Pc%xp,ppm_dim,Pc%Npart,Npart_new,info)
             or_fail("ppm_map_part_prop")

             ! Update states
             ! Number of particles on this processor
             Pc%Npart = Npart_new
             Pc%Mpart = Pc%Npart

             ! This is the active topology for these particles
             IF (.NOT.partial) Pc%active_topoid = topoid

             ! Particles are now mapped on the active topology
             Pc%flags(ppm_part_partial) = .TRUE.
             ! Pc have been re-indexed and ghosts have not been computed
             Pc%flags(ppm_part_ghosts) = .FALSE.

             !values for poperty arrays have been mapped and ghosts
             !are no longer up-to-date
             prop => Pc%props%begin()
             DO WHILE (ASSOCIATED(prop))
                prop%flags(ppm_ppt_ghosts) = .FALSE.
                prop => Pc%props%next()
             ENDDO

             ! particles have been re-indexed and neighbour lists not updated
             nl => Pc%neighs%begin()
             DO WHILE (ASSOCIATED(nl))
                nl%uptodate = .FALSE.
                nl => Pc%neighs%next()
             ENDDO
             Pc%flags(ppm_part_neighlists) = .FALSE.

             ! particles have been re-indexed and operators need be recomputed
             IF (ASSOCIATED(Pc%ops)) THEN
                op => Pc%ops%begin()
                DO WHILE (ASSOCIATED(op))
                   op%flags(ppm_ops_iscomputed) = .FALSE.
                   op => Pc%ops%next()
                ENDDO
             ENDIF

          ENDIF !(partial .AND. Pc%flags(ppm_part_partial))

          CALL ppm_time(t2,info)

          IF (partial) THEN
             Pc%stats%t_part_map = Pc%stats%t_part_map + (t2-t1)
          ELSE
             Pc%stats%t_global_map = Pc%stats%t_global_map + (t2-t1)
          ENDIF

          end_subroutine()
      END SUBROUTINE DTYPE(part_map)

      SUBROUTINE DTYPE(part_map_positions)(Pc,info,global,topoid)
          !!!  Partial/Global mapping for particles. Only computes the mappings
          !!! and puts them into the send buffers.
          !!!  Assumptions:
          !!! * All the particles have to be inside the domain
          !!!   (otherwise -> "unassigned particle error")
          !-------------------------------------------------------------------------
          !  Modules
          !-------------------------------------------------------------------------
          USE ppm_module_map
          IMPLICIT NONE

          DEFINE_MK()
          !-------------------------------------------------------------------------
          ! Include
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(DTYPE(ppm_t_particles))    :: Pc

          !!! Data structure containing the particles
          INTEGER,           INTENT(  OUT) :: info
          !!! Return status, on success 0.
          !-------------------------------------------------------------------------
          !  Optional Arguments
          !-------------------------------------------------------------------------
          LOGICAL, OPTIONAL, INTENT(IN   ) :: global
          !!! does a global mapping. Default is false (i.e. partial mapping)
          INTEGER, OPTIONAL, INTENT(IN   ) :: topoid
          !!! topology id
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          INTEGER :: ltopoid
          !!! index variable

          LOGICAL :: partial

          start_subroutine("part_map_positions")

          !-----------------------------------------------------------------
          !  Checks
          !-----------------------------------------------------------------
          check_associated(<#Pc%xp#>,&
          & "Particles structure had not been defined. Call allocate first")
          !check that particles are inside the domain
          check_true(<#Pc%flags(ppm_part_areinside)#>,&
          & "some particles may be outside the domain. Apply BC first")

          partial = .TRUE.
          IF (PRESENT(global)) THEN
             IF (global) THEN
                IF (.NOT.PRESENT(topoid)) THEN
                   fail("need the topoid parameter for global mapping")
                ENDIF
                partial = .FALSE.
             ENDIF
          ENDIF

          ltopoid=MERGE(Pc%active_topoid,topoid,partial)

          !-----------------------------------------------------------------------
          !  Map the particles onto the topology
          !-----------------------------------------------------------------------
          IF (partial .AND. Pc%flags(ppm_part_partial)) THEN
              !Particles have already been mapped onto this topology
              !nothing to do
          ELSE
             IF (partial) THEN
                CALL ppm_map_part_partial(ltopoid,Pc%xp,Pc%Npart,info)
                or_fail("ppm_map_part_partial")

                Pc%stats%nb_part_map = Pc%stats%nb_part_map + 1
             ELSE
                CALL ppm_map_part_global(ltopoid,Pc%xp,Pc%Npart,info)
                or_fail("ppm_map_part_global")

                Pc%stats%nb_global_map = Pc%stats%nb_global_map + 1
             ENDIF
          ENDIF

          Pc%active_topoid = ltopoid

          end_subroutine()
      END SUBROUTINE DTYPE(part_map_positions)


      SUBROUTINE DTYPE(part_map_push)(Pc,info,Field)
          !!! Push particles properties into the send buffer after a partial mapping
          !!!  Assumptions:
          !!! * All the particles have to be inside the domain
          !!!   (otherwise -> "unassigned particle error")
          !-------------------------------------------------------------------------
          !  Modules
          !-------------------------------------------------------------------------
          USE ppm_module_map
          IMPLICIT NONE

          DEFINE_MK()
          !-------------------------------------------------------------------------
          ! Include
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(DTYPE(ppm_t_particles))                    :: Pc

          !!! Data structure containing the particles
          INTEGER,                           INTENT(  OUT) :: info
          !!! Return status, on success 0.
          !-------------------------------------------------------------------------
          !  Optional Arguments
          !-------------------------------------------------------------------------
          CLASS(ppm_t_main_abstr), OPTIONAL, TARGET        :: Field
          !!! Push only this field. Default is to push all of them.
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          CLASS(ppm_t_discr_data), POINTER :: prop

          start_subroutine("part_map_push")

          !-----------------------------------------------------------------
          !  Checks
          !-----------------------------------------------------------------
          check_associated(<#Pc%xp#>,&
          & "Particles structure had not been defined. Call allocate first")
          !check that particles are inside the domain
          check_true(<#Pc%flags(ppm_part_areinside)#>,&
          & "some particles may be outside the domain. Apply BC first")

          check_true(<#(ppm_map_type_isactive(ppm_param_map_global).OR.ppm_map_type_isactive(ppm_param_map_partial))#>, &
          & "Need to call map_positions before map_push")

          !-----------------------------------------------------------------------
          !  push properties into the send buffer
          !-----------------------------------------------------------------------
          IF (PRESENT(Field)) THEN
             SELECT TYPE (Field)
             CLASS IS (ppm_t_field_)
                NULLIFY(prop)
                !Get discretization
                !Points the iterator (props%iter_id) to the discretization
                !of Field on the particle set Pc.
                CALL Pc%get_discr(Field,prop,info)
                or_fail("could not get discretization for this Field")

                IF (.NOT.prop%flags(ppm_ppt_map_parts)) THEN
                   fail("Cannot map a property for which the flag ppm_ppt_map_parts is not true")
                ENDIF

             CLASS IS (DTYPE(ppm_t_part_prop)_)
                Pc%props%iter_id=Pc%props%get_id(Field)

                IF (.NOT.Field%flags(ppm_ppt_map_parts)) THEN
                   fail("Cannot map a property for which the flag ppm_ppt_map_parts is not true")
                ENDIF

             END SELECT

             CALL Pc%map_part_push_legacy(Pc%props%iter_id,info)
             or_fail("map_part_push")
          ELSE
             !push all properties into the send buffer
             prop => Pc%props%begin()
             DO WHILE (ASSOCIATED(prop))
                IF (prop%flags(ppm_ppt_map_parts)) THEN
                   CALL Pc%map_part_push_legacy(Pc%props%iter_id,info)
                   or_fail("map_part_push")
                ENDIF
                prop => Pc%props%next()
             ENDDO
          ENDIF

          end_subroutine()
      END SUBROUTINE DTYPE(part_map_push)

      SUBROUTINE DTYPE(part_map_send)(Pc,info)

          !!!  Send buffer during partial mappings
          !!!  Assumptions:
          !!! * Particles positions need to have been mapped onto the topology
          !!!
          USE ppm_module_map
          IMPLICIT NONE

          DEFINE_MK()
          !-------------------------------------------------------------------------
          ! Include
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(DTYPE(ppm_t_particles)) :: Pc

          !!! Data structure containing the particles
          INTEGER,        INTENT(  OUT) :: info
          !!! Return status, on success 0.

          start_subroutine("part_map_send")

          !-----------------------------------------------------------------
          !  Checks
          !-----------------------------------------------------------------
          !check that particles are allocated
          check_associated(<#Pc%xp#>,&
          & "Particles structure had not been defined. Call allocate first")

          !-----------------------------------------------------------------
          !  Send the buffer
          !-----------------------------------------------------------------
          CALL ppm_map_part_send(Pc%Npart,Pc%NewNpart,info)
          or_fail("ppm_map_part_send")

          end_subroutine()
      END SUBROUTINE DTYPE(part_map_send)

      SUBROUTINE DTYPE(part_map_isend)(Pc,info,sendrecv)
          !!!  Send buffer during partial mappings
          !!!  Assumptions:
          !!! * Particles positions need to have been mapped onto the topology
          !!!
          USE ppm_module_map
          IMPLICIT NONE

          DEFINE_MK()
          !-------------------------------------------------------------------------
          ! Include
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(DTYPE(ppm_t_particles))    :: Pc

          !!! Data structure containing the particles
          INTEGER,           INTENT(  OUT) :: info
          !!! Return status, on success 0.

          LOGICAL, OPTIONAL, INTENT(IN   ) :: sendrecv

          start_subroutine("part_map_isend")

          IF (PRESENT(sendrecv)) THEN
             IF (sendrecv) THEN
                !-----------------------------------------------------------------
                !  Checks
                !-----------------------------------------------------------------
                !check that particles are allocated
                check_associated(<#Pc%xp#>,&
                & "Particles structure had not been defined. Call allocate first")

                !-----------------------------------------------------------------
                !  Send the buffer
                !-----------------------------------------------------------------
                CALL ppm_map_part_isend(Pc%Npart,Pc%NewNpart,info,sendrecv=.TRUE.)
                or_fail("ppm_map_part_isend")
             ELSE
                !-----------------------------------------------------------------
                !  Send the buffer
                !-----------------------------------------------------------------
                CALL ppm_map_part_isend(Pc%Npart,Pc%NewNpart,info,sendrecv=.FALSE.)
                or_fail("ppm_map_part_isend")
             ENDIF
          ELSE
             !-----------------------------------------------------------------
             !  Checks
             !-----------------------------------------------------------------
             !check that particles are allocated
             check_associated(<#Pc%xp#>,&
             & "Particles structure had not been defined. Call allocate first")

             !-----------------------------------------------------------------
             !  Send the buffer
             !-----------------------------------------------------------------
             CALL ppm_map_part_isend(Pc%Npart,Pc%NewNpart,info)
             or_fail("ppm_map_part_isend")
          ENDIF

          end_subroutine()
      END SUBROUTINE DTYPE(part_map_isend)


      SUBROUTINE DTYPE(part_map_pop)(Pc,info,Field)
          !-------------------------------------------------------------------------
          !  Modules
          !-------------------------------------------------------------------------
          USE ppm_module_map
          IMPLICIT NONE

          DEFINE_MK()
          !-------------------------------------------------------------------------
          ! Include
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(DTYPE(ppm_t_particles))                    :: Pc

          !!! Data structure containing the particles
          INTEGER,                           INTENT(  OUT) :: info
          !!! Return status, on success 0.
          !-------------------------------------------------------------------------
          !  Optional Arguments
          !-------------------------------------------------------------------------
          CLASS(ppm_t_main_abstr), OPTIONAL, TARGET        :: Field
          !!! Pop only this field. Default is to pop all of them
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          CLASS(ppm_t_discr_data), POINTER :: prop

          start_subroutine("part_map_pop")

          !-----------------------------------------------------------------
          !  Checks
          !-----------------------------------------------------------------
          check_associated(<#Pc%xp#>,&
          & "Particles structure had not been defined. Call allocate first")
          !check that particles are inside the domain
          check_true(<#Pc%flags(ppm_part_areinside)#>,&
          & "some particles may be outside the domain. Apply BC first")

          !-----------------------------------------------------------------------
          !  Pop properties out of the buffer
          !-----------------------------------------------------------------------
          IF (PRESENT(Field)) THEN
             SELECT TYPE (Field)
             CLASS IS (ppm_t_field_)
                NULLIFY(prop)
                !Get discretization
                !Points the iterator (props%iter_id) to the discretization
                !of Field on the particle set Pc.
                CALL Pc%get_discr(Field,prop,info)
                or_fail("could not get discretization for this Field")

                IF (.NOT.prop%flags(ppm_ppt_map_parts)) THEN
                   fail("Cannot map a property for which the flag ppm_ppt_map_parts is not true")
                ENDIF

             CLASS IS (DTYPE(ppm_t_part_prop)_)
                Pc%props%iter_id=Pc%props%get_id(Field)

                IF (.NOT.Field%flags(ppm_ppt_map_parts)) THEN
                   fail("Cannot map a property for which the flag ppm_ppt_map_parts is not true")
                ENDIF

                prop => Pc%props%at(Pc%props%iter_id)

             END SELECT

             CALL Pc%map_part_pop_legacy(Pc%props%iter_id,Pc%NewNpart,info)
             or_fail("map_part_pop")

             prop%flags(ppm_ppt_partial) = .TRUE.
          ELSE
             prop => Pc%props%last()
             DO WHILE (ASSOCIATED(prop))
                IF (prop%flags(ppm_ppt_map_parts)) THEN
                   CALL Pc%map_part_pop_legacy(Pc%props%iter_id,Pc%NewNpart,info)
                   or_fail("map_part_pop_legacy")

                   prop%flags(ppm_ppt_partial) = .TRUE.
                ENDIF
                prop => Pc%props%prev()
             ENDDO
          ENDIF

          end_subroutine()
      END SUBROUTINE DTYPE(part_map_pop)


      SUBROUTINE DTYPE(part_map_pop_positions)(Pc,info)
          !-------------------------------------------------------------------------
          !  Modules
          !-------------------------------------------------------------------------
          USE ppm_module_map
          IMPLICIT NONE

          DEFINE_MK()
          !-------------------------------------------------------------------------
          ! Include
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(DTYPE(ppm_t_particles)) :: Pc

          !!! Data structure containing the particles
          INTEGER,        INTENT(  OUT) :: info
          !!! Return status, on success 0.
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          CLASS(DTYPE(ppm_t_part_prop)_), POINTER :: prop
          CLASS(DTYPE(ppm_t_neighlist)_), POINTER :: nl
          CLASS(ppm_t_operator_discr_),   POINTER :: op

          start_subroutine("part_map_pop_positions")

          !-----------------------------------------------------------------
          !  Checks
          !-----------------------------------------------------------------
          check_associated(<#Pc%xp#>,&
          & "Particles structure had not been defined. Call allocate first")

          check_true(<#(ppm_map_type_isactive(ppm_param_map_global).OR.ppm_map_type_isactive(ppm_param_map_partial))#>, &
          & "Need to call map_positions before map_push")

          CALL ppm_map_part_pop(Pc%xp,ppm_dim,Pc%Npart,Pc%NewNpart,info)
          or_fail("map_part_pop")

          ! Update states
          ! Number of particles on this processor
          Pc%Npart = Pc%NewNpart
          Pc%Mpart = Pc%Npart

          ! Particles are now mapped on the active topology
          Pc%flags(ppm_part_partial) = .TRUE.
          ! Pc have been re-indexed and ghosts have not been computed
          Pc%flags(ppm_part_ghosts) = .FALSE.

          !   values for poperty arrays have been mapped and ghosts
          !   are no longer up-to-date
          prop => Pc%props%begin()
          DO WHILE (ASSOCIATED(prop))
             prop%flags(ppm_ppt_ghosts) = .FALSE.
             prop => Pc%props%next()
          ENDDO

          ! particles have been re-indexed and neighbour lists not updated
          nl => Pc%neighs%begin()
          DO WHILE (ASSOCIATED(nl))
             nl%uptodate = .FALSE.
             nl => Pc%neighs%next()
          ENDDO
          Pc%flags(ppm_part_neighlists) = .FALSE.

          ! particles have been re-indexed and operators need be recomputed
          IF (ASSOCIATED(Pc%ops)) THEN
             op => Pc%ops%begin()
             DO WHILE (ASSOCIATED(op))
                op%flags(ppm_ops_iscomputed) = .FALSE.
                op => Pc%ops%next()
             ENDDO
          ENDIF

          end_subroutine()
      END SUBROUTINE DTYPE(part_map_pop_positions)

      SUBROUTINE DTYPE(part_map_store)(Pc,info)
          !!! This routine stores the state of the current/active particle mapping.
          !!!
          !!! This is especially useful for alternating ghost-get/ghost-put
          !!! mappings, which are used for implementing the communication in
          !!! symmetric particle interactions. The `ppm_map_part_ghost_get` mapping
          !!! can then be stored after its first call and reloaded after each
          !!! ghost-put mapping.
          !TODO
          !This routine should be moved to the mapping when it is complete
          !-------------------------------------------------------------------------
          !  Modules
          !-------------------------------------------------------------------------
          USE ppm_module_map
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          ! Include
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(DTYPE(ppm_t_particles)) :: Pc
          !!! Data structure containing the particles

          INTEGER,        INTENT(  OUT) :: info
          !!! Return status, on success 0.
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------

          start_subroutine("part_map_store")

          CALL ppm_map_part_store(info)
          or_fail("ppm_map_part_store")

          end_subroutine()
      END SUBROUTINE DTYPE(part_map_store)

      SUBROUTINE DTYPE(part_map_load)(Pc,info)
          !!! This routine loads the internally stored particle mapping.
          !!!
          !!! See `ppm_map_part_store` for more information.
          !TODO
          !This routine should be moved to the mapping when it is complete
          !-------------------------------------------------------------------------
          !  Modules
          !-------------------------------------------------------------------------
          USE ppm_module_map
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          ! Include
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(DTYPE(ppm_t_particles)) :: Pc
          !!! Data structure containing the particles

          INTEGER,        INTENT(  OUT) :: info
          !!! Return status, on success 0.
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------

          start_subroutine("part_map_load")

          CALL ppm_map_part_load(info)
          or_fail("ppm_map_part_load")

          end_subroutine()
      END SUBROUTINE DTYPE(part_map_load)

      !Yaser
      !TOCHECK
      !destroy the link in particle data
      !DESTROY ENTRY
      SUBROUTINE DTYPE(neigh_destroy)(neigh,info)
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(DTYPE(ppm_t_neighlist)) :: neigh

          INTEGER,       INTENT(  OUT) :: info
          !!! Returns status, 0 upon success.
          !-------------------------------------------------------------------------
          ! local variables
          !-------------------------------------------------------------------------
          INTEGER :: ldu(2)

          start_subroutine("neigh_destroy")

          IF (ASSOCIATED(neigh%nvlist)) THEN
             CALL ppm_alloc(neigh%nvlist,ldu,ppm_param_dealloc,info)
             or_fail_dealloc("neigh%nvlist")
          ENDIF
          IF (ASSOCIATED(neigh%vlist)) THEN
             CALL ppm_alloc(neigh%vlist,ldu,ppm_param_dealloc,info)
             or_fail_dealloc("neigh%vlist")
          ENDIF

          NULLIFY(neigh%nvlist)
          NULLIFY(neigh%vlist)
          neigh%name=''
          NULLIFY(neigh%Part)
          neigh%uptodate = .FALSE.
          neigh%nneighmin=0
          neigh%nneighmax=0

          end_subroutine()
      END SUBROUTINE DTYPE(neigh_destroy)
