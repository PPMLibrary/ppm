      MODULE pointer_tracker
         ! Hashmap
         ! transate incomming data into a character array
         ! hash, check if associated (linked list)
         ! return identifier
         USE ppm_module_core
         TYPE abstr_ptr
            CHARACTER(LEN=10) :: hash
            INTEGER, DIMENSION(:), POINTER :: ptr
            TYPE(abstr_ptr), POINTER :: next => null()
         END TYPE abstr_ptr
         TYPE ptr_map
            INTEGER ::  nextpointer = 1
            TYPE(abstr_ptr) :: map(100)
            !CONTAINS
            !   PROCEDURE :: insert
            !   PROCEDURE :: lookup
         END TYPE ptr_map
         INCLUDE 'pointers/ppm_pointers_interface.f'
         INCLUDE 'maps/pmap_typedef.f'
         INCLUDE 'maps/primitive_maps.f'
         !INTERFACE get_pointer
         !   MODULE PROCEDURE get_ppm_pointer
            !MODULE PROCEDURE get_abstr_pointer, get_mapping_pointer, &
            !      get_particles_stats_pointer, &
            !      get_integer1d_pointer, &
            !      get_integer2d_pointer
         !END INTERFACE get_pointer
         CONTAINS
            INCLUDE 'pointers/ppm_pointers.f'
            INCLUDE 'maps/primitive_maps_procedures.f'
            INTEGER FUNCTION FNVHash(int_str, numwords)
               IMPLICIT NONE
               !INTEGER(8), PARAMETER :: PRIME_8 = 14695981039346656037
               !INTEGER(8), PARAMETER :: OFFSET_8 = 1099511628211
               INTEGER, PARAMETER :: PRIME = 16777619
               INTEGER, PARAMETER :: OFFSET = -2128831035
               INTEGER, DIMENSION(:), INTENT(in) :: int_str
               INTEGER, INTENT(in):: numwords
               INTEGER :: x
               INTEGER :: i

               x = PRIME
               DO i=1, numwords
                  x = IEOR(x, int_str(i))
                  x = x * OFFSET
               ENDDO
               FNVHash = x

            END FUNCTION FNVHash
            INTEGER FUNCTION get_integer1d_pointer(some_ptr)
               INTEGER, DIMENSION(:), POINTER :: some_ptr
               INTEGER, DIMENSION(:), POINTER :: buffer
               INTEGER length

               length = size(transfer(some_ptr, buffer))
               ALLOCATE (buffer(length))
               buffer = transfer(some_ptr, buffer)

               WRITE (*,*) length
               IF (associated(buffer, TARGET=some_ptr)) THEN
                  WRITE(*,*) "We can use this"
               ELSE
                  WRITE(*,*) "We need to 1337 h4cks"
               ENDIF

               get_integer1d_pointer=FNVHash(buffer, length)
            END FUNCTION get_integer1d_pointer

            INTEGER FUNCTION get_integer2d_pointer(some_ptr)
               INTEGER, DIMENSION(:,:), POINTER :: some_ptr
               INTEGER, DIMENSION(:), POINTER :: buffer

               length = size(transfer(some_ptr,buffer))
               ALLOCATE (buffer(length))
               buffer = transfer(some_ptr, buffer)

               get_integer2d_pointer=FNVHash(buffer, size(buffer))
            END FUNCTION get_integer2d_pointer

      END MODULE
