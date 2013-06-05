      MODULE pointer_tracker
         ! Hashmap
         ! transate incomming data into a character array
         ! hash, check if associated (linked list)
         ! return identifier
         USE ppm_module_core
         TYPE abstr_ptr
            CHARACTER(LEN=10) :: hash
            CLASS(ppm_t_main_abstr), POINTER :: ptr
            TYPE(abstr_ptr), POINTER :: next => null()
         END TYPE abstr_ptr
         INTERFACE get_pointer
            MODULE PROCEDURE get_abstr_pointer, get_mapping_pointer, &
                  get_particles_stats_pointer
         END INTERFACE get_pointer
         CONTAINS
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
            CHARACTER(32) FUNCTION get_abstr_pointer(some_ptr)
               CLASS(ppm_t_main_abstr), INTENT(in) :: some_ptr
               CHARACTER(LEN=1), DIMENSION(:), ALLOCATABLE :: chars
               CHARACTER(LEN=32) :: phex
               CHARACTER(LEN=20) :: format_
               INTEGER length

               length = size(transfer(some_ptr, chars))
               chars = transfer(some_ptr, chars)

               WRITE(format_, '("(", I0, "Z2.2)")') length
               WRITE (phex, format_) chars
               !WRITE(*,*) trim(phex)
               get_abstr_pointer=phex
            END FUNCTION get_abstr_pointer
            CHARACTER(32) FUNCTION get_mapping_pointer(some_ptr)
               CLASS(ppm_t_mapping_d_), INTENT(in) :: some_ptr
               CHARACTER(LEN=1), DIMENSION(:), ALLOCATABLE :: chars
               CHARACTER(LEN=32) :: phex
               CHARACTER(LEN=20) :: format_
               INTEGER length

               length = size(transfer(some_ptr, chars))
               ALLOCATE (chars(length))
               chars = transfer(some_ptr, chars)

               WRITE(format_, '("(", I0, "Z2.2)")') length
               WRITE (phex, format_) chars
               !WRITE(*,*) trim(phex)
               get_mapping_pointer=phex
            END FUNCTION get_mapping_pointer
            CHARACTER(32) FUNCTION get_particles_stats_pointer(some_ptr)
               CLASS(particles_stats_d_), INTENT(in) :: some_ptr
               CHARACTER(LEN=1), DIMENSION(:), ALLOCATABLE :: chars
               CHARACTER(LEN=32) :: phex
               CHARACTER(LEN=20) :: format_
               INTEGER length

               length = size(transfer(some_ptr, chars))
               ALLOCATE (chars(length))
               chars = transfer(some_ptr, chars)

               WRITE(format_, '("(", I0, "Z2.2)")') length
               WRITE (phex, format_) chars
               !WRITE(*,*) trim(phex)
               get_particles_stats_pointer=phex
            END FUNCTION get_particles_stats_pointer
            CHARACTER(32) FUNCTION get_real1_pointer(some_ptr)
               REAL, DIMENSION(:), POINTER :: some_ptr
               CHARACTER(LEN=1), DIMENSION(:), ALLOCATABLE :: chars
               CHARACTER(LEN=32) :: phex
               CHARACTER(LEN=20) :: format_
               INTEGER length

               length = size(transfer(some_ptr, chars))
               ALLOCATE (chars(length))
               chars = transfer(some_ptr, chars)

               WRITE(format_, '("(", I0, "Z2.2)")') length
               WRITE (phex, format_) chars
               !WRITE(*,*) trim(phex)
               get_real1_pointer=phex
            END FUNCTION get_real1_pointer

      END MODULE
