      MODULE pointer_tracker
         ! Hashmap
         ! transate incomming data into a character array
         ! hash, check if associated (linked list)
         ! return identifier
         USE ppm_module_core
         INCLUDE 'pointers/ppm_pointers_interface.f'
         INCLUDE 'trees/tree_abstract_typedef.f'
         TYPE(pointer_trees), SAVE :: pointer_data
         CONTAINS
            INCLUDE 'pointers/ppm_pointers.f'
            INCLUDE 'trees/tree_abstract.f'
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
            INCLUDE 'intrinsic/get_pointers.f'

      END MODULE
