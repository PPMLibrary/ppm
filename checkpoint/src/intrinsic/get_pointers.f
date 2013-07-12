
            CHARACTER(LEN=32) FUNCTION get_integer1d_pointer(some_ptr)
               INTEGER, DIMENSION(:), POINTER :: some_ptr
               INTEGER, DIMENSION(:), POINTER :: buffer
               INTEGER length
               TYPE(integer1d_tree), POINTER :: nodelett

               length = size(transfer(some_ptr, buffer))
               ALLOCATE (buffer(length))
               buffer = transfer(some_ptr, buffer)

               ALLOCATE (nodelett)
               nodelett%hash = FNVHash(buffer, length)
               nodelett%val => some_ptr
               DEALLOCATE (buffer)

               CALL pointer_insert(pointer_data, nodelett, &
                  get_integer1d_pointer)
            END FUNCTION get_integer1d_pointer

            CHARACTER(LEN=32) FUNCTION get_integer2d_pointer(some_ptr)
               INTEGER, DIMENSION(:,:), POINTER :: some_ptr
               INTEGER, DIMENSION(:), POINTER :: buffer
               INTEGER length
               TYPE(integer2d_tree), POINTER :: nodelett

               length = size(transfer(some_ptr,buffer))
               ALLOCATE (buffer(length))
               buffer = transfer(some_ptr, buffer)

               ALLOCATE(nodelett)
               nodelett%hash = FNVHash(buffer, length)
               nodelett%val => some_ptr
               DEALLOCATE (buffer)

               CALL pointer_insert(pointer_data, nodelett, &
                  get_integer2d_pointer)
            END FUNCTION get_integer2d_pointer
            CHARACTER(LEN=32) FUNCTION get_real1d_pointer(some_ptr)
               REAL(ppm_kind_double), DIMENSION(:), POINTER :: some_ptr
               INTEGER, DIMENSION(:), POINTER :: buffer
               INTEGER length
               TYPE(real1d_tree), POINTER :: nodelett

               length = size(transfer(some_ptr, buffer))
               ALLOCATE (buffer(length))
               buffer = transfer(some_ptr, buffer)

               ALLOCATE (nodelett)
               nodelett%hash = FNVHash(buffer, length)
               nodelett%val => some_ptr
               DEALLOCATE (buffer)

               CALL pointer_insert(pointer_data, nodelett, &
                  get_real1d_pointer)
            END FUNCTION get_real1d_pointer

            CHARACTER(LEN=32) FUNCTION get_real2d_pointer(some_ptr)
               REAL(ppm_kind_double), DIMENSION(:,:), POINTER :: some_ptr
               INTEGER, DIMENSION(:), POINTER :: buffer
               INTEGER length
               TYPE(real2d_tree), POINTER :: nodelett
               !TYPE(real2d_tree) :: nodelett

               length = size(transfer(some_ptr,buffer))
               ALLOCATE (buffer(length))
               buffer = transfer(some_ptr, buffer)

               ALLOCATE(nodelett)
               nodelett%hash = FNVHash(buffer, length)
               nodelett%val => some_ptr
               DEALLOCATE (buffer)

               CALL pointer_insert(pointer_data, nodelett, &
                  get_real2d_pointer)
            END FUNCTION get_REAL2d_pointer

            CHARACTER(LEN=32) FUNCTION get_logical1d_pointer(some_ptr)
               LOGICAL, DIMENSION(:), POINTER :: some_ptr
               INTEGER, DIMENSION(:), POINTER :: buffer
               INTEGER length
               TYPE(logical1d_tree), POINTER :: nodelett

               length = size(transfer(some_ptr,buffer))
               ALLOCATE (buffer(length))
               buffer = transfer(some_ptr, buffer)

               ALLOCATE(nodelett)
               nodelett%hash = FNVHash(buffer, length)
               nodelett%val => some_ptr
               DEALLOCATE (buffer)

               CALL pointer_insert(pointer_data, nodelett, &
                  get_logical1d_pointer)
            END FUNCTION get_logical1d_pointer

            CHARACTER(LEN=32) FUNCTION get_logical2d_pointer(some_ptr)
               LOGICAL, DIMENSION(:,:), POINTER :: some_ptr
               INTEGER, DIMENSION(:), POINTER :: buffer
               INTEGER length
               TYPE(logical2d_tree), POINTER :: nodelett

               length = size(transfer(some_ptr,buffer))
               ALLOCATE (buffer(length))
               buffer = transfer(some_ptr, buffer)

               ALLOCATE(nodelett)
               nodelett%hash = FNVHash(buffer, length)
               nodelett%val => some_ptr
               DEALLOCATE (buffer)

               CALL pointer_insert(pointer_data, nodelett, &
                  get_logical2d_pointer)
            END FUNCTION get_logical2d_pointer

            CHARACTER(LEN=32) FUNCTION get_integer64_1d_pointer(some_ptr)
               INTEGER(8), DIMENSION(:), POINTER :: some_ptr
               INTEGER, DIMENSION(:), POINTER :: buffer
               INTEGER length
               TYPE(integer64_1d_tree), POINTER :: nodelett

               length = size(transfer(some_ptr,buffer))
               ALLOCATE (buffer(length))
               buffer = transfer(some_ptr, buffer)

               ALLOCATE(nodelett)
               nodelett%hash = FNVHash(buffer, length)
               nodelett%val => some_ptr
               DEALLOCATE (buffer)

               CALL pointer_insert(pointer_data, nodelett, &
                  get_integer64_1d_pointer)
            END FUNCTION get_integer64_1d_pointer

            CHARACTER(LEN=32) FUNCTION get_integer64_2d_pointer(some_ptr)
               INTEGER(8), DIMENSION(:,:), POINTER :: some_ptr
               INTEGER, DIMENSION(:), POINTER :: buffer
               INTEGER length
               TYPE(integer64_2d_tree), POINTER :: nodelett

               ALLOCATE(nodelett)
               length = size(transfer(some_ptr,buffer))
               ALLOCATE (buffer(length))
               buffer = transfer(some_ptr, buffer)

               nodelett%hash = FNVHash(buffer, length)
               nodelett%val => some_ptr
               DEALLOCATE (buffer)

               CALL pointer_insert(pointer_data, nodelett, &
                  get_integer64_2d_pointer)
            END FUNCTION get_integer64_2d_pointer

            CHARACTER(LEN=32) FUNCTION get_complex1d_pointer(some_ptr)
               COMPLEX(8), DIMENSION(:), POINTER :: some_ptr
               INTEGER, DIMENSION(:), POINTER :: buffer
               INTEGER length
               TYPE(complex1d_tree), POINTER :: nodelett

               length = size(transfer(some_ptr,buffer))
               ALLOCATE (buffer(length))
               buffer = transfer(some_ptr, buffer)

               ALLOCATE(nodelett)
               nodelett%hash = FNVHash(buffer, length)
               nodelett%val => some_ptr
               DEALLOCATE (buffer)

               CALL pointer_insert(pointer_data, nodelett, &
                  get_complex1d_pointer)
            END FUNCTION get_complex1d_pointer

            CHARACTER(LEN=32) FUNCTION get_complex2d_pointer(some_ptr)
               COMPLEX(8), DIMENSION(:,:), POINTER :: some_ptr
               INTEGER, DIMENSION(:), POINTER :: buffer
               INTEGER length
               TYPE(complex2d_tree), POINTER :: nodelett

               length = size(transfer(some_ptr,buffer))
               ALLOCATE (buffer(length))
               buffer = transfer(some_ptr, buffer)

               ALLOCATE(nodelett)
               nodelett%hash = FNVHash(buffer, length)
               nodelett%val => some_ptr
               DEALLOCATE (buffer)

               CALL pointer_insert(pointer_data, nodelett, &
                  get_complex2d_pointer)
            END FUNCTION get_complex2d_pointer
