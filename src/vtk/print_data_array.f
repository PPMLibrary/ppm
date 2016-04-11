#ifndef VTK_FORMAT
#define VTK_FORMAT "ascii"
#endif
#ifndef VTK_TYPE
#define VTK_TYPE "Float64"
#endif
           WRITE(iUnit,'(A)',ADVANCE='NO') "        <DataArray"
#ifdef VTK_TYPE
           WRITE(iUnit,'(3A)',ADVANCE='NO') " type='", VTK_TYPE, "'"
#endif
#ifdef VTK_NAME
           WRITE(iUnit,'(3A)',ADVANCE='NO') " Name='", &
           & VTK_NAME(1:LEN_TRIM(VTK_NAME)), "'"
#endif
#ifdef VTK_NDIM
           WRITE(scratch, *) VTK_NDIM
           scratch = ADJUSTL(scratch)
           WRITE(iUnit,'(3A)',ADVANCE='NO') " NumberOfComponents='", &
           & scratch(1:LEN_TRIM(scratch)), "'"
#endif
           WRITE(iUnit,'(3A)',ADVANCE='NO') " format='", VTK_FORMAT, "'"
#ifdef VTK_OFFSET
           WRITE(scratch, *) VTK_OFFSET
           scratch = ADJUSTL(scratch)
           WRITE(iUnit,'(3A)',ADVANCE='NO') " offset='", &
           & scratch(1:LEN_TRIM(scratch)), "'"
#endif
           WRITE(iUnit,'(A)') ">"
#ifdef VTK_RANGE
#ifndef VTK_RANGE_START
#define VTK_RANGE_START 1
#endif
           DO i=VTK_RANGE_START,VTK_RANGE
              WRITE(iUnit,'(I0)',ADVANCE='NO') i
              IF (i .LT. VTK_RANGE) WRITE(iUnit,'(A)',ADVANCE='NO') " "
           ENDDO
#else
#ifdef VTK_MESH
#ifndef VTK_MESH_ILBOUND
#define VTK_MESH_ILBOUND 1
#endif
#ifndef VTK_MESH_IUBOUND
#define VTK_MESH_IUBOUND UBOUND(VTK_MESH,2)
#endif
#ifndef VTK_MESH_JLBOUND
#define VTK_MESH_JLBOUND 1
#endif
#ifndef VTK_MESH_JUBOUND
#define VTK_MESH_JUBOUND UBOUND(VTK_MESH,3)
#endif
#if __DIM == __3D
#ifndef VTK_MESH_KLBOUND
#define VTK_MESH_KLBOUND 1
#endif
#ifndef VTK_MESH_KUBOUND
#define VTK_MESH_KUBOUND UBOUND(VTK_MESH,4)
#endif
           DO k=VTK_MESH_KLBOUND,VTK_MESH_KUBOUND
#endif
               DO j=VTK_MESH_JLBOUND,VTK_MESH_JUBOUND
                   DO i=VTK_MESH_ILBOUND,VTK_MESH_IUBOUND
#ifndef VTK_MESH_COMPONENT
#if __DIM == __2D
                       WRITE(scratch, *) VTK_MESH(i,j)
#elif __DIM == __3D
                       WRITE(scratch, *) VTK_MESH(i,j,k)
#endif
#else
#if __DIM == __2D
                       WRITE(scratch, *) VTK_MESH(VTK_MESH_COMPONENT,i,j)
#elif __DIM == __3D
                       WRITE(scratch, *) VTK_MESH(VTK_MESH_COMPONENT,i,j,k)
#endif
#endif
                       scratch = ADJUSTL(scratch)
                       WRITE(iUnit, '(A)', ADVANCE='NO') &
                       scratch(1:LEN_TRIM(scratch))
#if __DIM == __2D
                       IF (i*j.LT.(VTK_MESH_IUBOUND)*(VTK_MESH_JUBOUND)) &
#elif __DIM == __3D
                       IF (i*j*k.LT.(VTK_MESH_IUBOUND)*(VTK_MESH_JUBOUND)*(VTK_MESH_KUBOUND)) &
#endif
                       & WRITE(iUnit, '(A)', ADVANCE='NO') " "
                   ENDDO
               ENDDO
#if __DIM == __3D
           ENDDO
#endif
#elif defined VTK_SCALAR
           DO i=LBOUND(VTK_SCALAR,1),UBOUND(VTK_SCALAR,1)
              WRITE(scratch, *) VTK_SCALAR(i)
              scratch = ADJUSTL(scratch)
              WRITE(iUnit, '(A)', ADVANCE='NO') scratch(1:LEN_TRIM(scratch))
              IF (i .LT. UBOUND(VTK_SCALAR,1)) WRITE(iUnit, '(A)', ADVANCE='NO') " "
           ENDDO
#elif defined VTK_INTEGER
           DO i=LBOUND(VTK_INTEGER,1),UBOUND(VTK_INTEGER,1)
              WRITE(scratch, *) REAL(VTK_INTEGER(i),8)
              scratch = ADJUSTL(scratch)
              WRITE(iUnit, '(A)', ADVANCE='NO') scratch(1:LEN_TRIM(scratch))
              IF (i .LT. UBOUND(VTK_INTEGER,1)) WRITE(iUnit, '(A)', ADVANCE='NO') " "
           ENDDO
#else
           DO j=LBOUND(VTK_VECTOR,2),UBOUND(VTK_VECTOR,2)
              DO i=LBOUND(VTK_VECTOR,1),UBOUND(VTK_VECTOR,1)
                 WRITE(scratch, *) VTK_VECTOR(i,j)
                 scratch = ADJUSTL(scratch)
                 WRITE(iUnit, '(A)', ADVANCE='NO') scratch(1:LEN_TRIM(scratch))
                 IF (i .LT. UBOUND(VTK_VECTOR,1) .OR. &
                 &   j .LT. UBOUND(VTK_VECTOR,2))     &
                 &  WRITE(iUnit, '(A)', ADVANCE='NO') " "
              ENDDO
#ifdef APPEND_ZEROS
              IF (nd .EQ. 2) THEN
                 IF (j .LT. UBOUND(VTK_VECTOR,2)) THEN
                    WRITE(iUnit,'(A)',ADVANCE='NO') "0 "
                 ELSE
                    WRITE(iUnit,'(A)',ADVANCE='NO') " 0"
                 ENDIF
              ENDIF
#endif
           ENDDO

#endif
#endif
           WRITE(iUnit,'(/A)') "        </DataArray>"
#undef VTK_NAME
#undef VTK_TYPE
#undef VTK_NDIM
#undef VTK_OFFSET
#undef VTK_RANGE
#undef VTK_RANGE_START
#undef VTK_SCALAR
#undef VTK_INTEGER
#undef VTK_VECTOR
#undef APPEND_ZEROS
#undef VTK_MESH
#undef VTK_MESH_ILBOUND
#undef VTK_MESH_IUBOUND
#undef VTK_MESH_JLBOUND
#undef VTK_MESH_JUBOUND
#undef VTK_MESH_KLBOUND
#undef VTK_MESH_KUBOUND
#undef VTK_MESH_COMPONENT
