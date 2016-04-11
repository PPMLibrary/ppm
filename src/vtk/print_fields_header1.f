#ifndef VTK_PARALLEL
         WRITE(iUnit,'(A)', ADVANCE='NO') "    <Piece"
#ifdef VTK_EXTENT
         WRITE(iUnit, '(A)', ADVANCE='NO') " Extent='"
         DO i=LBOUND(VTK_EXTENT,1),UBOUND(VTK_EXTENT,1)
            WRITE(scratch, *) VTK_EXTENT(i)
            scratch = ADJUSTL(scratch)
            WRITE(iUnit, '(A)', ADVANCE='NO') scratch(1:LEN_TRIM(scratch))
            IF (i .LT. UBOUND(VTK_EXTENT,1)) WRITE(iUnit, '(A)', ADVANCE='NO') " "
         ENDDO
         WRITE(iUnit, '(A)', ADVANCE='NO') "'"
#else
#ifdef VTK_NPOINTS
         WRITE(scratch, *) VTK_NPOINTS
         scratch = ADJUSTL(scratch)
         WRITE(iUnit, '(3A)', ADVANCE='NO') " NumberOfPoints='", &
         & scratch(1:LEN_TRIM(scratch)), "'"
#endif
#ifdef VTK_NVERTS
         WRITE(scratch, *) VTK_NVERTS
         scratch = ADJUSTL(scratch)
         WRITE(iUnit, '(3A)', ADVANCE='NO') " NumberOfVerts='", &
         & scratch(1:LEN_TRIM(scratch)), "'"
#endif
#ifdef VTK_NLINES
         WRITE(scratch, *) VTK_NLINES
         scratch = ADJUSTL(scratch)
         WRITE(iUnit, '(3A)', ADVANCE='NO') " NumberOfLines='", &
         & scratch(1:LEN_TRIM(scratch)), "'"
#endif
#ifdef VTK_NSTRIPS
         WRITE(scratch, *) VTK_NSTRIPS
         scratch = ADJUSTL(scratch)
         WRITE(iUnit, '(3A)', ADVANCE='NO') " NumberOfStrips='", &
         & scratch(1:LEN_TRIM(scratch)), "'"
#endif
#ifdef VTK_NPOLYS
         WRITE(scratch, *) VTK_NPOLYS
         scratch = ADJUSTL(scratch)
         WRITE(iUnit, '(3A)', ADVANCE='NO') " NumberOfPolys='", &
         & scratch(1:LEN_TRIM(scratch)), "'"
#endif
#endif
         WRITE(iUnit, '(A)') ">"
#endif
#undef VTK_EXTENT
#undef VTK_NPOINTS
#undef VTK_NVERTS
#undef VTK_NLINES
#undef VTK_NSTRIPS
#undef VTK_NPOLYS
