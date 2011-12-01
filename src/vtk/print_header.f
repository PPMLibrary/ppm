#ifndef VTK_VERSION
#define VTK_VERSION "0.1"
#endif
#ifndef VTK_BYTE_ORDER
#define VTK_BYTE_ORDER "LittleEndian"
#endif
           WRITE(iUnit,'(A)')  "<?xml version='1.0'?>"
           WRITE(iUnit,'(7A)') "<VTKFile type='", VTK_FILE_TYPE, &
                               "' version='", VTK_VERSION,  &
                               "' byte_order='", VTK_BYTE_ORDER, "'>"
           WRITE(iUnit,'(2A)',advance='no') "  <", VTK_FILE_TYPE
#ifdef VTK_WHOLE_EXTENT
           WRITE(iUnit, '(A)', advance='no') " WholeExtent='"
           DO i=LBOUND(VTK_WHOLE_EXTENT,1),UBOUND(VTK_WHOLE_EXTENT,1)
              WRITE(scratch, *) VTK_WHOLE_EXTENT(i)
              scratch = ADJUSTL(scratch)
              WRITE(iUnit, '(A)', advance='no') &
                   scratch(1:LEN_TRIM(scratch))
              IF (i .LT. UBOUND(VTK_WHOLE_EXTENT,1)) &
                   WRITE(iUnit, '(A)', advance='no') " "
           END DO
           WRITE(iUnit, '(A)', advance='no') "'"
#endif
#ifdef VTK_GHOSTLEVEL
           WRITE(scratch, *) VTK_GHOSTLEVEL
           scratch = ADJUSTL(scratch)
           WRITE(iUnit, '(3A)', advance='no') " GhostLevel='", &
                scratch(1:LEN_TRIM(scratch)), "'"
#endif
#ifdef VTK_ORIGIN
           WRITE(iUnit, '(A)', advance='no') " Origin='"
           DO i=LBOUND(VTK_ORIGIN,1),UBOUND(VTK_ORIGIN,1)
              WRITE(scratch, *) VTK_ORIGIN(i)
              scratch = ADJUSTL(scratch)
              WRITE(iUnit, '(A)', advance='no') &
                   scratch(1:LEN_TRIM(scratch))
              IF (i .LT. UBOUND(VTK_ORIGIN,1)) &
                   WRITE(iUnit, '(A)', advance='no') " "
           END DO
           WRITE(iUnit, '(A)', advance='no') "'"
#endif
#ifdef VTK_SPACING
           WRITE(iUnit, '(A)', advance='no') " Spacing='"
           DO i=LBOUND(VTK_SPACING,1),UBOUND(VTK_SPACING,1)
              WRITE(scratch, *) VTK_SPACING(i)
              scratch = ADJUSTL(scratch)
              WRITE(iUnit, '(A)', advance='no') &
                   scratch(1:LEN_TRIM(scratch))
              IF (i .LT. UBOUND(VTK_SPACING,1)) &
                   WRITE(iUnit, '(A)', advance='no') " "
           END DO
           WRITE(iUnit, '(A)', advance='no') "'"
#endif
           WRITE(iUnit,'(A)') ">"
#ifndef VTK_PARALLEL
           WRITE(iUnit,'(A)', advance='no')     "    <Piece"
#ifdef VTK_EXTENT
           WRITE(iUnit, '(A)', advance='no') " Extent='"
           DO i=LBOUND(VTK_EXTENT,1),UBOUND(VTK_EXTENT,1)
              WRITE(scratch, *) VTK_EXTENT(i)
              scratch = ADJUSTL(scratch)
              WRITE(iUnit, '(A)', advance='no') &
                   scratch(1:LEN_TRIM(scratch))
              IF (i .LT. UBOUND(VTK_EXTENT,1)) &
                   WRITE(iUnit, '(A)', advance='no') " "
           END DO
           WRITE(iUnit, '(A)', advance='no') "'"
#else
#ifdef VTK_NPOINTS
           WRITE(scratch, *) VTK_NPOINTS
           scratch = ADJUSTL(scratch)
           WRITE(iUnit, '(3A)', advance='no') " NumberOfPoints='", &
                scratch(1:LEN_TRIM(scratch)), "'"
#endif
#ifdef VTK_NVERTS
           WRITE(scratch, *) VTK_NVERTS
           scratch = ADJUSTL(scratch)
           WRITE(iUnit, '(3A)', advance='no') " NumberOfVerts='", &
                scratch(1:LEN_TRIM(scratch)), "'"
#endif
#ifdef VTK_NLINES
           WRITE(scratch, *) VTK_NLINES
           scratch = ADJUSTL(scratch)
           WRITE(iUnit, '(3A)', advance='no') " NumberOfLines='", &
                scratch(1:LEN_TRIM(scratch)), "'"
#endif
#ifdef VTK_NSTRIPS
           WRITE(scratch, *) VTK_NSTRIPS
           scratch = ADJUSTL(scratch)
           WRITE(iUnit, '(3A)', advance='no') " NumberOfStrips='", &
                scratch(1:LEN_TRIM(scratch)), "'"
#endif
#ifdef VTK_NPOLYS
           WRITE(scratch, *) VTK_NPOLYS
           scratch = ADJUSTL(scratch)
           WRITE(iUnit, '(3A)', advance='no') " NumberOfPolys='", &
                scratch(1:LEN_TRIM(scratch)), "'"
#endif
#endif
           WRITE(iUnit, '(A)') ">"
#endif
#undef VTK_VERSION
#undef VTK_WHOLE_EXTENT
#undef VTK_ORIGIN
#undef VTK_GHOSTLEVEL
#undef VTK_SPACING
#undef VTK_EXTENT
#undef VTK_NPOINTS
#undef VTK_NVERTS
#undef VTK_NLINES
#undef VTK_NSTRIPS
#undef VTK_NPOLYS
