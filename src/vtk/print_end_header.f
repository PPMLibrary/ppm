#ifndef VTK_PARALLEL
         WRITE(iUnit,'(A)')  "    </Piece>"
#endif
         WRITE(iUnit,'(3A)') "  </", VTK_FILE_TYPE(1:LEN_TRIM(VTK_FILE_TYPE)), ">"
         WRITE(iUnit,'(A)')  "</VTKFile>"
#undef VTK_FILE_TYPE
#undef VTK_PARALLEL
