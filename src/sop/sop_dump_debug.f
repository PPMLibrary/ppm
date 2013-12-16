
      SUBROUTINE DTYPE(sop_dump_2d)(field,n1,n2,iunit,info,dir)

          IMPLICIT NONE

          DEFINE_MK()
          ! arguments
          REAL(MK),  DIMENSION(:,:),    INTENT(IN   )   :: field
          INTEGER,                      INTENT(IN   )   :: n1
          INTEGER,                      INTENT(IN   )   :: n2
          INTEGER,                      INTENT(IN   )   :: iunit
          INTEGER,                      INTENT(  OUT)   :: info

          CHARACTER(LEN=*),OPTIONAL,    INTENT(IN)      :: dir

          ! local variable
          INTEGER                                       :: i,j
          CHARACTER(LEN=256)                            :: myformat,filename
          CHARACTER(LEN=256)                            :: filenum
          CHARACTER(LEN=256)                            :: caller='sop_dump_debug'

          !!--------------------------------------------------------------------
          !! Initialize
          !!--------------------------------------------------------------------
          info = 0

          !!--------------------------------------------------------------------
          !! Write output file
          !!--------------------------------------------------------------------
          WRITE(myformat,'(A,I1,A)') '(',n1,'E25.15E3,2X)'

          CLOSE(iunit)


          WRITE(filenum,*) iunit
          IF (PRESENT(dir)) THEN
              WRITE(filename,'(A,A,A)') TRIM(dir),'/fort.',TRIM(ADJUSTL(filenum))
          ELSE
              WRITE(filename,'(A,A)') './fort.',TRIM(ADJUSTL(filenum))
          ENDIF
          OPEN(UNIT=iunit,FILE=filename,STATUS='REPLACE')

          DO j=1,n2
              WRITE(iunit,myformat) (field(i,j), i=1,n1)
          ENDDO
          CLOSE(iunit)

      END SUBROUTINE DTYPE(sop_dump_2d)


      SUBROUTINE DTYPE(sop_dump_1d)(field,n1,iunit,info,dir)


          IMPLICIT NONE

          DEFINE_MK()
          ! arguments
          REAL(MK),  DIMENSION(:),      INTENT(IN   )   :: field
          INTEGER,                      INTENT(IN   )   :: n1
          INTEGER,                      INTENT(IN   )   :: iunit
          INTEGER,                      INTENT(  OUT)   :: info

          CHARACTER(LEN=*),OPTIONAL,    INTENT(IN)      :: dir

          ! local variable
          INTEGER                                       :: i
          CHARACTER(LEN=256)                            :: myformat,filename
          CHARACTER(LEN=256)                            :: filenum
          CHARACTER(LEN=256)                            :: caller='sop_dump_debug'

          !!--------------------------------------------------------------------
          !! Initialize
          !!--------------------------------------------------------------------
          info = 0

          !!--------------------------------------------------------------------
          !! Write output file
          !!--------------------------------------------------------------------
          WRITE(myformat,'(A)') '(E25.15E3)'

          CLOSE(iunit)
          WRITE(filenum,*) iunit
          IF (PRESENT(dir)) THEN
              WRITE(filename,'(A,A,A)') TRIM(dir),'/fort.',TRIM(ADJUSTL(filenum))
          ELSE
              WRITE(filename,'(A,A)') './fort.',TRIM(ADJUSTL(filenum))
          ENDIF
          OPEN(UNIT=iunit,FILE=filename,STATUS='REPLACE')
          DO i=1,n1
              WRITE(iunit,myformat) field(i)
          ENDDO
          CLOSE(iunit)

      END SUBROUTINE DTYPE(sop_dump_1d)

#if __KIND == __SINGLE_PRECISION
      SUBROUTINE sop_dump_1di(field,n1,iunit,info,dir)

          IMPLICIT NONE

          ! arguments
          INTEGER,  DIMENSION(:),       INTENT(IN   )   :: field
          INTEGER,                      INTENT(IN   )   :: n1
          INTEGER,                      INTENT(IN   )   :: iunit
          INTEGER,                      INTENT(  OUT)   :: info

          CHARACTER(LEN=*),OPTIONAL,    INTENT(IN)      :: dir

          ! local variable
          INTEGER                                       :: i
          CHARACTER(LEN=256)                            :: myformat,filename
          CHARACTER(LEN=256)                            :: filenum
          CHARACTER(LEN=256)                            :: caller='sop_dump_debug'

          !!--------------------------------------------------------------------
          !! Initialize
          !!--------------------------------------------------------------------
          info = 0

          !!--------------------------------------------------------------------
          !! Write output file
          !!--------------------------------------------------------------------
          WRITE(myformat,'(A)') '(I7)'

          CLOSE(iunit)
          WRITE(filenum,*) iunit
          IF (PRESENT(dir)) THEN
              WRITE(filename,'(A,A,A)') TRIM(dir),'/fort.',TRIM(ADJUSTL(filenum))
          ELSE
              WRITE(filename,'(A,A)') './fort.',TRIM(ADJUSTL(filenum))
          ENDIF
          OPEN(UNIT=iunit,FILE=filename,STATUS='REPLACE')
          DO i=1,n1
              WRITE(iunit,myformat) field(i)
          ENDDO
          CLOSE(iunit)

      END SUBROUTINE sop_dump_1di

      SUBROUTINE sop_dump_2di(field,n1,n2,iunit,info,dir)


          IMPLICIT NONE

          ! arguments
          INTEGER,  DIMENSION(:,:),     INTENT(IN   )   :: field
          INTEGER,                      INTENT(IN   )   :: n1
          INTEGER,                      INTENT(IN   )   :: n2
          INTEGER,                      INTENT(IN   )   :: iunit
          INTEGER,                      INTENT(  OUT)   :: info

          CHARACTER(LEN=*),OPTIONAL,    INTENT(IN)      :: dir

          ! local variable
          INTEGER                                       :: i,j
          CHARACTER(LEN=256)                            :: myformat,filename
          CHARACTER(LEN=256)                            :: filenum
          CHARACTER(LEN=256)                            :: caller='sop_dump_debug'

          !!--------------------------------------------------------------------
          !! Initialize
          !!--------------------------------------------------------------------
          info = 0

          !!--------------------------------------------------------------------
          !! Write output file
          !!--------------------------------------------------------------------
          IF (n1 .LT. 100) THEN
              WRITE(myformat,'(A,I2,A)') '(',n1,'I6,2X)'
          ELSE IF (n1 .LT. 1000) THEN
              WRITE(myformat,'(A,I3,A)') '(',n1,'I6,2X)'
          ELSE IF (n1 .LT. 10000) THEN
              WRITE(myformat,'(A,I4,A)') '(',n1,'I6,2X)'
          ENDIF

          CLOSE(iunit)
          WRITE(filenum,*) iunit
          IF (PRESENT(dir)) THEN
              WRITE(filename,'(A,A,A)') TRIM(dir),'/fort.',TRIM(ADJUSTL(filenum))
          ELSE
              WRITE(filename,'(A,A)') './fort.',TRIM(ADJUSTL(filenum))
          ENDIF
          OPEN(UNIT=iunit,FILE=filename,STATUS='REPLACE')
          DO j=1,n2
              WRITE(iunit,myformat) (field(i,j), i=1,n1)
          ENDDO
          CLOSE(iunit)

      END SUBROUTINE sop_dump_2di
#endif
