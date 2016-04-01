      FUNCTION particles_dflt_partname(i)
          !!! Default name for a set of particles
          IMPLICIT NONE

          CHARACTER(LEN=ppm_char)          :: particles_dflt_partname
          INTEGER, OPTIONAL, INTENT(IN   ) :: i

          CHARACTER(LEN=ppm_char) :: buf

          IF (PRESENT(i)) THEN
             WRITE(buf,*) i
             WRITE(particles_dflt_partname,*) 'Particles_',ADJUSTL(TRIM(buf))
          ELSE
             WRITE(particles_dflt_partname,*) 'Particles'
          ENDIF
          RETURN
      END FUNCTION

      FUNCTION particles_dflt_pptname(i,ndim)
          !!! Default name for a scalar or vector property
          IMPLICIT NONE

          CHARACTER(LEN=ppm_char) :: particles_dflt_pptname
          INTEGER,  INTENT(IN   ) :: i,ndim

          CHARACTER(LEN=ppm_char) :: buf

          WRITE(buf,*) i
          IF (ndim .EQ. 1) THEN
             WRITE(particles_dflt_pptname,*) 'property_s',ADJUSTL(TRIM(buf))
          ELSE
             WRITE(particles_dflt_pptname,*) 'property_v',ADJUSTL(TRIM(buf))
          ENDIF
          RETURN
      END FUNCTION

      FUNCTION particles_dflt_nlname(i)
          !!! Default name for a neighlist
          IMPLICIT NONE

          CHARACTER(LEN=ppm_char) :: particles_dflt_nlname
          INTEGER,  INTENT(IN   ) :: i

          CHARACTER(LEN=ppm_char)   :: buf

          WRITE(buf,*) i
          WRITE(particles_dflt_nlname,*) 'neighlist_',ADJUSTL(TRIM(buf))
          RETURN
      END FUNCTION

      FUNCTION particles_dflt_opname(i)
          !!! Default name for an operator
          IMPLICIT NONE

          CHARACTER(LEN=ppm_char) :: particles_dflt_opname
          INTEGER,  INTENT(IN   ) :: i

          CHARACTER(LEN=ppm_char) :: buf

          WRITE(buf,*) i
          WRITE(particles_dflt_opname,*) 'operator_',ADJUSTL(TRIM(buf))
          RETURN
      END FUNCTION


