FUNCTION particles_dflt_partname(i)
    !!! Default name for a set of particles
    CHARACTER(LEN=ppm_char)   :: particles_dflt_partname
    INTEGER,OPTIONAL          :: i
    CHARACTER(LEN=ppm_char)   :: buf

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
    CHARACTER(LEN=ppm_char)   :: particles_dflt_pptname
    INTEGER                   :: i,ndim
    CHARACTER(LEN=ppm_char)   :: buf

    WRITE(buf,*) i
    IF (ndim .EQ. 1) THEN
        WRITE(particles_dflt_pptname,*) 'property_s',ADJUSTL(TRIM(buf))
    ELSE
        WRITE(particles_dflt_pptname,*) 'property_v',ADJUSTL(TRIM(buf))
    ENDIF
    RETURN
END FUNCTION

FUNCTION particles_dflt_opname(i)
    !!! Default name for an operator
    CHARACTER(LEN=ppm_char)   :: particles_dflt_opname
    INTEGER                   :: i
    CHARACTER(LEN=ppm_char)   :: buf

    WRITE(buf,*) i
    WRITE(particles_dflt_opname,*) 'operator_',ADJUSTL(TRIM(buf))
    RETURN
END FUNCTION


