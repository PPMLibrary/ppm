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

FUNCTION particles_dflt_nlname(i)
    !!! Default name for an operator
    CHARACTER(LEN=ppm_char)   :: particles_dflt_nlname
    INTEGER                   :: i
    CHARACTER(LEN=ppm_char)   :: buf

    WRITE(buf,*) i
    WRITE(particles_dflt_nlname,*) 'neighlist_',ADJUSTL(TRIM(buf))
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

FUNCTION color_print(text,color)
    !!! Return a string that will appear in color if printed out to the terminal
    CHARACTER(LEN=ppm_char)   :: color_print
    CHARACTER(LEN=*)          :: text
    CHARACTER(LEN=10)         :: mycolor
    INTEGER                   :: color

    write(mycolor,'(A,I0,A)') '[',color,'m'
    color_print = achar(27)//TRIM(mycolor)//TRIM(ADJUSTL(text))//achar(27)//'[0m'
    RETURN

END FUNCTION


FUNCTION factorial_m(multi_ind,ndim)
    INTEGER                 :: factorial_m,i,ndim
    INTEGER,DIMENSION(ndim) :: multi_ind

    factorial_m = 1
    DO i=1,ndim
        factorial_m = factorial_m * factorial(multi_ind(i))
    ENDDO
END FUNCTION factorial_m

FUNCTION factorial(n)
    INTEGER    :: n,i
    INTEGER    :: factorial
    factorial = 1
    DO i=1,n
        factorial = i*factorial
    ENDDO
    RETURN
END FUNCTION factorial

FUNCTION binomial(n,k)
    INTEGER    :: n,k
    INTEGER    :: binomial
    binomial = factorial(n)/(factorial(k)*factorial(n-k))
END FUNCTION binomial


