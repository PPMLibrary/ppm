      !--*- f90 -*--------------------------------------------------------------
      !  Module       :              ppm_module_util_functions
      !-------------------------------------------------------------------------
      ! Copyright (c) 2010 CSE Lab (ETH Zurich), MOSAIC Group (ETH Zurich), 
      !                    Center for Fluid Dynamics (DTU)
      !
      !
      ! This file is part of the Parallel Particle Mesh Library (PPM).
      !
      ! PPM is free software: you can redistribute it and/or modify
      ! it under the terms of the GNU Lesser General Public License 
      ! as published by the Free Software Foundation, either 
      ! version 3 of the License, or (at your option) any later 
      ! version.
      !
      ! PPM is distributed in the hope that it will be useful,
      ! but WITHOUT ANY WARRANTY; without even the implied warranty of
      ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
      ! GNU General Public License for more details.
      !
      ! You should have received a copy of the GNU General Public License
      ! and the GNU Lesser General Public License along with PPM. If not,
      ! see <http://www.gnu.org/licenses/>.
      !
      ! Parallel Particle Mesh Library (PPM)
      ! ETH Zurich
      ! CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Define types
      !-------------------------------------------------------------------------
#define __SINGLE_PRECISION 1
#define __DOUBLE_PRECISION 2
#define __SCALAR 3
#define __VECTOR 4

      MODULE ppm_module_util_functions
      !!! This module provides various utility functions

      USE ppm_module_data

      CONTAINS

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


      END MODULE ppm_module_util_functions
