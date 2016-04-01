      SUBROUTINE DTYPE(sop_create)(Pc,Npart,info,name)
          !!! create a set of particles
          !!! This allocates the particle positions.
          IMPLICIT NONE

          DEFINE_MK()
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(DTYPE(ppm_t_sop))                   :: Pc
          !!! Data structure containing the particles
          INTEGER,                    INTENT(IN   ) :: Npart
          !!! Number of particles
          INTEGER,                    INTENT(  OUT) :: info
          !!! Returns status, 0 upon success.
          !-------------------------------------------------------------------------
          !  Optional arguments
          !-------------------------------------------------------------------------
          CHARACTER(LEN=*), OPTIONAL, INTENT(IN   ) :: name
          !!! give a name to this Particle set
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          INTEGER :: i

          LOGICAL :: lalloc,ldealloc

          start_subroutine("sop_create")

          !Call the parent function
          CALL Pc%DTYPE(ppm_t_vbp)%create(Npart,info,name)
          or_fail("failed to initialize non-sop particle set")

          !and update the few fields that are specific to SOP
          Pc%adapt_wp   => NULL()
          Pc%D          => NULL()
          Pc%Dtilde     => NULL()
          ! Particles do not represent a level-set function
          Pc%level_set  = .FALSE.
          Pc%level      => NULL()
          !        Pc%level_old_id = 0
          Pc%level_grad => NULL()
          !        Pc%level_grad_old_id = 0
          ! Particles are by default isotropic
          Pc%anisotropic= .FALSE.

          end_subroutine()
      END SUBROUTINE DTYPE(sop_create)

