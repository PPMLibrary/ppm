test_suite ppm_module_container_typedef

use ppm_module_interfaces
use ppm_module_field_typedef
use ppm_module_mesh_typedef

#ifdef __MPI
    INCLUDE "mpif.h"
#endif

integer                          :: i,j,k
integer                          :: info

!---------------- init -----------------------

    init

    end init

!----------------------------------------------

!---------------- finalzie --------------------


    finalize

    end finalize

!----------------------------------------------

!------------- setup --------------------------

    setup

    end setup
!----------------------------------------------
        

!--------------- teardown ---------------------
    teardown

    end teardown
!----------------------------------------------

    test collections_basic_features
       type(ppm_c_equi_mesh) :: c_M
       class(ppm_t_equi_mesh_),pointer :: M => NULL()
       integer                      :: count_el,total_el=1000

       ppm_debug = 0

       do i=1,total_el
           allocate(ppm_t_equi_mesh::M,STAT=info)
           Assert_Equal(info,0)
           call c_M%push(M,info)
           Assert_Equal(info,0)
           M => NULL()
       enddo
       call c_M%destroy(info)
       Assert_Equal(info,0)

       do i=1,total_el
           allocate(ppm_t_equi_mesh::M,STAT=info)
           Assert_Equal(info,0)
           call c_M%push(M,info)
           Assert_Equal(info,0)
           M => NULL()
       enddo

       count_el = 0
       M => c_M%begin()
       do while (associated(M))
           count_el = count_el + 1
           M => c_M%next()
       enddo
       Assert_Equal(count_el,total_el)

       count_el = 0
       M => c_M%begin()
       do while (associated(M))
           count_el = count_el + 1
           if (mod(count_el,3).eq.0) then
               call c_M%remove(info)
               Assert_Equal(info,0)
           endif
           M => c_M%next()
       enddo

       count_el = 0
       M => c_M%last()
       do while (associated(M))
           count_el = count_el + 1
           M => c_M%prev()
       enddo
       Assert_Equal(count_el,total_el-total_el/3)

       M => c_M%begin()
       do while (associated(M))
           call c_M%remove(info)
           Assert_Equal(info,0)
           M => c_M%next()
       enddo
       Assert_Equal(c_M%nb,0)
    end test


    test nested_collections
       type(ppm_c_equi_mesh) :: c_M
       class(ppm_t_equi_mesh_),pointer :: M => NULL()
       class(ppm_t_subpatch_),  pointer :: p => NULL()
       integer                      :: count_el,total_el=100

       ! build a collection of meshes, each of which has
       ! a collection of subpatches
       do i=1,total_el
           allocate(ppm_t_equi_mesh::M,STAT=info)
           Assert_Equal(info,0)
           allocate(ppm_c_subpatch::M%subpatch)
           M%ID = i
           do j=1,total_el
               allocate(ppm_t_subpatch::p,STAT=info)
               Assert_Equal(info,0)
               p%meshID = j
               call M%subpatch%push(p,info)
               Assert_Equal(info,0)
           enddo
           Assert_True(associated(M%subpatch))
           call c_M%push(M,info)
           Assert_Equal(info,0)
       enddo

       !check that the subpatches are still allocated for each mesh
       i = 1
       M => c_M%begin()
       do while (associated(M))
           Assert_Equal(i,M%ID)
           j = 1
           Assert_True(associated(M%subpatch))
           p => M%subpatch%begin()
           do while (associated(p))
               Assert_Equal(j,p%meshID)
               p => M%subpatch%next()
               j = j+1
           enddo
           Assert_Equal(j-1,total_el)
           M => c_M%next()
           i = i +1
       enddo

       call c_M%destroy(info)
       Assert_Equal(info,0)

    end test




end test_suite
