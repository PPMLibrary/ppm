test_suite ppm_module_container_typedef

  USE ppm_module_interfaces
  USE ppm_module_field_typedef
  USE ppm_module_mesh_typedef

  INTEGER, PARAMETER :: debug = 0
  INTEGER, PARAMETER :: MK = KIND(1.0D0) !KIND(1.0E0)
  INTEGER, PARAMETER :: ndim=2
  INTEGER            :: tolexp
  INTEGER            :: info
  INTEGER            :: i,j,k

  !---------------- init -----------------------
  init
    USE ppm_module_init

    tolexp = INT(LOG10(EPSILON(1.0_MK)))
    CALL ppm_init(ndim,MK,tolexp,0,debug,info,99)
  end init
  !----------------------------------------------

  !------------- setup --------------------------
  setup
  end setup
  !----------------------------------------------

  !--------------- tearDOwn ---------------------
  tearDOwn
  end tearDOwn
  !----------------------------------------------

  !---------------- finalzie --------------------
  finalize
    USE ppm_module_finalize

    CALL ppm_finalize(info)
  end finalize
  !----------------------------------------------

  test collections_basic_features
    IMPLICIT NONE
    TYPE(ppm_c_equi_mesh)            :: c_M
    CLASS(ppm_t_equi_mesh_), POINTER :: M => NULL()
    INTEGER                          :: count_el,total_el=1000

    ppm_debug = 0

    DO i=1,total_el
       ALLOCATE(ppm_t_equi_mesh::M,STAT=info)
       Assert_Equal(info,0)

       CALL c_M%push(M,info)
       Assert_Equal(info,0)
    ENDDO

    CALL c_M%destroy(info)
    Assert_Equal(info,0)

    DO i=1,total_el
       ALLOCATE(ppm_t_equi_mesh::M,STAT=info)
       Assert_Equal(info,0)

       M%topoid = i
       CALL c_M%push(M,info)
       Assert_Equal(info,0)
    ENDDO

    count_el = 0
    M => c_M%begin()
    DO WHILE (ASSOCIATED(M))
       count_el = count_el + 1
       M => c_M%next()
    ENDDO
    Assert_Equal(count_el,total_el)

    ! check that the accessor (at) works
    M => c_M%at(17)
    Assert_Equal(M%topoid,17)
    M => c_M%at(total_el+1)
    Assert_False(ASSOCIATED(M))
    M => c_M%at(1)
    Assert_Equal(M%topoid,1)
    M => c_M%at(0)
    Assert_False(ASSOCIATED(M))

    count_el = 0
    M => c_M%begin()
    DO WHILE (ASSOCIATED(M))
        count_el = count_el + 1
        IF (mod(count_el,3).EQ.0) then
         CALL c_M%remove(info)
         Assert_Equal(info,0)
        ENDIF
        M => c_M%next()
    ENDDO

    count_el = 0
    M => c_M%last()
    DO WHILE (ASSOCIATED(M))
        count_el = count_el + 1
        M => c_M%prev()
    ENDDO
    Assert_Equal(count_el,total_el-total_el/3)

    M => c_M%begin()
    DO WHILE (ASSOCIATED(M))
       CALL c_M%remove(info)
       Assert_Equal(info,0)
       M => c_M%next()
    ENDDO
    Assert_Equal(c_M%nb,0)
  end test

  test nested_collections
    IMPLICIT NONE
    TYPE(ppm_c_equi_mesh)            :: c_M
    CLASS(ppm_t_equi_mesh_), POINTER :: M => NULL()
    CLASS(ppm_t_subpatch_),  POINTER :: p => NULL()

    INTEGER :: count_el,total_el=100

    ! build a collection of meshes, each of which has
    ! a collection of subpatches
    DO i=1,total_el
       ALLOCATE(ppm_t_equi_mesh::M,STAT=info)
       Assert_Equal(info,0)

       ALLOCATE(ppm_c_subpatch::M%subpatch)
       M%ID = i

       DO j=1,total_el
          ALLOCATE(ppm_t_subpatch::p,STAT=info)
          Assert_Equal(info,0)
          p%meshID = j
          CALL M%subpatch%push(p,info)
          Assert_Equal(info,0)
       ENDDO
       Assert_True(ASSOCIATED(M%subpatch))
       CALL c_M%push(M,info)
       Assert_Equal(info,0)
    ENDDO

    !check that the subpatches are still ALLOCATEd for each mesh
    i = 1
    M => c_M%begin()
    DO WHILE (ASSOCIATED(M))
       Assert_Equal(i,M%ID)
       j = 1
       Assert_True(ASSOCIATED(M%subpatch))
       p => M%subpatch%begin()
       DO WHILE (ASSOCIATED(p))
          Assert_Equal(j,p%meshID)
          p => M%subpatch%next()
          j = j+1
       ENDDO
       Assert_Equal(j-1,total_el)
       M => c_M%next()
       i = i +1
    ENDDO

    CALL c_M%destroy(info)
    Assert_Equal(info,0)

  end test

end test_suite
