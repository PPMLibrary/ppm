
#if   __KIND == __SINGLE_PRECISION
subroutine ppm_dbg_print_s(topoid,xp,np,ghostlayer,step,colortag,info)
#elif __KIND == __DOUBLE_PRECISION
subroutine ppm_dbg_print_d(topoid,xp,np,ghostlayer,step,colortag,info)
#endif

    use ppm_module_data
    use ppm_module_topo
    USE ppm_module_substart
    USE ppm_module_substop

    implicit none
#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
    
    ! arguments
    integer,            intent(in)    :: topoid
    real(mk), dimension(:,:), pointer :: xp
    integer,            intent(in)    :: np
    real(mk),           intent(in)    :: ghostlayer
    integer,            intent(in)    :: step
    integer,            intent(in)    :: colortag
    integer,            intent(out)   :: info
    ! local vars
    character(64)                     :: sfmt,pfmt
    character(64)                     :: sfname,pfname
    type(ppm_t_topo), pointer         :: topo
    integer                           :: i
    integer                           :: iunit
    real(mk)                          :: t0
      
    CALL substart('ppm_dbg_print',t0,info)

    iunit = 2342

    write(sfmt,'(A,I1,A)') '(',ppm_dim*2,'F12.8,I4)'
    write(pfmt,'(A,I1,A)') '(',ppm_dim,'E18.9,I4)'
    write(sfname,'(A,I5.5,A)') 'ppm_dbg_',step,'.sub'
    write(pfname,'(A,I5.5,A)') 'ppm_dbg_',step,'.dat'

    call ppm_topo_get(topoid,topo,info)
    open(iunit,file=sfname)
   
    write(iunit,'(I)') ppm_dim
    write(iunit,'(F12.8)') ghostlayer
    do i=1,topo%nsubs
        write(iunit,sfmt) topo%min_subd(:,i),topo%max_subd(:,i),topo%sub2proc(i)
    enddo
    close(iunit)
    open(iunit,file=pfname,access='append')
    do i=1,np
        write(iunit,pfmt) xp(:,i),colortag
    enddo
    close(iunit)

      CALL substop('ppm_dbg_print',t0,info)

#if   __KIND == __SINGLE_PRECISION
end subroutine ppm_dbg_print_s
#elif __KIND == __DOUBLE_PRECISION
end subroutine ppm_dbg_print_d
#endif
