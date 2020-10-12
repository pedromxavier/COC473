program main5
    use Func
    use Calc
    use Util
    implicit none

!   Command-line Args
    integer :: argc

!   character(len=32) :: 
!   Random seed definition
    call init_random_seed()

!   Get Command-Line Args
    argc = iargc()

    if (argc /= 0) then
        goto 101
    else
        goto 200
    end if

!   ====== Success ===================================
100 call info(':: Sucesso ::')
    goto 1
!   ====== Errors ====================================
101 call error('Este programa não aceita parâmetros.')
    goto 1
!   ====== Finish ====================================
1   stop
!   ==================================================

200 goto 100
end program main5