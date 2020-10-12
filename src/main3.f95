program main3
    use Matrix

    implicit none

!   Command-line Args
    integer :: argc

    character(len=32) :: vecx_fname, vecy_fname

!   Vectors
    integer :: m, n
    double precision, allocatable :: x(:)
    double precision, allocatable :: y(:)
    double precision, allocatable :: s(:)

!   Get Command-Line Args
    vecx_fname = 'vecx.txt'
    vecy_fname = 'vecy.txt'

    argc = iargc()

    if (argc == 1) then
        call getarg(1, vecx_fname)
    elseif(argc == 2) then
        call getarg(1, vecx_fname)
        call getarg(2, vecy_fname)
    elseif (argc > 2) then
        goto 92
    end if

    call read_vector(vecx_fname, x, m)
    call read_vector(vecy_fname, y, n)

    if (m /= n) then
        goto 90
    end if

!   Print Matrix
!   Print Matrix Name
    write(*, *) 'x:'
    call print_vector(x, m)

!   Print Matrix
!   Print Matrix Name
    write(*, *) 'y:'
    call print_vector(y, n)

!   Least Squares
    allocate(s(2))
    if (.not. least_squares(x, y, s, n)) then
        goto 80
    end if 

    write(*, *) 'a:'
    write(*, *) s(2)
    write(*, *) 'b:'
    write(*, *) s(1)

    goto 85

80  call error('Falho no método de mínimos quadrados.')
    goto 100
85  call info('Sucesso!')
    goto 100

90  call error('Os vetores x e y devem ter o mesmo tamnho.')
    goto 100
92  call error('Parâmetros em excesso: O programa espera apenas o nome do arquivo da matriz.')
    goto 100

100 deallocate(x)
    deallocate(y)
    deallocate(s)
    stop

end program main3