program main3
    use Func
    use Calc

    implicit none

!   Command-line Args
    integer :: argc

!   character(len=32) :: 
    double precision :: x, y, a, b, x0

!   Get Command-Line Args
    argc = iargc()

    if (argc /= 0) then
        goto 90
    end if

    a = -1.0D0
    b =  1.0D0
    x0 = 2.0D0 * RAN(0) - 1.0D0

    call info("Bisseção:")

    x = bissection(f2, a, b)
    y = f2(x)

    write(*, *) "x:"
    write(*, *) x

    write(*, *) "y:"
    write(*, *) y

    call info("Newton:")

    x = newton(f2, df2, x0)
    y = f2(x)

    write(*, *) "x:"
    write(*, *) x

    write(*, *) "y:"
    write(*, *) y

    call info("Secante:")

    x = secant(f2, x0)
    y = f2(x)

    write(*, *) "x:"
    write(*, *) x

    write(*, *) "y:"
    write(*, *) y

    goto 100

90  call error('Este programa não aceita parâmetros.')
    goto 101

100 call info('Sucesso!')
101 stop

end program main3