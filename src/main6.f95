program main6
    use Func
    use Calc
    use Util
    use Plot
    implicit none

    integer :: i, j

    integer :: n = 10
    double precision :: a, b, y0, dt

    double precision, dimension(:), allocatable :: t, y
    double precision, dimension(:, :), allocatable :: P

    type(StringArray), dimension(:), allocatable :: legend

    character(len=32) :: n_chr

!   Command-line Args
    integer :: argc

!    ENABLE_DEBUG = .TRUE.

!   Random seed definition
    call init_random_seed()

!   Get Command-Line Args
    argc = iargc()

    if (argc == 0) then
        goto 099
    else if (argc == 1) then
        call getarg(1, n_chr)
        read (unit=n_chr, fmt=*) n
        goto 099
    else
        goto 101
    end if

!   ====== Begin =====================================    
099 call info("Nº de pontos de integração: "//STR(n))
    goto 200

!   ====== Success ===================================
100 call info(':: Sucesso ::')
    goto 1
!   ====== Errors ====================================
101 call error('Este programa não aceita mais do que um parâmetro (número de pontos de integração).')
    goto 1
!   ====== Finish ====================================
1   stop
!   ==================================================

!   ===============================
201 i = i + 1
    go to(210, 220, 230, 299), i
!   ===============================
301 i = i + 1
    go to(310, 400), i
!   ===============================
400 call info(ENDL//"3) "//F78_NAME)
    i = 0
    goto 401

401 i = i + 1
    go to(410, 500), i
!   ===============================   
500 goto 100

200 call info(ENDL//"1) "//F13_NAME)
    i = 0
    call linspace(F13_A, F13_B, F13_DT, n, t)
    allocate(P(n, 4))
    allocate(y(n))
    allocate(legend(4))
    call show('dt', dt)
    goto 201

210 call info(":: Método de Euler ::")
    y(:) = euler(df13, y0, t, n)
    P(:, 1) = y(:)
    legend(1)%str = "Euler"
    goto 201

220 call info(":: Método de Runge-Kutta de 2ª ordem ::")
    y(:) = runge_kutta2(df13, y0, t, n)
    P(:, 2) = y(:)
    legend(2)%str = "Runge-Kutta II"
    goto 201

230 call info(":: Método de Runge-Kutta de 4ª ordem ::")
    y(:) = runge_kutta4(df13, y0, t, n)
    P(:, 3) = y(:)
    legend(3)%str = "Runge-Kutta IV"
    goto 201

299 call info(":: Plot 1 ::")
    P(:, 4) = (/ (f13(t(j)), j=1, n) /)
    legend(4)%str = "y(t)"
    call xy_multiplot(t, P, n, 4, fname='ex1', title="y'(t) = -2 t (y)²; y(0) = 1", xlabel='t', ylabel='y(t)', legend=legend)
    deallocate(P)
    deallocate(t)
    deallocate(y)
    deallocate(legend)
    goto 300

300 call info(ENDL//"2) "//F14_NAME)
    i = 0
    goto 301

310 goto 301

399 call info(":: Plot 2 ::")

    goto 400

410 goto 401
end program main6