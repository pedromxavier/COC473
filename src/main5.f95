program main5
    use Util
    use Func
    use Matrix
    use Calc
    use Plotlib
    implicit none

    double precision :: XMIN, XMAX, YMIN, YMAX

!   Command-line Args
    integer :: argc

!    ENABLE_DEBUG = .TRUE.

!   Random seed definition
    call init_random_seed()

!   Get Command-Line Args
    argc = iargc()

    if (argc == 0) then
        goto 100
    else
        goto 11
    end if

!   ====== Success ===================================
10  call info(':: Sucesso ::')
    goto 1
!   ====== Errors ====================================
11  call error('Este programa não aceita parâmetross.')
    goto 1
!   ====== Finish ====================================
1   stop
!   ==================================================

100 goto 200

200 call Q2
    goto 300

300 call Q3
    goto 400

400 call Q4
    goto 500

500 call Q5
    goto 600

600 call Q6
    goto 700

700 call Q7
    goto 800

800 call warn(ENDL//":: Complmento ::"//ENDL)
    call QE1; call QE2; call QE3;
    goto 10

!   ===============================

    contains

    subroutine Q2
        implicit none
        integer :: n = 10
        double precision :: a, b, s

        call info("2)"//ENDL//F6_NAME)
        a = 0.0D0
        b = 1.0D0
        call info("[a, b] = ["//DSTR(a)//", "//DSTR(b)//"]")
        call info(":: Integração Polinomial ::")
        s = num_int(f6, a, b, n, kind="polynomial")
        call blue("I1 = ∫f(x) dx ≈ "//DSTR(s))
        call info(":: Quadratura de Gauss-Legendre ::")
        s = num_int(f6, a, b, n, kind="gauss-legendre")
        call blue("I1 = ∫f(x) dx ≈ "//DSTR(s))
        call info(":: Método de Romberg ::")
        s = num_int(f6, a, b, n, kind="romberg")
        call blue("I1 = ∫f(x) dx ≈ "//DSTR(s))

        a = 0.0D0
        b = 5.0D0
        call info("[a, b] = ["//DSTR(a)//", "//DSTR(b)//"]")
        call info(":: Integração Polinomial ::")
        s = num_int(f6, a, b, n, kind="polynomial")
        call blue("I2 = ∫f(x) dx ≈ "//DSTR(s))
        call info(":: Quadratura de Gauss-Legendre ::")
        s = num_int(f6, a, b, n, kind="gauss-legendre")
        call blue("I2 = ∫f(x) dx ≈ "//DSTR(s))
        call info(":: Método de Romberg ::")
        s = num_int(f6, a, b, n, kind="romberg")
        call blue("I2 = ∫f(x) dx ≈ "//DSTR(s))

    end subroutine

    subroutine Q3
        implicit none
        integer :: n
        double precision :: a, b, r
        double precision, dimension(INT_N) :: x
        double precision, dimension(4, INT_N) :: y

        type(StringArray), dimension(:), allocatable :: legend, with

        allocate(legend(4), with(4))

        legend(1)%str = 'Polinomial'
        legend(2)%str = 'Gauss-Legendre'
        legend(3)%str = 'Romberg'
        legend(4)%str = 'Adaptativo (Gauss)'

        with(1)%str = 'linespoints'
        with(2)%str = 'linespoints'
        with(3)%str = 'linespoints'
        with(4)%str = 'lines'

        a = 0.00D0
        b = 10.0D0

        call info(ENDL//"3)"//ENDL//F7_NAME)

        call info("[a, b] = ["//DSTR(a)//", "//DSTR(b)//"]")

        INT_N = 128

        XMIN = 1.0D0
        XMAX = INT_N
        YMIN = -100.0D0
        YMAX = 300.0D0

        x = (/ (n, n=1, INT_N) /)

        call begin_plot(fname='L5-Q3')

        call subplots(2, 1)

        call info(ENDL//"m0 ~ "//F7a_NAME//ENDL)

        r = adapt_int(f7a, a, b, INT_N, tol=1.0D-8, kind="gauss-legendre")

        call info(":: Valor de referência (Integração Adaptativa) tol = 1E-8 ::")
        call blue("m0 = ∫Sσ(ω) dω ≈ "//DSTR(r))

        do n = 1, INT_N
            y(:, n) = (/ &
                num_int(f7a, a, b, n, kind="polynomial"), &
!           
                num_int(f7a, a, b, n, kind="gauss-legendre"), &
!
                num_int(f7a, a, b, n, kind="romberg"), &
!
                r &
            /)
        end do

        call info(":: Integração Polinomial ::")
        call blue("m0 = ∫Sσ(ω) dω ≈ "//DSTR(y(1, 10)))
        call info(":: Quadratura de Gauss-Legendre ::")
        call blue("m0 = ∫Sσ(ω) dω ≈ "//DSTR(y(2, 10)))
        call info(":: Método de Romberg ::")
        call blue("m0 = ∫Sσ(ω) dω ≈ "//DSTR(y(3, 10)))

        do n = 1, 4
            call subplot(1, 1, x, y(n, :), INT_N)
        end do

        call subplot_config(1, 1, title='m0 = ∫Sσ(ω) dω ≈ '//DSTR(r), xlabel='n', ylabel='m0', grid=.TRUE., &
            legend=legend, with=with, xmin=XMIN, xmax=XMAX, ymin=YMIN, ymax=YMAX)

        call info(ENDL//"m2 ~ "//F7b_NAME//ENDL)

        r = adapt_int(f7b, a, b, INT_N, tol=1.0D-8, kind="gauss-legendre")

        call info(":: Valor de referência (Integração Adaptativa) tol = 1E-8 ::")
        call blue("m2 = ∫Sσ(ω) dω ≈ "//DSTR(r))

        do n = 1, INT_N
            y(:, n) = (/ &
                num_int(f7b, a, b, n, kind="polynomial"), &
!           
                num_int(f7b, a, b, n, kind="gauss-legendre"), &
!
                num_int(f7b, a, b, n, kind="romberg"), &
!
                r &
            /)
        end do

        call info(":: Integração Polinomial ::")
        call blue("m2 = ∫Sσ(ω) dω ≈ "//DSTR(y(1, 10)))
        call info(":: Quadratura de Gauss-Legendre ::")
        call blue("m2 = ∫Sσ(ω) dω ≈ "//DSTR(y(2, 10)))
        call info(":: Método de Romberg ::")
        call blue("m2 = ∫Sσ(ω) dω ≈ "//DSTR(y(3, 10)))

        do n = 1, 4
            call subplot(2, 1, x, y(n, :), INT_N)
        end do

        call subplot_config(2, 1, title='m2 = ∫ω² Sσ(ω) dω ≈ '//DSTR(r), xlabel='n', ylabel='m2', grid=.TRUE., &
            legend=legend, with=with, xmin=XMIN, xmax=XMAX, ymin=YMIN, ymax=YMAX)

        call render_plot(clean=.TRUE.)
    end subroutine

    subroutine Q4
        implicit none
        integer :: n
        double precision :: a, b, r
        double precision, dimension(INT_N) :: x
        double precision, dimension(4, INT_N) :: y

        type(StringArray), dimension(:), allocatable :: legend, with

        allocate(legend(4), with(4))

        legend(1)%str = 'Polinomial'
        legend(2)%str = 'Gauss-Legendre'
        legend(3)%str = 'Romberg'
        legend(4)%str = 'Adaptativo (Gauss)'

        with(1)%str = 'linespoints'
        with(2)%str = 'linespoints'
        with(3)%str = 'linespoints'
        with(4)%str = 'lines'

        a = 0.00D0
        b = 10.0D0

        call info(ENDL//"4)"//ENDL//F8_NAME)

        call info("[a, b] = ["//DSTR(a)//", "//DSTR(b)//"]")

        INT_N = 128

        XMIN = 1.0D0
        XMAX = INT_N
        YMIN = -10.0D0
        YMAX = 100.0D0

        x = (/ (n, n=1, INT_N) /)

        call begin_plot(fname='L5-Q4')

        call subplots(2, 1)

        call info(ENDL//"m0 ~ "//F8a_NAME//ENDL)

        r = adapt_int(f8a, a, b, INT_N, tol=1.0D-8, kind="gauss-legendre")

        call info(":: Valor de referência (Integração Adaptativa) tol = 1E-8 ::")
        call blue("m0 = ∫Sσ(ω) dω ≈ "//DSTR(r))

        do n = 1, INT_N
            y(:, n) = (/ &
                num_int(f8a, a, b, n, kind="polynomial"), &
!           
                num_int(f8a, a, b, n, kind="gauss-legendre"), &
!
                num_int(f8a, a, b, n, kind="romberg"), &
!
                r &
            /)
        end do

        call info(":: Integração Polinomial ::")
        call blue("m0 = ∫Sσ(ω) dω ≈ "//DSTR(y(1, 10)))
        call info(":: Quadratura de Gauss-Legendre ::")
        call blue("m0 = ∫Sσ(ω) dω ≈ "//DSTR(y(2, 10)))
        call info(":: Método de Romberg ::")
        call blue("m0 = ∫Sσ(ω) dω ≈ "//DSTR(y(3, 10)))

        do n = 1, 4
            call subplot(1, 1, x, y(n, :), INT_N)
        end do

        call subplot_config(1, 1, title='m0 = ∫Sσ(ω) dω ≈ '//DSTR(r), xlabel='n', ylabel='m0', grid=.TRUE., &
            legend=legend, with=with, xmin=XMIN, xmax=XMAX, ymin=YMIN, ymax=YMAX)

        call info(ENDL//"m2 ~ "//F8b_NAME//ENDL)

        r = adapt_int(f8b, a, b, INT_N, tol=1.0D-8, kind="gauss-legendre")

        call info(":: Valor de referência (Integração Adaptativa) tol = 1E-8 ::")
        call blue("m2 = ∫Sσ(ω) dω ≈ "//DSTR(r))

        do n = 1, INT_N
            y(:, n) = (/ &
                num_int(f8b, a, b, n, kind="polynomial"), &
!           
                num_int(f8b, a, b, n, kind="gauss-legendre"), &
!
                num_int(f8b, a, b, n, kind="romberg"), &
!
                r &
            /)
        end do

        call info(":: Integração Polinomial ::")
        call blue("m2 = ∫Sσ(ω) dω ≈ "//DSTR(y(1, 10)))
        call info(":: Quadratura de Gauss-Legendre ::")
        call blue("m2 = ∫Sσ(ω) dω ≈ "//DSTR(y(2, 10)))
        call info(":: Método de Romberg ::")
        call blue("m2 = ∫Sσ(ω) dω ≈ "//DSTR(y(3, 10)))

        do n = 1, 4
            call subplot(2, 1, x, y(n, :), INT_N)
        end do

        call subplot_config(2, 1, title='m2 = ∫ω² Sσ(ω) dω ≈ '//DSTR(r), xlabel='n', ylabel='m2', grid=.TRUE., &
            legend=legend, with=with, xmin=XMIN, xmax=XMAX, ymin=YMIN, ymax=YMAX)

        call render_plot(clean=.TRUE.)
    end subroutine

    subroutine Q5
        implicit none
        integer :: n
        double precision :: a, b, s

        call info(ENDL//"5) "//F9_NAME)
        a = 0.0D0
        b = 4.0D0
        call info("[a, b] = ["//DSTR(a)//", "//DSTR(b)//"]")
        
        n = 4
        call info(":: Integração Polinomial ::")
        call blue('n = '//STR(n))
        s = num_int(f9, a, b, n, kind="polynomial")
        call blue("A = ∫f(x) dx ≈ "//DSTR(s))

        n = 2
        call info(":: Quadratura de Gauss-Legendre ::")
        call blue('n = '//STR(n))
        s = num_int(f9, a, b, n, kind="gauss-legendre")
        call blue("A = ∫f(x) dx ≈ "//DSTR(s))
    end subroutine

    subroutine Q6
        implicit none
        integer :: n
        double precision :: a, b, s

        n = 10

        call blue('n = '//STR(n))

        call info(ENDL//"6) "//F10_NAME)
        a = 0.0D0
        b = 3.0D0
        call info("[a, b] = ["//DSTR(a)//", "//DSTR(b)//"]")
        call info(":: Integração Polinomial ::")
        s = num_int(f10, a, b, n, kind="polynomial")
        call blue("A = ∫f(x) dx ≈ "//DSTR(s))
        call info(":: Quadratura de Gauss-Legendre ::")
        s = num_int(f10, a, b, n, kind="gauss-legendre")
        call blue("A = ∫f(x) dx ≈ "//DSTR(s))    
        call info(":: Método de Romberg ::")
        s = num_int(f10, a, b, n, kind="romberg")
        call blue("A = ∫f(x) dx ≈ "//DSTR(s))     
    end subroutine

    subroutine Q7
        implicit none
        integer :: n
        double precision :: a, b, r, s

        call info(ENDL//"7)")

        n = 10

        call blue('n = '//STR(n))

        call info("A1 ~ "//F11_NAME)
        a = DNINF
        b = 1.0D0
        call info("[a, b] = ["//DSTR(a)//", "//DSTR(b)//"]")
        call info(":: Quadratura de Gauss-Hermite e de Gauss-Legendre ::")
        r = num_int(f11a, a, -a, n, kind="gauss-hermite")
        s = num_int(f11b, -b, b, n, kind="gauss-legendre")
        call blue("A1 = ∫f(x) dx ≈ "//DSTR(r+s))
    
        call info("A2 ~ "//F12_NAME)
        a = DNINF
        b = DINF
        call info("[a, b] = ["//DSTR(a)//", "//DSTR(b)//"]")
        call info(":: Quadratura de Gauss-Hermite ::")
        s = num_int(f12, a, b, n, kind="gauss-hermite")
        call blue("A2 = ∫f(x) dx ≈ "//DSTR(s))
    end subroutine

    subroutine QE1
        implicit none
        double precision :: x, y, dy

        x = 3.0D0

        call info(ENDL//'1)')

        call info(FL5_QE1_NAME)
        call info(DFL5_QE1_NAME)

        call info(':: Derivada Analítica ::')
        dy = DFL5_QE1(x)
        call blue("f'("//DSTR(x)//") = "//DSTR(dy))

        call info(":: Diferenças Finitas ::"//ENDL)

        call info(':: Diferença Central (Δx = 1E-2)::')
        y = d(FL5_QE1, x, dx=1.0D-2, kind='central')
        call blue("f'("//DSTR(x)//") ≈ "//DSTR(dy))
        call blue("|δy| = "//DSTR(DABS(y - dy)))

        call info(':: Passo à frente (Δx = 1E-2)::')
        y = d(FL5_QE1, x, dx=1.0D-2, kind='forward')
        call blue("f'("//DSTR(x)//") ≈ "//DSTR(dy))
        call blue("|δy| = "//DSTR(DABS(y - dy)))

        call info(':: Passo atrás (Δx = 1E-2)::')
        y = d(FL5_QE1, x, dx=1.0D-2, kind='backward')
        call blue("f'("//DSTR(x)//") ≈ "//DSTR(dy))
        call blue("|δy| = "//DSTR(DABS(y - dy)))

        call info(ENDL//":: Extrapolação de Richard ::")

        call info(':: Diferença Central (Δx = 1E-2, p = 1)::')
        y = richard(FL5_QE1, x, dx=1.0D-2, kind='central')
        call blue("f'("//DSTR(x)//") ≈ "//DSTR(dy))
        call blue("|δy| = "//DSTR(DABS(y - dy)))

        call info(':: Diferença Central (Δx = 1E-2, p = 2)::')
        y = richard(FL5_QE1, x, dx=1.0D-2, p=2.0D0, kind='central')
        call blue("f'("//DSTR(x)//") ≈ "//DSTR(dy))
        call blue("|δy| = "//DSTR(DABS(y - dy)))

        call info(':: Passo à frente (Δx = 1E-2, p = 1)::')
        y = richard(FL5_QE1, x, dx=1.0D-2, kind='forward')
        call blue("f'("//DSTR(x)//") ≈ "//DSTR(dy))
        call blue("|δy| = "//DSTR(DABS(y - dy)))

        call info(':: Passo à frente (Δx = 1E-2, p = 2)::')
        y = richard(FL5_QE1, x, dx=1.0D-2, p=2.0D0, kind='forward')
        call blue("f'("//DSTR(x)//") ≈ "//DSTR(dy))
        call blue("|δy| = "//DSTR(DABS(y - dy)))

        call info(':: Passo atrás (Δx = 1E-2, p = 1)::')
        y = richard(FL5_QE1, x, dx=1.0D-2, kind='backward')
        call blue("f'("//DSTR(x)//") ≈ "//DSTR(dy))
        call blue("|δy| = "//DSTR(DABS(y - dy)))

        call info(':: Passo atrás (Δx = 1E-2, p = 2)::')
        y = richard(FL5_QE1, x, dx=1.0D-2, p=2.0D0, kind='backward')
        call blue("f'("//DSTR(x)//") ≈ "//DSTR(dy))
        call blue("|δy| = "//DSTR(DABS(y - dy)))

    end subroutine

    subroutine QE2
        implicit none
        double precision :: x, y, dy

        x = 2.0D0

        call info(ENDL//'2)')

        call info(FL5_QE2_NAME)
        call info(DFL5_QE2_NAME)

        call info(':: Derivada Analítica ::')
        dy = DFL5_QE2(x)
        call blue("f'("//DSTR(x)//") = "//DSTR(dy))

        call info(":: Diferenças Finitas ::"//ENDL)

        call info(':: Diferença Central (Δx = 1E-2)::')
        y = d(FL5_QE2, x, dx=1.0D-2, kind='central')
        call blue("f'("//DSTR(x)//") ≈ "//DSTR(dy))
        call blue("|δy| = "//DSTR(DABS(y - dy)))

        call info(':: Passo à frente (Δx = 1E-2)::')
        y = d(FL5_QE2, x, dx=1.0D-2, kind='forward')
        call blue("f'("//DSTR(x)//") ≈ "//DSTR(dy))
        call blue("|δy| = "//DSTR(DABS(y - dy)))

        call info(':: Passo atrás (Δx = 1E-2)::')
        y = d(FL5_QE2, x, dx=1.0D-2, kind='backward')
        call blue("f'("//DSTR(x)//") ≈ "//DSTR(dy))
        call blue("|δy| = "//DSTR(DABS(y - dy)))

        call info(ENDL//":: Extrapolação de Richard ::")

        call info(':: Diferença Central (Δx = 1E-2, p = 1)::')
        y = richard(FL5_QE2, x, dx=1.0D-2, kind='central')
        call blue("f'("//DSTR(x)//") ≈ "//DSTR(dy))
        call blue("|δy| = "//DSTR(DABS(y - dy)))

        call info(':: Diferença Central (Δx = 1E-2, p = 2)::')
        y = richard(FL5_QE2, x, dx=1.0D-2, p=2.0D0, kind='central')
        call blue("f'("//DSTR(x)//") ≈ "//DSTR(dy))
        call blue("|δy| = "//DSTR(DABS(y - dy)))

        call info(':: Passo à frente (Δx = 1E-2, p = 1)::')
        y = richard(FL5_QE2, x, dx=1.0D-2, kind='forward')
        call blue("f'("//DSTR(x)//") ≈ "//DSTR(dy))
        call blue("|δy| = "//DSTR(DABS(y - dy)))

        call info(':: Passo à frente (Δx = 1E-2, p = 2)::')
        y = richard(FL5_QE2, x, dx=1.0D-2, p=2.0D0, kind='forward')
        call blue("f'("//DSTR(x)//") ≈ "//DSTR(dy))
        call blue("|δy| = "//DSTR(DABS(y - dy)))

        call info(':: Passo atrás (Δx = 1E-2, p = 1)::')
        y = richard(FL5_QE2, x, dx=1.0D-2, kind='backward')
        call blue("f'("//DSTR(x)//") ≈ "//DSTR(dy))
        call blue("|δy| = "//DSTR(DABS(y - dy)))

        call info(':: Passo atrás (Δx = 1E-2, p = 2)::')
        y = richard(FL5_QE2, x, dx=1.0D-2, p=2.0D0, kind='backward')
        call blue("f'("//DSTR(x)//") ≈ "//DSTR(dy))
        call blue("|δy| = "//DSTR(DABS(y - dy)))
    end subroutine

    subroutine QE3
        implicit none
        double precision :: x, y, dy

        x = 6.0D0

        call info(ENDL//'3)')

        call info(FL5_QE3_NAME)
        call info(DFL5_QE3_NAME)

        call info(':: Derivada Analítica ::')
        dy = DFL5_QE3(x)
        call blue("f'("//DSTR(x)//") = "//DSTR(dy))

        call info(":: Diferenças Finitas ::"//ENDL)

        call info(':: Diferença Central (Δx = 1E-2)::')
        y = d(FL5_QE3, x, dx=1.0D-2, kind='central')
        call blue("f'("//DSTR(x)//") ≈ "//DSTR(dy))
        call blue("|δy| = "//DSTR(DABS(y - dy)))

        call info(':: Passo à frente (Δx = 1E-2)::')
        y = d(FL5_QE3, x, dx=1.0D-2, kind='forward')
        call blue("f'("//DSTR(x)//") ≈ "//DSTR(dy))
        call blue("|δy| = "//DSTR(DABS(y - dy)))

        call info(':: Passo atrás (Δx = 1E-2)::')
        y = d(FL5_QE3, x, dx=1.0D-2, kind='backward')
        call blue("f'("//DSTR(x)//") ≈ "//DSTR(dy))
        call blue("|δy| = "//DSTR(DABS(y - dy)))

        call info(ENDL//":: Extrapolação de Richard ::")

        call info(':: Diferença Central (Δx = 1E-2, p = 1)::')
        y = richard(FL5_QE3, x, dx=1.0D-2, kind='central')
        call blue("f'("//DSTR(x)//") ≈ "//DSTR(dy))
        call blue("|δy| = "//DSTR(DABS(y - dy)))

        call info(':: Diferença Central (Δx = 1E-2, p = 2)::')
        y = richard(FL5_QE3, x, dx=1.0D-2, p=2.0D0, kind='central')
        call blue("f'("//DSTR(x)//") ≈ "//DSTR(dy))
        call blue("|δy| = "//DSTR(DABS(y - dy)))

        call info(':: Passo à frente (Δx = 1E-2, p = 1)::')
        y = richard(FL5_QE3, x, dx=1.0D-2, kind='forward')
        call blue("f'("//DSTR(x)//") ≈ "//DSTR(dy))
        call blue("|δy| = "//DSTR(DABS(y - dy)))

        call info(':: Passo à frente (Δx = 1E-2, p = 2)::')
        y = richard(FL5_QE3, x, dx=1.0D-2, p=2.0D0, kind='forward')
        call blue("f'("//DSTR(x)//") ≈ "//DSTR(dy))
        call blue("|δy| = "//DSTR(DABS(y - dy)))

        call info(':: Passo atrás (Δx = 1E-2, p = 1)::')
        y = richard(FL5_QE3, x, dx=1.0D-2, kind='backward')
        call blue("f'("//DSTR(x)//") ≈ "//DSTR(dy))
        call blue("|δy| = "//DSTR(DABS(y - dy)))

        call info(':: Passo atrás (Δx = 1E-2, p = 2)::')
        y = richard(FL5_QE3, x, dx=1.0D-2, p=2.0D0, kind='backward')
        call blue("f'("//DSTR(x)//") ≈ "//DSTR(dy))
        call blue("|δy| = "//DSTR(DABS(y - dy)))
    end subroutine
end program main5