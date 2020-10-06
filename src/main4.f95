program main4
    use Func
    use Calc

    implicit none

!   Command-line Args
    integer :: argc

!   character(len=32) :: 
    integer :: i = 0

!   R -> R
    double precision :: x, y, a, b, x0
    double precision, dimension(:), allocatable :: x123

!   R^3 -> R^3
    double precision, dimension(:), allocatable :: xx(:), yy(:), xx0(:)
    double precision, dimension(:, :), allocatable :: B0(:, :)

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

200 call info(':: Zeros de funções ::')
!   Renew go to counter
    i = 0
    goto 210

210 i = i + 1
    call info('')
    go to (211, 212, 300), i

300 call info(":: Sistemas não-lineares ::")
!   Renew go to counter
    i = 0
    goto 310

310 i = i + 1
    write(*, *)
    go to (311, 400), i

400 call info(":: Integração Numérica ::")
!   Renew go to counter
    i = 0
    goto 401

401 i = i + 1
    write(*, *)
    go to (500), i

500 goto 100
!   Bissection

!   ======= Function zeros =============
211 call info("1) "//F1_NAME)
    call info(": Método da Bissecção :")
!   == Bounds definition ==
    a = -1000.0D0
    b =  1000.0D0
!   == Algorithm run ======
    x = bissection(f1, a, b)
    y = f1(x)
!   == Results ============    
    call show('x', x)
    call show('y', y)
!   =======================
    
    call info(": Método de Newton (Zero de função) :")
!   == Bounds definition ==
    x0 = DRAND(0.0D0, b)
!   == Algorithm run ======
    x = newton(f1, df1, x0)
    y = f1(x)
!   == Results ============    
    call show('x', x)
    call show('y', y)
!   =======================

    call info("Secante:")
!   == Bounds definition ==
    x0 = DRAND(0.0D0, b)
!   == Algorithm run ======
    x = secant(f1, x0)
    y = f1(x)
!   == Results ============    
    call show('x', x)
    call show('y', y)
!   =======================

    call info(": Método de Newton (Zero de função) :")
!   == Bounds definition ==
    x0 = DRAND(0.0D0, b)
!   == Algorithm run ======
    x = newton(f1, df1, x0)
    y = f1(x)
!   == Results ============    
    call show('x', x)
    call show('y', y)
!   =======================

    call info(": Método da Secante :")
!   == Bounds definition ==
    x0 = DRAND(0.0D0, b)
!   == Algorithm run ======
    x = secant(f1, x0)
    y = f1(x)
!   == Results ============    
    call show('x', x)
    call show('y', y)
!   =======================

    call info("Interpolação Inversa:")
!   == Bounds definition ==
    allocate(x123(3))
    x123 = (/-200.0D0, -100.0D0, 50.0D0/)
!   == Algorithm run ======
    x = inv_interp(f1, x123)
    y = f1(x)
!   == Results ============    
    call show('x', x)
    call show('y', y)
!   =======================
    deallocate(x123)
    goto 210

212 call info("2) "//F2_NAME)
    call info(": Método da Bissecção :")
!   == Bounds definition ==
    a = -1000.0D0
    b =  1000.0D0
!   == Algorithm run ======
    x = bissection(f2, a, b)
    y = f2(x)
!   == Results ============    
    call show('x', x)
    call show('y', y)
!   =======================
    
    call info(": Método de Newton (Zero de função) :")
!   == Bounds definition ==
    x0 = DRAND(0.0D0, b)
!   == Algorithm run ======
    x = newton(f2, df2, x0)
    y = f2(x)
!   == Results ============    
    call show('x', x)
    call show('y', y)
!   =======================

    call info("Secante:")
!   == Bounds definition ==
    x0 = DRAND(0.0D0, b)
!   == Algorithm run ======
    x = secant(f2, x0)
    y = f2(x)
!   == Results ============    
    call show('x', x)
    call show('y', y)
!   =======================

    call info(": Método de Newton (Zero de função) :")
!   == Bounds definition ==
    x0 = DRAND(0.0D0, b)
!   == Algorithm run ======
    x = newton(f2, df2, x0)
    y = f2(x)
!   == Results ============    
    call show('x', x)
    call show('y', y)
!   =======================

    call info(": Método da Secante :")
!   == Bounds definition ==
    x0 = DRAND(0.0D0, b)
!   == Algorithm run ======
    x = secant(f2, x0)
    y = f2(x)
!   == Results ============    
    call show('x', x)
    call show('y', y)
!   =======================

    call info("Interpolação Inversa:")
!   == Bounds definition ==
    allocate(x123(3))
    x123 = (/2.0D0, 5.0D0, 15.0D0/)
!   == Algorithm run ======
    x = inv_interp(f2, x123)
    y = f2(x)
!   == Results ============    
    call show('x', x)
    call show('y', y)
!   =======================
    deallocate(x123)
    goto 210

311 call info(": Método de Newton (Sistemas Não-Lineares) [Derivadas Parciais Analíticas] :")
!   == Bounds definition ==
    allocate(xx0(3))
    allocate(xx(3))
    allocate(yy(3))
    xx0 = (/DRAND(-1.0D0, 1.0D0), DRAND(-1.0D0, 1.0D0), DRAND(-1.0D0, 1.0D0)/)
!   == Algorithm run ======
    xx = sys_newton(ff1, dff1, xx0, 3)
    yy = ff1(xx, 3)
!   == Results ============    
    call info('')
    call info('x  = ')
    call print_vector(xx, 3)
    call info('')
    call info('y  = ')
    call print_vector(yy, 3)
    call info('')
!   =======================
    deallocate(xx0)
    deallocate(xx)
    deallocate(yy)
!   =======================
    
    call info(": Método de Newton (Sistemas Não-Lineares) [Derivadas Parciais Numéricas] :")
!   == Bounds definition ==
    allocate(xx0(3))
    allocate(xx(3))
    allocate(yy(3))
    xx0 = (/DRAND(-1.0D0, 1.0D0), DRAND(-1.0D0, 1.0D0), DRAND(-1.0D0, 1.0D0)/)
!   == Algorithm run ======
    xx = sys_newton_num(ff1, xx0, 3)
    yy = ff1(xx, 3)
!   == Results ============    
    call info('')
    call info('x  = ')
    call print_vector(xx, 3)
    call info('')
    call info('y  = ')
    call print_vector(yy, 3)
    call info('')
!   =======================
    deallocate(xx0)
    deallocate(xx)
    deallocate(yy)
!   ======================= 

    call info(": Método de Broyden :")
!   == Bounds definition ==
    allocate(B0(3, 3))
    allocate(xx0(3))
    allocate(xx(3))
    allocate(yy(3))
    B0 = id_matrix(3)
    xx0 = (/DRAND(-1.0D0, 1.0D0), DRAND(-1.0D0, 1.0D0), DRAND(-1.0D0, 1.0D0)/)
!   == Algorithm run ======
    xx = sys_broyden(ff1, xx0, B0, 3)
    yy = ff1(xx, 3)
!   == Results ============    
    call info('')
    call info('x  = ')
    call print_vector(xx, 3)
    call info('')
    call info('y  = ')
    call print_vector(yy, 3)
    call info('')
!   =======================
    deallocate(B0)
    deallocate(xx0)
    deallocate(xx)
    deallocate(yy)
!   =======================
    goto 310
end program main4