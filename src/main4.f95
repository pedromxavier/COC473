program main4
    use Func
    use Calc
    use Util
    implicit none

!   Command-line Args
    integer :: argc

!   character(len=32) :: 
    integer :: i = 0, j = 0, k

    logical :: ok

!   R -> R
    double precision :: x, y, a, b, x0
    double precision, dimension(:), allocatable :: x123

!   R^3 -> R^3
    double precision, dimension(:), allocatable :: xx(:), yy(:), bb(:), xx0(:), bb0(:)
    double precision, dimension(:, :), allocatable :: J0(:, :)

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
    goto 201

201 i = i + 1
    call info('')
    go to (210, 220, 300), i

300 call info(":: Sistemas não-lineares ::")
!   Renew go to counter
    i = 0
    goto 301

301 i = i + 1
    write(*, *)
    go to (310, 320, 330, 400), i

400 goto 100
!   Bissection

!   ======= Function zeros =============
210 call info("1) "//F1_NAME)
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
    ok = .FALSE.
    do while(.NOT. ok)
        x0 = DRAND(0.0D0, b)
!   == Algorithm run ======
        x = newton(f1, df1, x0, ok)
        y = f1(x)
    end do
!   == Results ============
    call show('x0', x0)    
    call show('x', x)
    call show('y', y)
!   =======================

    call info(": Método da Secante :")
!   == Bounds definition ==
    ok = .FALSE.
    do while(.NOT. ok)
        x0 = DRAND(0.0D0, b)
!   == Algorithm run ======
        x = secant(f1, x0, ok)
        y = f1(x)
    end do
!   == Results ============
    call show('x0', x0)    
    call show('x', x)
    call show('y', y)
!   =======================

    call info("Interpolação Inversa:")
!   == Bounds definition ==
    allocate(x123(3))
    ok = .FALSE.
    do while(.NOT. ok)
        x123 = (/ (DRAND(0.0D0, b), k=1,3)/)
!   == Algorithm run ======
        x = inv_interp(f1, x123, ok)
        y = f1(x)
    end do
!   == Results ============
    call show('x1', x123(1))
    call show('x2', x123(2))
    call show('x3', x123(3))
    call show('x', x)
    call show('y', y)
!   =======================
!   deallocate(x123)
    goto 201

220 call info("2) "//F2_NAME)
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
    ok = .FALSE.
    do while(.NOT. ok)
        x0 = DRAND(0.0D0, b)
!   == Algorithm run ======
        x = newton(f2, df2, x0, ok)
        y = f2(x)
    end do
!   == Results ============
    call show('x0', x0)
    call show('x', x)
    call show('y', y)
!   =======================

    call info(": Método da Secante :")
!   == Bounds definition ==
    ok = .FALSE.
    do while(.NOT. ok)
        x0 = DRAND(0.0D0, b)
!   == Algorithm run ======
        x = secant(f2, x0, ok)
        y = f2(x)
    end do
!   == Results ============
    call show('x0', x0)
    call show('x', x)
    call show('y', y)
!   =======================

    call info("Interpolação Inversa:")
!   == Bounds definition ==
!   allocate(x123(3))
    ok = .FALSE.
    do while(.NOT. ok)
        x123 = (/ (DRAND(0.0D0, b), k=1,3)/)
!   == Algorithm run ======
        x = inv_interp(f2, x123, ok)
        y = f2(x)
    end do
!   == Results ============
    call show('x1', x123(1))
    call show('x2', x123(2))
    call show('x3', x123(3))   
    call show('x', x)
    call show('y', y)
!   =======================
    deallocate(x123)
    goto 201

310 call info('3) '//F3_NAME)
    call info(": Método de Newton (Sistemas Não-Lineares) [Derivadas Parciais Analíticas] :")
!   == Bounds definition ==
    allocate(xx0(3))
    allocate(xx(3))
    allocate(yy(3))
    ok = .FALSE.
    do while (.NOT. ok)
        xx0 = (/ (DRAND(0.0D0, 1.0D0), k=1,3) /)
!   == Algorithm run ======
        xx = sys_newton(ff1, dff1, xx0, 3, ok)
        yy = ff1(xx, 3)
    end do
!   == Results ============
    call show_vector('x0', xx0, 3)
    call show_vector('x', xx, 3)
    call show_vector('y', yy, 3)
!   =======================
!   deallocate(xx0)
!   deallocate(xx)
!   deallocate(yy)
!   =======================
    
    call info(": Método de Newton (Sistemas Não-Lineares) [Derivadas Parciais Numéricas] :")
!   == Bounds definition ==
!   allocate(xx0(3))
!   allocate(xx(3))
!   allocate(yy(3))
    ok = .FALSE.
    do while (.NOT. ok)
        xx0 = (/ (DRAND(0.0D0, 1.0D0), k=1,3) /)
!   == Algorithm run ======
        xx = sys_newton_num(ff2, xx0, 3, ok)
    end do
    yy = ff2(xx, 3)
!   == Results ============    
    call show_vector('x0', xx0, 3)
    call show_vector('x', xx, 3)
    call show_vector('y', yy, 3)
!   =======================
!   deallocate(xx0)
    deallocate(xx)
    deallocate(yy)
!   ======================= 

    call info(": Método de Broyden :")
!   == Bounds definition ==
    allocate(J0(3, 3))
!   allocate(xx0(3))
!   allocate(xx(3))
!   allocate(yy(3))
    J0 = id_matrix(3)
    ok = .FALSE.
    do while (.NOT. ok)
        xx0 = (/ (DRAND(0.0D0, 1.0D0), k=1,3) /)
!   == Algorithm run ======
        xx = sys_broyden(ff2, xx0, J0, 3, ok)
    end do
    yy = ff2(xx, 3)
!   == Results ============    
    call show_vector('x0', xx0, 3)
    call show_vector('x', xx, 3)
    call show_vector('y', yy, 3)
!   =======================
!   deallocate(J0)
!   deallocate(xx0)
!   deallocate(xx)
!   deallocate(yy)
!   =======================
    goto 301

320 call info('4) '//F4_NAME)
    j = 0
    goto 321

321 j = j + 1
    goto (322, 323, 324, 301), j

322 t1 = 0.00D0
    t2 = 3.00D0
    call show('t1', t1)
    call show('t2', t2)
    goto 325

323 t1 = 0.75D0
    t2 = 6.50D0
    call show('t1', t1)
    call show('t2', t2)
    goto 325

324 t1 =  0.000D0
    t2 = 11.667D0
    call show('t1', t1)
    call show('t2', t2)
    goto 325

325 call info(": Método de Newton (Sistemas Não-Lineares) [Derivadas Parciais Analíticas] :")
!   == Bounds definition ==
!   allocate(xx0(3))
!   allocate(xx(3))
!   allocate(yy(3))
    ok = .FALSE.
    do while(.NOT. ok)
        xx0 = (/ (DRAND(0.0D0, 1.0D0), k=1,3) /)
!   == Algorithm run ======
        xx = sys_newton(ff2, dff2, xx0, 3, ok)
        yy = ff2(xx, 3)
    end do
!   == Results ============    
    call show_vector('x0', xx0, 3)
    call show_vector('x', xx, 3)
    call show_vector('y', yy, 3)
!   =======================
!   deallocate(xx0)
!   deallocate(xx)
!   deallocate(yy)
!   =======================
    
    call info(": Método de Newton (Sistemas Não-Lineares) [Derivadas Parciais Numéricas] :")
!   == Bounds definition ==
!   allocate(xx0(3))
!   allocate(xx(3))
!   allocate(yy(3))
    ok = .FALSE.
    do while (.NOT. ok)
        xx0 = (/ (DRAND(0.0D0, 1.0D0), k=1,3) /)
!   == Algorithm run ======
        xx = sys_newton_num(ff2, xx0, 3, ok)
    end do
    yy = ff2(xx, 3)
!   == Results ============    
    call show_vector('x0', xx0, 3)
    call show_vector('x', xx, 3)
    call show_vector('y', yy, 3)
!   =======================
!   deallocate(xx0)
!   deallocate(xx)
!   deallocate(yy)
!   ======================= 

    call info(": Método de Broyden :")
!   == Bounds definition ==
!   allocate(J0(3, 3))
!   allocate(xx0(3))
!   allocate(xx(3))
!   allocate(yy(3))
    J0 = id_matrix(3)
    ok = .FALSE.
    do while (.NOT. ok)
        xx0 =  (/ (DRAND(0.0D0, 1.0D0), k=1,3) /)
!   == Algorithm run ======
        xx = sys_broyden(ff2, xx0, J0, 3, ok)
    end do
    yy = ff2(xx, 3)
!   == Results ============    
    call show_vector('x0', xx0, 3)
    call show_vector('x', xx, 3)
    call show_vector('y', yy, 3)
!   =======================
    deallocate(J0)
    deallocate(xx0)
    deallocate(xx)
    deallocate(yy)
!   =======================
    goto 321

330 call info("5) "//F5_NAME)
!   ===============
    allocate(xx(3))
    allocate(yy(3))
    allocate(bb(3))
    allocate(bb0(3))
!   ===============================
    xx(:) = (/1.0D0, 2.0D0, 3.0D0/)
    yy(:) = (/1.0D0, 2.0D0, 9.0D0/)
 
    call info(": Método não-linear de Mínimos Quadrados :")

    ok = .FALSE.

    bb0(:) = (/ 0.293477D0, &
                0.978723D0, &
                0.842565D0 /)
    bb(:) = sys_least_squares(ff5, dff5, xx, yy, bb0, 3, 3, ok)
    
    call show_vector('x', xx, 3)
    call show_vector('y', yy, 3)
    call show_vector('b0', bb0, 3)
    call show_vector('b', bb, 3)

    call info(": Método não-linear de Mínimos Quadrados [Derivadas Numéricas]:")

    ok = .FALSE.

    bb0(:) = (/ 0.293477D0, &
                0.978723D0, &
                0.842565D0 /)
    bb(:) = sys_least_squares_num(ff5, xx, yy, bb0, 3, 3, ok)

    call show_vector('x', xx, 3)
    call show_vector('y', yy, 3)
    call show_vector('b0', bb0, 3)
    call show_vector('b', bb, 3)
!   =========================
    deallocate(xx)
    deallocate(yy)
    deallocate(bb)
    deallocate(bb0)
!   =========================

    goto 301
end program main4