program main4
    use Func
    use Calc
    use Util
    implicit none

!   Command-line Args
    integer :: argc

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
11  call error('Este programa não aceita parâmetros.')
    goto 1
!   ====== Finish ====================================
1   stop
!   ==================================================

100 call Q1
    goto 200

200 call Q2
    goto 300

300 call Q3
    goto 400

400 call Q4
    goto 500

500 call Q5
    goto 10
!   ==================================================
    contains

    subroutine Q1
        implicit none
        logical :: ok
        integer :: k
        double precision :: a, b, x, y, x0
        double precision :: xx(3)
        
        call blue("1) "//F1_NAME)
        call info(": Método da Bissecção :")
    !   == Bounds definition ==
        a = -1000.0D0
        b =  1000.0D0
        call blue("[a, b] = ["//DSTR(a)//", "//DSTR(b)//"]")
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
        ok = .FALSE.
        do while(.NOT. ok)
            xx = (/ (DRAND(0.0D0, b), k=1,3)/)
    !   == Algorithm run ======
            x = inv_interp(f1, xx, ok)
            y = f1(x)
        end do
    !   == Results ============
        call blue('x1 = '//DSTR(xx(1))//'; x2 = '//DSTR(xx(2))//'; x3 = '//DSTR(xx(3))//';')
        call show('x', x)
        call show('y', y)
    !   =======================
    end subroutine

    subroutine Q2
        implicit none
        logical :: ok
        integer :: k
        double precision :: a, b, x, y, x0
        double precision :: xx(3)

        call info(ENDL//"2) "//F2_NAME)
        call info(": Método da Bissecção :")
    !   == Bounds definition ==
        a = -1000.0D0
        b =  1000.0D0
        call blue("[a, b] = ["//DSTR(a)//", "//DSTR(b)//"]")
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
            xx = (/ (DRAND(0.0D0, b), k=1,3)/)
    !   == Algorithm run ======
            x = inv_interp(f2, xx, ok)
            y = f2(x)
        end do
    !   == Results ============
        call blue('x1 = '//DSTR(xx(1))//'; x2 = '//DSTR(xx(2))//'; x3 = '//DSTR(xx(3))//';')
        call show('x', x)
        call show('y', y)
    !   =======================
    end subroutine

    subroutine Q3
        implicit none
        logical :: ok
        integer :: k
        double precision, dimension(F3_N) :: x, y, x0
        double precision, dimension(F3_N, F3_N) :: J0
        
        call blue(ENDL//'3) '//F3_NAME)
        x0 = rand_vector(F3_N, 0.0D0, 1.0D0)
        call info(": Método de Newton (Sistemas Não-Lineares) [Derivadas Parciais Analíticas] :")
    !   == Bounds definition ==
        do k=1, D_MAX_ITER
    !   == Algorithm run ======
            x = sys_newton(f3, df3, x0, F3_N, ok)
            if (.NOT. ok) then
                x0 = rand_vector(F3_N, 0.0D0, 1.0D0)
            else
                exit
            end if
        end do
        if (.NOT. ok) then
            call error('Este método não convergiu.')
        else
            y = f3(x, F3_N)
        !   == Results ============    
            call show_vector('x0', x0, F3_N)
            call show_vector('x', x, F3_N)
            call show_vector('y', y, F3_N)
        !   =======================
        end if
        
        call info(": Método de Newton (Sistemas Não-Lineares) [Derivadas Parciais Numéricas] :")
    !   == Bounds definition ==
        do k=1, D_MAX_ITER
    !   == Algorithm run ======
            x = sys_newton_num(f3, x0, F3_N, ok)
            if (.NOT. ok) then
                x0 = rand_vector(F3_N, 0.0D0, 1.0D0)
            else
                exit
            end if
        end do
        if (.NOT. ok) then
            call error('Este método não convergiu.')
        else
            y = f3(x, F3_N)
        !   == Results ============    
            call show_vector('x0', x0, F3_N)
            call show_vector('x', x, F3_N)
            call show_vector('y', y, F3_N)
        !   =======================
        end if

        call info(": Método de Broyden :")
    !   == Bounds definition ==
        J0 = id_matrix(F3_N)
        do k=1, D_MAX_ITER
    !   == Algorithm run ======
            x = sys_broyden(f3, x0, J0, F3_N, ok)
            if (.NOT. ok) then
                x0 = rand_vector(F3_N, 0.0D0, 1.0D0)
            else
                exit
            end if
        end do
        if (.NOT. ok) then
            call error('Este método não convergiu.')
        else
            y = f3(x, F3_N)
        !   == Results ============    
            call show_vector('x0', x0, F3_N)
            call show_vector('x', x, F3_N)
            call show_vector('y', y, F3_N)
        !   =======================
        end if
    end subroutine

    subroutine Q4
        implicit none
        logical :: ok
        integer :: i, k
        double precision, dimension(3) :: x, y, x0
        double precision, dimension(3, 3) :: J0

        call info('4) '//F4_NAME)

        do i=1, F4_N
            F4_T1 = F4_TT1(i)
            F4_T2 = F4_TT2(i)

            call show('θ1', F4_T1)
            call show('θ2', F4_T2)

            call info(": Método de Newton (Sistemas Não-Lineares) [Derivadas Parciais Analíticas] :")
            x0 = rand_vector(F4_N, 0.0D0, 1.0D0)
            do k=1, 5 * D_MAX_ITER
        !   == Algorithm run ======
                x = sys_newton(f4, df4, x0, F4_N, ok)
                if (.NOT. ok) then
                    x0 = rand_vector(F4_N, 0.0D0, 1.0D0)
                else
                    exit
                end if
            end do
            if (.NOT. ok) then
                call error('Este método não convergiu.')
            else
                y = f4(x, F4_N)
            !   == Results ============    
                call show_vector('x0', x0, F4_N)
                call show_vector('x', x, F4_N)
                call show_vector('y', y, F4_N)
            !   =======================
            end if
            
            call info(": Método de Newton (Sistemas Não-Lineares) [Derivadas Parciais Numéricas] :")
            do k=1, 5 * D_MAX_ITER
        !   == Algorithm run ======
                x = sys_newton_num(f4, x0, F4_N, ok)
                if (.NOT. ok) then
                    x0 = rand_vector(F4_N, 0.0D0, 1.0D0)
                else
                    exit
                end if
            end do
            if (.NOT. ok) then
                call error('Este método não convergiu.')
            else
                y = f4(x, F4_N)
            !   == Results ============    
                call show_vector('x0', x0, F4_N)
                call show_vector('x', x, F4_N)
                call show_vector('y', y, F4_N)
            !   =======================
            end if

            call info(": Método de Broyden :")
            J0 = id_matrix(F4_N)
            do k=1, 5 * D_MAX_ITER
        !   == Algorithm run ======
                x = sys_broyden(f4, x0, J0, F4_N, ok)
                if (.NOT. ok) then
                    x0 = rand_vector(F4_N, -1.0D0, 1.0D0)
                else
                    exit
                end if
            end do
            if (.NOT. ok) then
                call error('Este método não convergiu.')
            else
                y = f4(x, F4_N)
            !   == Results ============    
                call show_vector('x0', x0, F4_N)
                call show_vector('x', x, F4_N)
                call show_vector('y', y, F4_N)
            !   =======================
            end if
        end do
    end subroutine

    subroutine Q5
        implicit none
        logical :: ok
        integer :: k
        double precision, dimension(3) :: x, y, b, b0
        
        call info("5) "//F5_NAME)

    !   ===============================
        x = (/1.0D0, 2.0D0, 3.0D0/)
        y = (/1.0D0, 2.0D0, 9.0D0/)

        b0 = rand_vector(F5_N, -1.0D0, 1.0D0)
    
        call info(": Método não-linear de Mínimos Quadrados :")
        do k=1, D_MAX_ITER
    !   == Algorithm run ======
            b = sys_least_squares(f5, df5, x, y, b0, F5_N, F5_N, ok)
            if (.NOT. ok) then
                b0 = rand_vector(F5_N, 0.8D0, 1.0D0)
            else
                exit
            end if
        end do
        if (.NOT. ok) then
            call error('Este método não convergiu.')
        else
            y = f5(x, b, F5_N, F5_N)
        !   == Results ============    
            call show_vector('b0', b0, F5_N)
            call show_vector('b', b, F5_N)
            call show_vector('x', x, F5_N)
            call show_vector('y', y, F5_N)
        !   =======================
        end if

        call info(": Método não-linear de Mínimos Quadrados [Derivadas Numéricas]:")
        do k=1, D_MAX_ITER
    !   == Algorithm run ======
            b = sys_least_squares(f5, df5, x, y, b0, F5_N, F5_N, ok)
            if (.NOT. ok) then
                b0 = rand_vector(F5_N, 0.8D0, 1.0D0)
            else
                exit
            end if
        end do
        if (.NOT. ok) then
            call error('Este método não convergiu.')
        else
            y = f5(x, b, F5_N, F5_N)
        !   == Results ============    
            call show_vector('b0', b0, F5_N)
            call show_vector('b', b, F5_N)
            call show_vector('x', x, F5_N)
            call show_vector('y', y, F5_N)
        !   =======================
        end if
    end subroutine
end program main4