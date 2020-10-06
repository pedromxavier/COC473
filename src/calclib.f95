!   Calc Module

    module Calc
        implicit none
        integer :: MAX_ITER = 100000
        double precision :: h = 1.0D-5
        double precision :: TOL = 1.0D-7
    contains
        subroutine init_random_seed()
            integer :: i, n, clock
            integer, allocatable :: seed(:)
        
            call RANDOM_SEED(SIZE=n)
            allocate(seed(n))
            call SYSTEM_CLOCK(COUNT=clock)
            seed = clock + 37 * (/ (i - 1, i = 1, n) /)
            call RANDOM_SEED(PUT=seed)
            deallocate(seed)
        end subroutine

!       ===== I/O Metods =====
        subroutine error(text)
!           Red Text
            implicit none
            character(len=*) :: text
            write (*, *) ''//achar(27)//'[31m'//text//''//achar(27)//'[0m'
        end subroutine

        subroutine warn(text)
!           Yellow Text
            implicit none
            character(len=*) :: text
            write (*, *) ''//achar(27)//'[93m'//text//''//achar(27)//'[0m'
        end subroutine

        subroutine info(text)
!           Green Text
            implicit none
            character(len=*) :: text
            write (*, *) ''//achar(27)//'[32m'//text//''//achar(27)//'[0m'
        end subroutine

        subroutine show(var, value)
!           Violet Text
            implicit none
            character(len=*) :: var
            character(len=24) :: val
            double precision :: value

10          format(F24.16, '')

            write (val, 10) value
            write (*, *) ''//achar(27)//'[36m'//var//' = '//val//''//achar(27)//'[0m'
        end subroutine

        function id_matrix(n) result (A)
            implicit none

            integer :: n
            double precision :: A(n, n)
            integer :: j

            A(:, :) = 0.0D0

            do j = 1, n
                A(j, j) = 1.0D0
            end do
            return
        end function

        subroutine print_matrix(A, m, n)
            implicit none

            integer :: m, n
            double precision :: A(m, n)

            integer :: i, j

20          format(' |', F32.12, ' ')
21          format(F30.12, '|')
22          format(F30.12, ' ')

            do i = 1, m
                do j = 1, n
                    if (j == 1) then
                        write(*, 20, advance='no') A(i, j)
                    elseif (j == n) then
                        write(*, 21, advance='yes') A(i, j)
                    else
                        write(*, 22, advance='no') A(i, j)
                    end if
                end do
            end do
        end subroutine
        
        subroutine read_matrix(fname, A, m, n)
            implicit none
            character(len=*) :: fname
            integer :: m, n
            double precision, allocatable :: A(:, :)

            integer :: i

            open(unit=33, file=fname, status='old', action='read')
            read(33, *) m
            read(33, *) n
            allocate(A(m, n))

            do i = 1, m
                read(33,*) A(i,:)
            end do

            close(33)
        end subroutine

        subroutine print_vector(x, n)
            implicit none

            integer :: n
            double precision :: x(n)
            
            integer :: i

30          format(' |', F30.12, '|')

            do i = 1, n
                write(*, 30) x(i)
            end do
        end subroutine

        subroutine read_vector(fname, b, t)
            implicit none
            character(len=*) :: fname
            integer :: t
            double precision, allocatable :: b(:)

            open(unit=33, file=fname, status='old', action='read')
            read(33, *) t
            allocate(b(t))

            read(33,*) b(:)

            close(33)
        end subroutine

        function NORM(x, n) result (s)
            implicit none

            integer :: n
            double precision :: x(n)
            double precision :: s

            s = SQRT(DOT_PRODUCT(x, x))
            return
        end function

!       ================= Matrix Method ====================
        function inv(A, n) result (Ainv)
            integer :: n
            double precision :: A(n, n), Ainv(n, n)
            double precision :: work(n)
            integer :: ipiv(n)   ! pivot indices
            integer :: info
          
            ! External procedures defined in LAPACK
            external DGETRF
            external DGETRI
          
            ! Store A in Ainv to prevent it from being overwritten by LAPACK
            Ainv(:, :) = A(:, :)
          
            ! DGETRF computes an LU factorization of a general M-by-N matrix A
            ! using partial pivoting with row interchanges.
            call DGETRF(n, n, Ainv, n, ipiv, info)
          
            if (info /= 0) then
                call print_matrix(A, n, n)
                stop 'Matrix is numerically singular!'
            end if
          
            ! DGETRI computes the inverse of a matrix using the LU factorization
            ! computed by DGETRF.
            call DGETRI(n, Ainv, n, ipiv, work, n, info)
          
            if (info /= 0) then
                call print_matrix(A, n, n)
                stop 'Matrix inversion failed!'
            end if

            return
        end function

        recursive subroutine cross_quick_sort(x, y, u, v, n)
            integer :: n, i, j, u, v
            double precision :: p, aux, auy
            double precision :: x(n), y(n)

            i = u
            j = v

            p = x((u + v) / 2)

            do while (i <= j)
                do while (x(i) < p)
                    i = i + 1
                end do
                do while(x(j) > p)
                    j = j - 1
                end do
                if (i <= j) then
                    aux = x(i)
                    auy = y(i)
                    x(i) = x(j)
                    y(i) = y(j)
                    x(j) = aux
                    y(j) = auy
                    i = i + 1
                    j = j - 1
                end if
            end do

            if (u < j) then
                call cross_quick_sort(x, y, u, j, n)
            end if
            if (i < v) then
                call cross_quick_sort(x, y, i, v, n)
            end if

            return
        end subroutine

        subroutine cross_sort(x, y, n)
            implicit none
            integer :: n
            double precision :: x(n), y(n)
            
            call cross_quick_sort(x, y, 1, n, n)
        end subroutine

        function DRAND(a, b) result (y)
            implicit none
            double precision :: a, b, x, y
            ! x in [0, 1)            
            call RANDOM_NUMBER(x)
            y = (x * (b - a)) + a
            return
        end function

        function outer_product(x, y, n) result(A)
            implicit none
            integer :: n
            double precision, dimension(n), intent(in) :: x, y
            double precision, dimension(n, n) :: A
            integer :: i, j

            do i=1,n
                do j=1,n
                    A(i, j) = x(i) * y(j)
                end do
            end do
            return
        end function

!       ================= Numerical Mathods =================

        function d(f, x) result (y)
            implicit none
            double precision :: f
            double precision :: x, y
        
            y = (f(x + h) - f(x - h)) / (2 * h)
            return
        end function

        function dp(f, x, i, n) result (y)
            implicit none
            integer :: i, n
            double precision :: f
            double precision :: x(n), xh(n)
            double precision :: y

            xh(:) = 0.0D0
            xh(i) = h

            y = (f(x + xh) - f(x - xh)) / (2 * h)
            return
        end function

        function grad(f, x, n) result (y)
            implicit none
            integer :: i, n
            double precision :: f
            double precision :: xh(n), x(n), y(n)

            xh(:) = 0.0D0
            do i=1, n
!               Compute partial derivative with respect to x_i                
                xh(i) = h
                y(i) = (f(x + xh) - f(x - xh)) / (2 * h)
                xh(i) = 0.0D0
            end do
            return
        end function

!       =====================================================

        function lagrange(x0, y0, n, x) result (y)
            implicit none
            integer :: n
            double precision :: x0(n), y0(n)
            double precision :: x, y, yi
            integer :: i, j

            y = 0.0D0
            do i = 1, n
                yi = y0(i)
                do j = 1, n
                    if (i /= j) then
                        yi = yi * (x - x0(j)) / (x0(i) - x0(j))
                    end if
                end do
                y = y + yi
            end do

            return
        end function

        function bissection(f, aa, bb) result (x)
            implicit none
            double precision :: f
            double precision, intent(in) :: aa, bb
            double precision :: a, b, x

            if (bb < aa) then
                a = bb
                b = aa
            else
                a = aa
                b = bb
            end if

            do while (DABS(a - b) > TOL)
                x = (a + b) / 2
                if (f(a) > f(b)) then
                    if (f(x) > 0) then
                        a = x
                    else
                        b = x
                    end if
                else 
                    if (f(x) < 0) then
                        a = x
                    else
                        b = x
                    end if
                end if
            end do
            x = (a + b) / 2
            return
        end function

        function newton(f, df, x0) result (x)
            implicit none
            integer :: i

            double precision :: f, df
            double precision, intent(in) :: x0
            double precision :: x, xk

            xk = x0

            do i = 1, MAX_ITER
                x = xk - f(xk) / df(xk)
                if (DABS(x - xk) > TOL) then
                    xk = x
                else
                    return
                end if 
            end do
            call error("O método de Newton não convergiu.")
            return
        end function

        function secant(f, x0) result (x)
            implicit none
            integer :: i
            double precision :: xk(3), yk(2)
            double precision, intent(in) :: x0
            double precision :: x
            interface
                function f(x) result (y)
                    implicit none
                    double precision :: x, y
                end function
            end interface

            xk(1) = x0
            xk(2) = x0 + h
            yk(1) = f(xk(1))
            do i = 1, MAX_ITER
                yk(2) = f(xk(2))
                xk(3) = xk(2) - (yk(2) * (xk(2) - xk(1))) / (yk(2) - yk(1))             
                if (DABS(xk(3) - xk(2)) > TOL) then
                    xk(1:2) = xk(2:3)
                    yk(1) = yk(2)
                else
                    x = xk(3)
                    return
                end if 
            end do

            call error("O método da Secante não convergiu.")
            return
        end function

        function inv_interp(f, x00) result (x)
            implicit none
            double precision :: f
            double precision :: x, xk
            double precision, intent(in) :: x00(3)
            double precision :: x0(3), y0(3)
            integer :: i, k
            integer :: j(1)

            x0(:) = x00(:)
            xk = 1.0D+32

            do k = 1, MAX_ITER
                call cross_sort(x0, y0, 3)

!               Cálculo de y                
                do i = 1, 3
                    y0(i) = f(x0(i))
                end do

                x = lagrange(y0, x0, 3, 0.0D0)

                if (DABS(x - xk) > TOL) then
                    j(:) = MAXLOC(DABS(y0))
                    i = j(1)
                    x0(i) = x
                    y0(i) = f(x)
                    xk = x
                else
                    return
                end if
            end do

            call error("O Método de Interpolação Inversa não convergiu.")

            return
        end function

        function sys_newton(ff, dff, xx0, n) result (x)
            implicit none
            integer :: n
            double precision, dimension(n), intent(in) :: xx0
            double precision, dimension(n) :: x, xk, dx
            double precision :: J(n, n) 
            integer :: k
            interface
                function ff(x, n) result (y)
                    implicit none
                    integer :: n
                    double precision :: x(n), y(n)
                end function
            end interface

            interface
                function dff(x, n) result (J)
                    implicit none
                    integer :: n
                    double precision :: x(n), J(n, n)
                end function
            end interface

            xk(:) = xx0(:)

            do k=1, MAX_ITER
                J(:, :) = dff(xk, n)
                dx(:) = -MATMUL(inv(J, n), ff(xk, n))
                x(:) = xk(:) + dx(:)

                if ((NORM(dx, n) / NORM(x, n)) > TOL) then
                    xk(:) = x(:)
                else
                    return
                end if
            end do
            call error("O Método de Newton não convergiu.")
            return
        end function

        function sys_newton_num(ff, xx0, n) result (x)
!           Same as previous function, with numerical partial derivatives            
            implicit none
            integer :: n
            double precision, dimension(n), intent(in) :: xx0
            double precision, dimension(n):: x, xk, xh, dx, y
            double precision, dimension(n, n) :: J
            integer :: i, k
            interface
                function ff(x, n) result (y)
                    implicit none
                    integer :: n
                    double precision :: x(n), y(n)
                end function
            end interface

            xh(:) = 0.0D0
            xk(:) = xx0(:)

            do k=1, MAX_ITER
!               Compute Jacobian Matrix
                J(:, :) = 0.0D0
                do i=1, n
!                   Partial derivative with respect do the i-th coordinates                    
                    xh(i) = h
                    J(:, i) = (ff(x(:) + xh(:), n) - ff(x(:) - xh(:), n)) / (2 * h)
                    xh(i) = 0.0D0
                end do

                y(:) = ff(xk, n)
                dx(:) = -matmul(inv(J, n), y)
                x(:) = xk(:) + dx(:)

                if ((NORM(dx, n) / NORM(x, n)) > TOL) then
                    xk(:) = x(:)
                else
                    return
                end if
            end do
            call error("O Método de Newton não convergiu. (Derivada Numérica)")
            return
        end function

        function sys_broyden(ff, xx0, B0, n) result (x)
            implicit none
            integer :: n
            double precision, dimension(n), intent(in) :: xx0
            double precision, dimension(n, n), intent(in) :: B0
            double precision, dimension(n) :: x, xk, yk, dx
            double precision, dimension(n, n) :: J, Bk
            integer :: k
            interface
                function ff(x, n) result (y)
                    implicit none
                    integer :: n
                    double precision :: x(n), y(n)
                end function
            end interface

            xk(:) = xx0(:)
            Bk(:, :) = B0(:, :)

            do k=1, MAX_ITER
                J(:, :) = Bk(:, :)
                dx(:) = -MATMUL(inv(J, n), ff(xk, n))
                x(:) = xk(:) + dx(:)
                yk(:) = ff(x, n) - ff(xk, n)
                if ((norm(dx, n) / norm(x, n)) > TOL) then
                    Bk(:, :) = Bk(:, :) + (outer_product(yk(:) - MATMUL(Bk, dx), dx, n)  /  DOT_PRODUCT(dx, dx))
                    xk(:) = x(:)
                else
                    return
                end if
            end do
            call error("O Método de Broyden não convergiu")
            return
        end function

        function sys_least_squares(ff, dff, xx0, n) result (x)
            implicit none
            integer :: n
            double precision, dimension(n), intent(in) :: xx0
            double precision, dimension(n) :: x, xk, dx
            double precision :: J(n, n) 
            integer :: k
            interface
                function ff(x, n) result (y)
                    implicit none
                    integer :: n
                    double precision :: x(n), y(n)
                end function
            end interface

            interface
                function dff(x, n) result (J)
                    implicit none
                    integer :: n
                    double precision :: x(n), J(n, n)
                end function
            end interface

            xk(:) = xx0(:)

            do k=1, MAX_ITER
                J(:, :) = dff(xk, n)
                dx(:) = -MATMUL(inv(MATMUL(TRANSPOSE(J), J), n), MATMUL(TRANSPOSE(J), ff(xk, n)))
                x(:) = xk(:) + dx(:)

                if ((NORM(dx, n) / NORM(x, n)) > TOL) then
                    xk(:) = x(:)
                else
                    return
                end if
            end do
            call error("O Método (Não-linear) de Mínimos Quadrados não convergiu.")
            return
        end function

        function sys_least_squares_num(ff, xx0, n) result (x)
!           Same as previous function, with numerical partial derivatives            
            implicit none
            integer :: n
            double precision, dimension(n), intent(in) :: xx0
            double precision, dimension(n):: x, xk, xh, dx
            double precision, dimension(n, n) :: J
            integer :: i, k
            interface
                function ff(x, n) result (y)
                    implicit none
                    integer :: n
                    double precision :: x(n), y(n)
                end function
            end interface

            xh(:) = 0.0D0
            xk(:) = xx0(:)

            do k=1, MAX_ITER
!               Compute Jacobian Matrix
                J(:, :) = 0.0D0
                do i=1, n
!                   Partial derivative with respect do the i-th coordinates                    
                    xh(i) = h
                    J(:, i) = (ff(x(:) + xh(:), n) - ff(x(:) - xh(:), n)) / (2 * h)
                    xh(i) = 0.0D0
                end do
                dx(:) = -MATMUL(inv(MATMUL(TRANSPOSE(J), J), n), MATMUL(TRANSPOSE(J), ff(xk, n)))
                x(:) = xk(:) + dx(:)

                if ((NORM(dx, n) / NORM(x, n)) > TOL) then
                    xk(:) = x(:)
                else
                    return
                end if
            end do
            call error("O Método (Não-Linear) de Mínimos Quadrados não convergiu. (Derivada Numérica)")
            return
        end function
    end module Calc
