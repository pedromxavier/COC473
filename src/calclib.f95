!   Calc Module

    module Calc
        use Util
        use Matrix
        implicit none
        double precision :: h = 1.0D-5
    contains
!       ================= Numerical Mathods =================
        function d(f, x, kind) result (y)
            implicit none
            character (len=*), optional :: kind
            double precision, intent(in) :: x
            double precision :: y
            interface
                function f(x) result (y)
                    implicit none
                    double precision, intent(in) :: x
                    double precision :: y
                end function
            end interface

            if (.NOT. PRESENT(kind)) then
                kind = "central"
            end if

            if (kind == "central") then
                y = (f(x + h) - f(x - h)) / (2 * h)
            else if (kind == "forward") then
                y = (f(x + h) - f(x)) / h
            else if (kind == "backward") then
                y = (f(x) - f(x - h)) / h
            else 
                call error("Unexpected value `"//kind//" for derivative kind."// &
                "Options are: `central`, `forward` and `backward`.")
            end if
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
            logical :: ok = .TRUE.
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
                dx(:) = -MATMUL(inv(J, n, ok), ff(xk, n))
                x(:) = xk(:) + dx(:)
                
                if (.NOT. ok) then
                    call error("Tentativa de inversão de matriz singular.")
                    exit
                else if ((NORM(dx, n) / NORM(x, n)) > TOL) then
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
            logical :: ok = .TRUE.
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
                dx(:) = -matmul(inv(J, n, ok), y)
                x(:) = xk(:) + dx(:)

                if (.NOT. ok) then
                    call error("Tentativa de inversão de matriz singular.")
                    exit
                else if ((NORM(dx, n) / NORM(x, n)) > TOL) then
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
            logical :: ok = .TRUE.
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
                dx(:) = -MATMUL(inv(J, n, ok), ff(xk, n))
                x(:) = xk(:) + dx(:)
                yk(:) = ff(x, n) - ff(xk, n)
                if (.NOT. ok) then
                    call error("Tentativa de inversão de matriz singular.")
                    exit
                else if ((norm(dx, n) / norm(x, n)) > TOL) then
                    Bk(:, :) = Bk(:, :) + (outer_product(yk(:) - MATMUL(Bk, dx), dx, n)  /  DOT_PRODUCT(dx, dx))
                    xk(:) = x(:)
                else
                    return
                end if
            end do
            call error("O Método de Broyden não convergiu")
            return
        end function

        function sys_least_squares(ff, dff, x, y, b0, m, n) result (b)
            implicit none
            integer :: m, n
            double precision, dimension(n), intent(in) :: x, y, b0
            double precision, dimension(n) :: b, bk, db
            double precision :: J(n, n) 
            integer :: k
            logical :: ok = .TRUE.
            interface
                function ff(x, b, m, n) result (z)
                    implicit none
                    integer :: m, n
                    double precision, dimension(n), intent(in) :: x
                    double precision, dimension(m), intent(in) :: b
                    double precision, dimension(n) :: z
                end function
            end interface

            interface
                function dff(x, b, m, n) result (J)
                    implicit none
                    integer :: m, n
                    double precision, dimension(n), intent(in) :: x
                    double precision, dimension(m), intent(in) :: b
                    double precision, dimension(n, m) :: J
                end function
            end interface

            bk(:) = b0(:)

            do k=1, MAX_ITER
                J(:, :) = dff(x, b, m, n)
                db(:) = -MATMUL(inv(MATMUL(TRANSPOSE(J), J), n, ok), MATMUL(TRANSPOSE(J), ff(x, bk, m, n) - y))
                b(:) = bk(:) + db(:)

                if (.NOT. ok) then
                    call error("Tentativa de inversão de matriz singular.")
                    exit
                else if ((NORM(db, m) / NORM(b, m)) > TOL) then
                    bk(:) = b(:)
                else
                    return
                end if
            end do
            call error("O Método (Não-linear) de Mínimos Quadrados não convergiu.")
            return
        end function

        function sys_least_squares_num(ff, x, y, b0, m, n) result (b)
!           Same as previous function, with numerical partial derivatives
            implicit none
            integer :: m, n
            double precision, dimension(n), intent(in) :: x, y, b0
            double precision, dimension(n) :: b, bk, db, bh
            double precision :: J(n, n) 
            integer :: i, k
            logical :: ok = .TRUE.
            interface
                function ff(x, b, m, n) result (z)
                    implicit none
                    integer :: m, n
                    double precision, dimension(n), intent(in) :: x
                    double precision, dimension(m), intent(in) :: b
                    double precision, dimension(n) :: z
                end function
            end interface

            bh(:) = 0.0D0
            bk(:) = b0(:)

            do k=1, MAX_ITER
!               Compute Jacobian Matrix
                do i=1, m
!                   Partial derivative with respect do the i-th coordinates                    
                    bh(i) = h
                    J(:, i) = (ff(x, b(:) + bh(:), m, n) - ff(x, b(:) - bh(:), m, n)) / (2 * h)
                    bh(i) = 0.0D0
                end do
                db(:) = -MATMUL(inv(MATMUL(TRANSPOSE(J), J), n, ok), MATMUL(TRANSPOSE(J), ff(x, bk, m, n) - y))
                b(:) = bk(:) + db(:)

                if (.NOT. ok) then
                    call error("Tentativa de inversão de matriz singular.")
                    exit
                else if ((NORM(db, m) / NORM(b, m)) > TOL) then
                    bk(:) = b(:)
                else
                    return
                end if
            end do

            call error("O Método (Não-Linear) de Mínimos Quadrados não convergiu. (Derivada Numérica)")
            return
        end function

!       ============ Numerical Integration ========
        function polynomial_int(f, a, b) result (s)
            implicit none
            double precision :: a, b, s
            interface
                function f(x) result (y)
                    double precision :: x, y
                end function
            end interface


            
            return
        end function

        function gauss_quad_int(f, a, b) result (s)
            implicit none
            double precision :: a, b, s
            interface
                function f(x) result (y)
                    double precision :: x, y
                end function
            end interface

            return
        end function

        function gauss_hermite_quad_int(f, a, b) result (s)
            implicit none
            double precision :: a, b, s
            interface
                function f(x) result (y)
                    double precision :: x, y
                end function
            end interface

            return
        end function
    end module Calc
