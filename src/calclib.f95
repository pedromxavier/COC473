!   Calc Module

    module Calc
        use Util
        use Matrix
        implicit none
        double precision :: h = 1.0D-5

        character (len=*), parameter :: GAUSS_LEGENDRE_QUAD = "quadratures/gauss-legendre/gauss-legendre"
        character (len=*), parameter :: GAUSS_HERMITE_QUAD = "quadratures/gauss-hermite/gauss-hermite"
    contains
!       ================= Numerical Mathods =================
        function d(f, x, dx, kind) result (y)
            implicit none
            character (len=*), optional :: kind
            double precision, optional :: dx
            character (len=:), allocatable :: t_kind
            double precision :: x, y, t_dx
            
            interface
                function f(x) result (y)
                    implicit none
                    double precision :: x, y
                end function
            end interface

            if (.NOT. PRESENT(dx)) then
                t_dx = h
            else
                t_dx = dx
            end if

            if (.NOT. PRESENT(kind)) then
                t_kind = "central"
            else
                t_kind = kind
            end if

            if (t_kind == "central") then
                y = (f(x + t_dx) - f(x - t_dx)) / (2 * t_dx)
            else if (t_kind == "forward") then
                y = (f(x + t_dx) - f(x)) / t_dx
            else if (t_kind == "backward") then
                y = (f(x) - f(x - t_dx)) / t_dx
            else 
                call error("Unexpected value `"//t_kind//" for derivative kind."// &
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

        function newton(f, df, x0, ok) result (x)
            implicit none
            integer :: i
            double precision :: f, df
            double precision, intent(in) :: x0
            double precision :: x, xk
            logical :: ok

            ok = .TRUE.
            xk = x0

            do i = 1, MAX_ITER
                x = xk - f(xk) / df(xk)
                if (DABS(x - xk) > TOL) then
                    xk = x
                else
                    return
                end if 
            end do
            ok = .FALSE.
            return
        end function

        function secant(f, x0, ok) result (x)
            implicit none
            integer :: i
            double precision :: xk(3), yk(2)
            double precision, intent(in) :: x0
            double precision :: x
            logical :: ok
            interface
                function f(x) result (y)
                    implicit none
                    double precision :: x, y
                end function
            end interface

            ok = .TRUE.

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
            ok = .FALSE.
            return
        end function

        function inv_interp(f, x00, ok) result (x)
            implicit none
            double precision :: f
            double precision :: x, xk
            double precision, intent(in) :: x00(3)
            double precision :: x0(3), y0(3)
            integer :: i, k
            integer :: j(1)
            logical :: ok

            x0(:) = x00(:)
            xk = 1.0D+32

            ok = .TRUE.

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
            ok = .FALSE.
            return
        end function

        function sys_newton(ff, dff, xx0, n, ok) result (x)
            implicit none
            integer :: n
            double precision, dimension(n), intent(in) :: xx0
            double precision, dimension(n) :: x, xk, dx
            double precision :: J(n, n) 
            integer :: k
            logical :: ok
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

            ok = .TRUE.

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

        function sys_newton_num(ff, xx0, n, ok) result (x)
!           Same as previous function, with numerical partial derivatives            
            implicit none
            integer :: n
            double precision, dimension(n), intent(in) :: xx0
            double precision, dimension(n):: x, xk, xh, dx, y
            double precision, dimension(n, n) :: J
            integer :: i, k
            logical :: ok
            interface
                function ff(x, n) result (y)
                    implicit none
                    integer :: n
                    double precision :: x(n), y(n)
                end function
            end interface

            ok = .TRUE.

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

        function sys_broyden(ff, xx0, B0, n, ok) result (x)
            implicit none
            integer :: n
            double precision, dimension(n), intent(in) :: xx0
            double precision, dimension(n, n), intent(in) :: B0
            double precision, dimension(n) :: x, xk, yk, dx
            double precision, dimension(n, n) :: J, Bk
            integer :: k
            logical :: ok
            interface
                function ff(x, n) result (y)
                    implicit none
                    integer :: n
                    double precision :: x(n), y(n)
                end function
            end interface

            ok = .TRUE.

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

        function sys_least_squares(ff, dff, x, y, b0, m, n, ok) result (b)
            implicit none
            integer :: m, n
            double precision, dimension(n), intent(in) :: x, y, b0
            double precision, dimension(n) :: b, bk, db
            double precision :: J(n, n) 
            integer :: k
            logical :: ok
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

            ok = .TRUE.

            do k=1, MAX_ITER
                J(:, :) = dff(x, b, m, n)
                db(:) = -MATMUL(inv(MATMUL(TRANSPOSE(J), J), n, ok), MATMUL(TRANSPOSE(J), ff(x, bk, m, n) - y))
                b(:) = bk(:) + db(:)

                if (.NOT. ok) then
                    call error("Tentativa de inversão de matriz singular. k="//STR(k))
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

        function sys_least_squares_num(ff, x, y, b0, m, n, ok) result (b)
!           Same as previous function, with numerical partial derivatives
            implicit none
            integer :: m, n
            double precision, dimension(n), intent(in) :: x, y, b0
            double precision, dimension(n) :: b, bk, db, bh
            double precision :: J(n, n) 
            integer :: i, k
            logical :: ok
            interface
                function ff(x, b, m, n) result (z)
                    implicit none
                    integer :: m, n
                    double precision, dimension(n), intent(in) :: x
                    double precision, dimension(m), intent(in) :: b
                    double precision, dimension(n) :: z
                end function
            end interface

            ok = .TRUE.

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
        subroutine load_quad(x, w, k, fname)
!           Load Quadrature
            implicit none
            integer :: k, m, n
            character (len=*) :: fname
            double precision, dimension(k) :: x, w
            double precision, dimension(:, :), allocatable :: xw
            call read_matrix(fname, xw, m, n)
            if (n /= 2 .OR. m /= k) then
                call error("Invalid Matrix dimensions.")
                stop "ERROR"
            end if
            x(:) = xw(:, 1)
            w(:) = xw(:, 2)
            deallocate(xw)
        end subroutine

        function num_int(f, a, b, n, kind) result (s)
            implicit none
            integer :: n
            character (len=*), optional :: kind
            double precision :: a, b, s
            interface
                function f(x) result (y)
                    double precision :: x, y
                end function
            end interface

            if (.NOT. PRESENT(kind)) then
                kind = "polynomial"
            end if

            if (kind == "polynomial") then
                s = polynomial_int(f, a, b, n)
            else if (kind == "gauss-legendre") then
                s = gauss_legendre_int(f, a, b, n)
            else if (kind == "gauss-hermite") then
                s = gauss_hermite_int(f, a, b, n)
            else if (kind == "romberg") then
                s = romberg_int(f, a, b, n)
            else 
                call error("Unknown integration kind `"//kind//"."// &
                "Available options are: `polynomial`, `gauss-legendre`, `gauss-hermite` and `romberg`.")
            end if

        end function

        function polynomial_int(f, a, b, n) result (s)
            implicit none
            integer :: n, i
            double precision :: a, b, s
            double precision, dimension(n) :: x, y, w
            double precision, dimension(n, n) :: V
            interface
                function f(x) result (y)
                    double precision :: x, y
                end function
            end interface
            
            x(:) = ((b-a)/(n-1)) * (/ (i, i=0,n-1) /) + a
            y(:) = (/ ((b**i - a**i)/i, i=1, n) /)
            V(:, :) = vandermond_matrix(x, n)
            w(:) = solve(V, y, n)
            s = 0.0D0
            do i=1, n
                s = s + (w(i) * f(x(i)))
            end do
            return
        end function

        function gauss_legendre_int(f, a, b, n) result (s)
            implicit none
            integer, intent(in) :: n
            double precision, intent(in) :: a, b
            double precision :: s
            double precision, dimension(n) :: xx, ww
            integer :: k
            character(len=*), parameter :: fname = GAUSS_LEGENDRE_QUAD
            interface
                function f(x) result (y)
                    double precision :: x, y
                end function
            end interface

            call load_quad(xx, ww, n, fname//STR(n)//".txt")

            xx(:) = ((b - a) * xx(:) + (b + a)) / 2
            s = 0.0D0
            do k=1, n
                s = s + (ww(k) * f(xx(k)))
            end do
            s = s * ((b - a) / 2)
            return
        end function

        function gauss_hermite_int(f, a, b, n) result (s)
            implicit none
            integer, intent(in) :: n
            double precision, intent(in) :: a, b
            double precision :: s
            double precision, dimension(n) :: xx, ww
            integer :: k
            character(len=*), parameter :: fname = GAUSS_HERMITE_QUAD
            interface
                function f(x) result (y)
                    double precision :: x, y
                end function
            end interface

            call load_quad(xx, ww, n, fname//STR(n)//".txt")

            if (a /= DNINF .OR. b /= DINF) then
                call error("O Método de Gauss-Hermite deve ser usado no intervalo dos reais.")
                stop
            end if
            
            s = 0.0D0
            do k=1, n
                s = s + (ww(k) * f(xx(k)))
            end do

            return
        end function

        recursive function adapt_int(f, a, b, n, acc, kind) result (s)
            implicit none
            integer :: n
            character (len=*), optional :: kind
            double precision, intent(in) :: a, b
            double precision, optional :: acc
            interface
                function f(x) result (y)
                    double precision :: x, y
                end function
            end interface
            double precision :: p, q, e, r, s, t_tol

            if (.NOT. PRESENT(acc)) then
                t_tol = TOL
            else
                t_tol = acc
            end if

            p = num_int(f, a, b, n / 2, kind = kind)
            q = num_int(f, a, b, n, kind = kind)
            e = DABS(p - q)
            if (e <= t_tol) then
                s = q
            else
                r = (b + a) / 2
                s = adapt_int(f, a, r, n, acc = t_tol, kind = kind) + adapt_int(f, r, b, n, acc = t_tol, kind = kind)
            end if
            return
        end function

        function romberg_int(f, a, b, n, acc) result (s)
            implicit none
            integer, intent(in) :: n
            double precision, intent(in) :: a, b
            double precision, optional :: acc
            interface
                function f(x) result (y)
                    double precision :: x, y
                end function
            end interface
            integer :: i, j, k, t_n
            double precision :: p, r, s, dx, t_tol
!           Previous row, Current row and Temporary row            
            double precision, dimension(:), allocatable :: RP, RC, RT
            
            if (.NOT. PRESENT(acc)) then
                t_tol = TOL
            else
                t_tol = acc
            end if

            t_n = ILOG2(n)

            allocate(RP(t_n))
            allocate(RC(t_n))
            allocate(RT(t_n))

            dx = (b - a)

            RP(1) = (f(a) + f(b)) * (dx / 2)

            do i=2, t_n - 1
                dx = dx / 2

                r = 0
                k = 2 ** (i - 1)

                do j = 1, k
                    r = r + f(a + (2 * j - 1) * dx)
                end do

                RC(1) = (dx * r) + (RP(1) / 2)

                do j = 1, i
                    p = 4 ** j
                    RC(j + 1) = (p * RC(j) - RP(j)) / (p - 1)
                end do

                if ((i > 2) .AND. (DABS(RP(i - 1) - RC(i)) < acc)) then
                    RP(t_n) = RC(i - 1)
                    exit
                end if 

                RT(:) = RP(:)
                RP(:) = RC(:)
                RC(:) = RT(:)
            end do
            s = RP(t_n)

            deallocate(RP)
            deallocate(RC)
            deallocate(RT)
            return 
        end function

        function richard(f, x, p, q, dx, kind) result (y)
!           Richard Extrapolation
            implicit none
            double precision, optional :: dx
            character(len=*), optional :: kind
            double precision :: x, y, p, q, t_dx, dx1, dx2, d1, d2
            interface
                function f(x) result (y)
                    implicit none
                    double precision :: x, y
                end function
            end interface
            
            if (.NOT. PRESENT(dx)) then
                t_dx = h
            else
                t_dx = dx
            end if

            dx1 = t_dx
            d1 = d(f, x, dx1, kind = kind)
            dx2 = dx1 / q
            d2 = d(f, x, dx2, kind = kind)

            y = d1 + (d1 - d2) / ((q ** (-p)) - 1.0D0)
            return
        end function

!   ======== Ordinary Differential Equations ==========
    function ode_solve(df, y0, t, n, kind) result (y)
        implicit none
        integer :: n
        double precision, intent(in) :: y0
        double precision, dimension(n), intent(in) :: t
        double precision, dimension(n) :: y
        character(len=*), optional :: kind
        character(len=:), allocatable :: t_kind
        interface
            function df(t, y) result (u)
                implicit none
                double precision :: t, y, u
            end function
        end interface

        if (.NOT. PRESENT(kind)) then
            t_kind = 'euler'
        else
            t_kind = kind
        end if

        if (t_kind == 'euler') then
            y = euler(df, y0, t, n)
        else if (t_kind == 'runge-kutta2') then
            y = runge_kutta2(df, y0, t, n)
        else if (t_kind == 'runge-kutta4') then
            y = runge_kutta4(df, y0, t, n)
        else
            call error("As opções são: `euler`, `runge-kutta2` e `runge-kutta4`.")
            stop
        end if
        return
    end function


    function euler(df, y0, t, n) result (y)
        implicit none
        integer :: k, n
        double precision, intent(in) :: y0
        double precision :: dt
        double precision, dimension(n), intent(in) :: t
        double precision, dimension(n) :: y
        interface
            function df(t, y) result (u)
                implicit none
                double precision :: t, y, u
            end function
        end interface

        y(1) = y0
        do k=2, n
            dt = t(k) - t(k - 1)
            y(k) = y(k - 1) + df(t(k - 1), y(k - 1)) * dt
        end do
        return
    end function 

    function runge_kutta2(df, y0, t, n) result (y)
        implicit none
        integer :: k, n
        double precision, intent(in) :: y0
        double precision :: k1, k2, dt
        double precision, dimension(n), intent(in) :: t
        double precision, dimension(n) :: y
        interface
            function df(t, y) result (u)
                implicit none
                double precision :: t, y, u
            end function
        end interface

        y(1) = y0
        do k=2, n
            dt = t(k) - t(k - 1)
            k1 = df(t(k - 1), y(k - 1))
            k2 = df(t(k - 1) + dt, y(k - 1) + k1 * dt)
            y(k) = y(k - 1) + dt * (k1 + k2) / 2
        end do
        return
    end function 

    function runge_kutta4(df, y0, t, n) result (y)
        implicit none
        integer :: k, n
        double precision, intent(in) :: y0
        double precision :: k1, k2, k3, k4, dt
        double precision, dimension(n), intent(in) :: t
        double precision, dimension(n) :: y
        interface
            function df(t, y) result (u)
                implicit none
                double precision :: t, y, u
            end function
        end interface

        y(1) = y0
        do k=2, n
            dt = t(k) - t(k - 1)
            k1 = df(t(k - 1), y(k - 1))
            k2 = df(t(k - 1) + dt / 2, y(k - 1) + k1 * dt / 2)
            k3 = df(t(k - 1) + dt / 2, y(k - 1) + k2 * dt / 2)
            k4 = df(t(k - 1) + dt, y(k - 1) + dt * k3)
            y(k) = y(k - 1) + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6
        end do
        return
    end function

    function ode2_solve(d2f, y0, dy0, t, n, kind) result (y)
        implicit none
        integer :: n
        double precision, intent(in) :: y0, dy0
        double precision, dimension(n), intent(in) :: t
        double precision, dimension(n) :: y
        character(len=*), optional :: kind
        character(len=:), allocatable :: t_kind
        interface
            function d2f(t, y, dy) result (u)
                implicit none
                double precision :: t, y, dy, u
            end function
        end interface

        if (.NOT. PRESENT(kind)) then
            t_kind = 'taylor'
        else
            t_kind = kind
        end if

        if (t_kind == 'taylor') then
            y = taylor(d2f, y0, dy0, t, n)
        else if (t_kind == 'runge-kutta-nystrom') then
            y = runge_kutta_nystrom(d2f, y0, dy0, t, n)
        else
            call error("As opções são: `taylor`, `runge-kutta-nystrom`.")
            stop
        end if
        return
    end function

    function taylor(d2f, y0, dy0, t, n) result (y)
        implicit none
        integer :: k, n
        double precision, intent(in) :: y0, dy0
        double precision :: dt, dy, d2y
        double precision, dimension(n), intent(in) :: t
        double precision, dimension(n) :: y
        interface
            function d2f(t, y, dy) result (d2y)
                implicit none
                double precision :: t, y, dy, d2y
            end function
        end interface
!       Solution
        y(1) = y0
!       1st derivative
        dy = dy0
        do k=2, n
            dt = t(k) - t(k - 1)
            d2y = d2f(t(k - 1), y(k - 1), dy)
            y(k) = y(k - 1) + (dy * dt) + (d2y * dt ** 2) / 2
            dy = dy + d2y * dt
        end do
        return
    end function

    function runge_kutta_nystrom(d2f, y0, dy0, t, n) result (y)
        implicit none
        integer :: k, n
        double precision, intent(in) :: y0, dy0
        double precision :: k1, k2, k3, k4, dt, dy, l, q
        double precision, dimension(n), intent(in) :: t
        double precision, dimension(n) :: y
        interface
            function d2f(t, y, dy) result (u)
                implicit none
                double precision :: t, y, dy, u
            end function
        end interface

        y(1) = y0
        dy = dy0
        do k=2, n
            dt = t(k) - t(k - 1)
            k1 = (d2f(t(k - 1), y(k - 1), dy) * dt) / 2
            q = ((dy + k1 / 2) * dt) / 2
            k2 = (d2f(t(k - 1) + dt / 2, y(k - 1) + q, dy + k1) * dt) / 2
            k3 = (d2f(t(k - 1) + dt / 2, y(k - 1) + q, dy + k2) * dt) / 2
            l = (dy + k3) * dt
            k4 = (d2f(t(k - 1) + dt, y(k - 1) + l, dy + 2* k3) * dt) / 2

            y(k) = y(k - 1) + (dy + (k1 + k2 + k3) / 3) * dt
            dy = dy + (k1 + 2 * k2 + 2 * k3 + k4) / 3
        end do
        return
    end function
    end module Calc
