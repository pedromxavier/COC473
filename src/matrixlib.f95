!   Matrix Module

    module Matrix
        use Util
        implicit none
        integer :: D_MAX_ITER = 1000
        double precision :: D_TOL = 1.0D-5
    contains
        subroutine ill_cond()
!           Prompts the user with an ill-conditioning warning.
            implicit none
            call error('Matriz mal-condicionada: este método não irá convergir.')
        end subroutine

        subroutine show_matrix(var, A, m, n)
            implicit none
            integer :: m, n
            character(len=*) :: var
            double precision, dimension(m, n), intent(in) :: A
            write (*, *) ''//achar(27)//'[36m'//var//' = '
            call print_matrix(A, m, n)
            write (*, *) ''//achar(27)//'[0m'
        end subroutine

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
            double precision, dimension(:, :), allocatable :: A
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

        subroutine read_vector(fname, b, n)
            implicit none
            character(len=*) :: fname
            integer :: n
            double precision, allocatable :: b(:)

            open(unit=33, file=fname, status='old', action='read')
            read(33, *) n
            allocate(b(n))
            read(33, *) b(:)
            close(33)
        end subroutine

        subroutine show_vector(var, x, n)
            implicit none
            integer :: n
            character(len=*) :: var
            double precision :: x(n)
            write (*, *) ''//achar(27)//'[36m'//var//' = '
            call print_vector(x, n)
            write (*, *) ''//achar(27)//'[0m'
        end subroutine


!       =========== Matrix Methods ============

        function clip(x, n, a, b) result (y)
            integer, intent(in) :: n
            integer :: k
            double precision, intent(in) :: a, b
            double precision, dimension(n), intent(in) :: x
            double precision, dimension(n) :: y

            do k=1, n
                if ((a <= x(k)) .AND. (x(k) <= b)) then
                    y(k) = x(k)
                else
                    y(k) = DNAN
                end if
            end do
            return
        end function

        function rand_vector(n, a, b) result (r)
            implicit none
            integer :: n, i
            double precision, dimension(n) :: r
            double precision, optional :: a, b
            double precision :: t_a, t_b
            
            if (.NOT. PRESENT(a)) then
                t_a = -1.0D0
            else
                t_a = a
            end if

            if (.NOT. PRESENT(b)) then
                t_b = 1.0D0
            else
                t_b = b
            end if

            do i = 1, n
                r(i) = DRAND(t_a, t_b)
            end do
            return
        end function

        function rand_matrix(m, n, a, b) result (R)
            implicit none
            integer :: m, n, i
            double precision, dimension(m, n) :: R
            double precision, optional :: a, b
            
            do i = 1, m
                R(i, :) = rand_vector(n, a=a, b=b)
            end do
            return
        end function

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

        function given_matrix(A, n, i, j) result (G)
            implicit none

            integer :: n, i, j
            double precision :: A(n, n), G(n, n)
            double precision :: t, c, s

            G(:, :) = id_matrix(n)

            t = 0.5D0 * DATAN2(2.0D0 * A(i,j), A(i, i) - A(j, j))
            s = DSIN(t)
            c = DCOS(t)

            G(i, i) = c
            G(j, j) = c
            G(i, j) = -s
            G(j, i) = s

            return
        end function

        function vandermond_matrix(x, n) result (V)
            implicit none
            integer :: n, i
            double precision, dimension(n), intent(in) :: x
            double precision, dimension(n, n) :: V
            V(1, :) = 1.0D0
            do i=2, n
                V(i, :) = V(i-1, :) * x(:)
            end do
            return
        end function

        function diagonally_dominant(A, n) result (ok)
            implicit none

            integer :: n
            double precision :: A(n, n)

            logical :: ok
            integer :: i

            do i = 1, n
                if (DABS(A(i, i)) < SUM(DABS(A(i, :i-1))) + SUM(DABS(A(i, i+1:)))) then
                    ok = .FALSE.
                    return
                end if
            end do
            ok = .TRUE.
            return
        end function

        recursive function positive_definite(A, n) result (ok)
!       Checks wether a matrix is positive definite
!       according to Sylvester's criterion.
            implicit none

            integer :: n
            double precision A(n, n)

            logical :: ok

            if (n == 1) then
                ok = (A(1, 1) > 0)
                return
            else
                ok = positive_definite(A(:n-1, :n-1), n-1) .AND. (det(A, n) > 0)
                return
            end if
        end function

        function symmetrical(A, n) result (ok)
!           Check if the Matrix is symmetrical
            integer :: n

            double precision :: A(n, n)

            integer :: i, j
            logical :: ok
            
            do i = 1, n
                do j = 1, i-1
                    if (A(i, j) /= A(j, i)) then
                        ok = .FALSE.
                        return
                    end if
                end do
            end do
            ok = .TRUE.
            return
        end function

        subroutine swap_rows(A, i, j, n)
            implicit none

            integer :: n
            integer :: i, j
            double precision A(n, n)
            double precision temp(n)

            temp(:) = A(i, :)
            A(i, :) = A(j, :)
            A(j, :) = temp(:)
        end subroutine

        function outer_product(x, y, n) result (A)
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

!       ================= Matrix Method ====================
        function inv(A, n, ok) result (Ainv)
            integer :: n
            double precision :: A(n, n), Ainv(n, n)
            double precision :: work(n)
            integer :: ipiv(n)   ! pivot indices
            integer :: info

            logical :: ok
          
            ! External procedures defined in LAPACK
            external DGETRF
            external DGETRI
          
            ! Store A in Ainv to prevent it from being overwritten by LAPACK
            Ainv(:, :) = A(:, :)
          
            ! DGETRF computes an LU factorization of a general M-by-N matrix A
            ! using partial pivoting with row interchanges.
            call DGETRF(n, n, Ainv, n, ipiv, info)
          
            if (info /= 0) then
                ok = .FALSE.
                return
            end if
          
            ! DGETRI computes the inverse of a matrix using the LU factorization
            ! computed by DGETRF.
            call DGETRI(n, Ainv, n, ipiv, work, n, info)
          
            if (info /= 0) then
                ok = .FALSE.
                return
            end if

            return
        end function

        function row_max(A, j, n) result(k)
            implicit none

            integer :: n
            double precision A(n, n)
            
            integer :: i, j, k
            double precision :: s

            s = 0.0D0
            do i = j, n
                if (A(i, j) > s) then
                    s = A(i, j)
                    k = i
                end if
            end do 
            return 
        end function

        function pivot_matrix(A, n) result (P)
            implicit none

            integer :: n
            double precision :: A(n, n)

            double precision :: P(n, n)

            integer :: j, k

            P = id_matrix(n)

            do j = 1, n
                k = row_max(A, j, n)
                if (j /= k) then 
                    call swap_rows(P, j, k, n)
                end if
            end do
            return
        end function

        function vector_norm(x, n) result (s)
            implicit none
            integer :: n
            double precision :: x(n)
            double precision :: s
            s = sqrt(dot_product(x, x))
            return
        end function

        function NORM(x, n) result (s)
            implicit none
            integer :: n
            double precision :: x(n)
            double precision :: s
            s = SQRT(DOT_PRODUCT(x, x))
            return
        end function

        function matrix_norm(A, n) result (s)
!           Frobenius norm
            implicit none
            integer :: n
            double precision :: A(n, n)
            double precision :: s

            s = DSQRT(SUM(A * A))
            return 
        end function

        function spectral_radius(A, n) result (r)
            implicit none

            integer :: n
            double precision :: A(n, n), x(n)
            double precision :: r, l
            logical :: ok        
            ok = power_method(A, n, x, l)
            r = DABS(l)    
            return
        end function

        recursive function det(A, n) result (d)
            implicit none
            integer :: n
            double precision, dimension(n, n) :: A
            double precision, dimension(n-1, n-1) :: X
            integer :: i
            double precision :: d, s

            if (n == 1) then
                d = A(1, 1)
                return
            elseif (n == 2) then
                d = A(1, 1) * A(2, 2) - A(1, 2) * A(2, 1)
                return
            else
                d = 0.0D0
                s = 1.0D0
                do i = 1, n
!                   Compute submatrix X
                    X(:,  :i-1) = A(2:,    :i-1)
                    X(:, i:   ) = A(2:, i+1:   )
                    d = s * det(X, n-1) * A(1, i) + d
                    s = -s
                end do
            end if
            return
        end function

        function LU_det(A, n) result (d)
            implicit none

            integer :: n
            integer :: i
            double precision :: A(n, n), L(n, n), U(n, n)
            double precision :: d

            d = 0.0D0

            if (.NOT. LU_decomp(A, L, U, n)) then
                call ill_cond()
                return
            end if

            do i = 1, n
                d = d * L(i, i) * U(i, i)
            end do

            return
        end function

        subroutine LU_matrix(A, L, U, n)
!           Splits Matrix in Lower and Upper-Triangular
            implicit none

            integer :: n
            double precision :: A(n, n), L(n, n), U(n, n)

            integer :: i

            L(:, :) = 0.0D0
            U(:, :) = 0.0D0
 
            do i = 1, n
                L(i, i) = 1.0D0
                L(i,  :i-1) = A(i,  :i-1)
                U(i, i:   ) = A(i, i:   )
            end do
        end subroutine

!       === Matrix Factorization Conditions ===
        function Cholesky_cond(A, n) result (ok)
            implicit none
            integer :: n
            double precision :: A(n, n)
            logical :: ok
            ok = symmetrical(A, n) .AND. positive_definite(A, n)
            return
        end function

        function PLU_cond(A, n) result (ok)
            implicit none
            integer :: n
            double precision A(n, n)
            integer :: i, j
            double precision :: s
            logical :: ok
            do j = 1, n
                s = 0.0D0
                do i = 1, j
                    if (A(i, j) > s) then
                        s = A(i, j)
                    end if
                end do
            end do
            ok = (s < 0.01D0)
            return
        end function

        function LU_cond(A, n) result (ok)
            implicit none
            integer :: n
            double precision A(n, n)
            logical :: ok
            ok = positive_definite(A, n)
            return
        end function
!        _      _____  _____ _______         __ 
!       | |    |_   _|/ ____|__   __|/\     /_ |
!       | |      | | | (___    | |  /  \     | |
!       | |      | |  \___ \   | | / /\ \    | |
!       | |____ _| |_ ____) |  | |/ ____ \   | |
!       |______|_____|_____/   |_/_/    \_\  |_|
!       ========================================
                                                
!       ======= Matrix Factorization Methods ========
        function PLU_decomp(A, P, L, U, n) result (ok)
            implicit none
            integer :: n
            double precision :: A(n,n), P(n,n), L(n,n), U(n,n)
            logical :: ok
!           Permutation Matrix
            P = pivot_matrix(A, n)
!           Decomposition over Row-Swapped Matrix
            ok = LU_decomp(matmul(P, A), L, U, n)
            return
        end function
            
        function LU_decomp(A, L, U, n) result (ok)
            implicit none
            integer :: n
            double precision :: A(n, n), L(n, n), U(n,n), M(n, n)
            logical :: ok
            integer :: i, j, k
!           Results Matrix
            M(:, :) = A(:, :)
            if (.NOT. LU_cond(A, n)) then
                call ill_cond()
                ok = .FALSE.
                return
            end if
            do k = 1, n-1
                do i = k+1, n
                    M(i, k) = M(i, k) / M(k, k)
                end do
                do j = k+1, n
                    do i = k+1, n
                        M(i, j) = M(i, j) - M(i, k) * M(k, j)
                    end do
                end do
            end do

!           Splits M into L & U
            call LU_matrix(M, L, U, n)

            ok = .TRUE.
            return

        end function

        function Cholesky_decomp(A, L, n) result (ok)
            implicit none

            integer :: n
            double precision :: A(n, n), L(n, n)

            logical :: ok

            integer :: i, j

            if (.NOT. Cholesky_cond(A, n)) then
                call ill_cond()
                ok = .FALSE.
                return
            end if

            do i = 1, n
                L(i, i) = sqrt(A(i, i) - sum(L(i, :i-1) * L(i, :i-1)))
                do j = 1 + 1, n
                    L(j, i) = (A(i, j) - sum(L(i, :i-1) * L(j, :i-1))) / L(i, i)
                end do
            end do

            ok = .TRUE.
            return
        end function

        function Jacobi_cond(A, n) result (ok)
            implicit none

            integer :: n

            double precision :: A(n, n)

            logical :: ok

            if (.NOT. spectral_radius(A, n) < 1.0D0) then
                ok = .FALSE.
                call ill_cond()
                return
            else
                ok = .TRUE.
                return
            end if
        end function

        function Jacobi(A, x, b, e, n, tol, max_iter) result (ok)
            implicit none
            
            logical :: ok

            integer :: n, i, k, t_max_iter
            integer, optional :: max_iter

            double precision :: A(n, n)
            double precision :: b(n), x(n), x0(n)
            double precision :: e, t_tol
            double precision, optional :: tol

            if (.NOT. PRESENT(tol)) then
                t_tol = D_TOL
            else
                t_tol = tol
            end if

            if (.NOT. PRESENT(max_iter)) then
                t_max_iter = D_MAX_ITER
            else
                t_max_iter = max_iter
            end if

            x0 = rand_vector(n)

            ok = Jacobi_cond(A, n)

            if (.NOT. ok) then
                return
            end if

            do k = 1, t_max_iter
                do i = 1, n
                    x(i) = (b(i) - dot_product(A(i, :), x0)) / A(i, i)
                end do
                x0(:) = x(:)
                e = vector_norm(matmul(A, x) - b, n)
                if (e < t_tol) then
                    return
                end if
            end do
            call error('Erro: Esse método não convergiu.')
            ok = .FALSE.
            return
        end function

        function Gauss_Seidel_cond(A, n) result (ok)
            implicit none

            integer :: n

            double precision :: A(n, n)

            logical :: ok

            integer :: i

            do i = 1, n
                if (A(i, i) == 0.0D0) then
                    ok = .FALSE.
                    call ill_cond()
                    return
                end if
            end do

            if (symmetrical(A, n) .AND. positive_definite(A, n)) then
                ok = .TRUE.
                return
            else
                call warn('Aviso: Esse método pode não convergir.')
                return
            end if
        end function

        function Gauss_Seidel(A, x, b, e, n, tol, max_iter) result (ok)
            implicit none
            logical :: ok
            integer :: n, i, j, k, t_max_iter
            integer, optional :: max_iter
            double precision :: A(n, n)
            double precision :: b(n), x(n)
            double precision :: e, s, t_tol
            double precision, optional :: tol

            if (.NOT. PRESENT(tol)) then
                t_tol = D_TOL
            else
                t_tol = tol
            end if

            if (.NOT. PRESENT(max_iter)) then
                t_max_iter = D_MAX_ITER
            else
                t_max_iter = max_iter
            end if

            ok = Gauss_Seidel_cond(A, n)

            if (.NOT. ok) then
                return
            end if

            do k = 1, t_max_iter
                do i = 1, n
                    s = 0.0D0
                    do j = 1, n
                        if (i /= j) then
                            s = s + A(i, j) * x(j)
                        end if
                    end do
                    x(i) = (b(i) - s) / A(i, i)
                end do
                e = vector_norm(matmul(A, x) - b, n)
                if (e < t_tol) then
                    return
                end if
            end do
            call error('Erro: Esse método não convergiu.')
            ok = .FALSE.
            return
        end function

!       Decomposição LU e afins
        subroutine LU_backsub(L, U, x, y, b, n)
            implicit none
            integer :: n
            double precision :: L(n, n), U(n, n)
            double precision :: b(n), x(n), y(n)
            integer :: i
!           Ly = b (Forward Substitution)
            do i = 1, n
                y(i) = (b(i) - SUM(L(i, 1:i-1) * y(1:i-1))) / L(i, i)
            end do
!           Ux = y (Backsubstitution)
            do i = n, 1, -1
                x(i) = (y(i) - SUM(U(i,i+1:n) * x(i+1:n))) / U(i, i)
            end do
        end subroutine

        function LU_solve(A, x, y, b, n) result (ok)
            implicit none

            integer :: n

            double precision :: A(n, n), L(n, n), U(n, n)
            double precision :: b(n), x(n), y(n)

            logical :: ok

            ok = LU_decomp(A, L, U, n)

            if (.NOT. ok) then
                return
            end if

            call LU_backsub(L, U, x, y, b, n)

            return
        end function

        function PLU_solve(A, x, y, b, n) result (ok)
            implicit none

            integer :: n

            double precision :: A(n, n), P(n,n), L(n, n), U(n, n)
            double precision :: b(n), x(n), y(n)

            logical :: ok

            ok = PLU_decomp(A, P, L, U, n)

            if (.NOT. ok) then
                return
            end if

            call LU_backsub(L, U, x, y, matmul(P, b), n)

            x(:) = matmul(P, x)

            return
        end function

        function Cholesky_solve(A, x, y, b, n) result (ok)
            implicit none

            integer :: n

            double precision :: A(n, n), L(n, n), U(n, n)
            double precision :: b(n), x(n), y(n)

            logical :: ok

            ok = Cholesky_decomp(A, L, n)

            if (.NOT. ok) then
                return
            end if

            U = transpose(L)

            call LU_backsub(L, U, x, y, b, n)

            return
        end function

!        _      _____  _____ _______         ___  
!       | |    |_   _|/ ____|__   __|/\     |__ \
!       | |      | | | (___    | |  /  \       ) |
!       | |      | |  \___ \   | | / /\ \     / / 
!       | |____ _| |_ ____) |  | |/ ____ \   / /_  
!       |______|_____|_____/   |_/_/    \_\ |____|
!       ==========================================

!       ============ Power Method ============
        function power_method(A, n, x, l, tol, max_iter) result (ok)
            implicit none
            logical :: ok
            integer :: n, k, t_max_iter
            integer, optional :: max_iter
            double precision :: A(n, n)
            double precision :: x(n)
            double precision :: l, ll, t_tol
            double precision, optional :: tol

            if (.NOT. PRESENT(tol)) then
                t_tol = D_TOL
            else
                t_tol = tol
            end if

            if (.NOT. PRESENT(max_iter)) then
                t_max_iter = D_MAX_ITER
            else
                t_max_iter = max_iter
            end if

!           Begin with random normal vector and set 1st component to zero
            x(:) = rand_vector(n)
            x(1) = 1.0D0

!           Initialize Eigenvalues
            l = 0.0D0

!           Checks if error tolerance was reached          
            do k=1, t_max_iter
                ll = l

                x(:) = matmul(A, x)                

!               Retrieve Eigenvalue
                l = x(1)

!               Retrieve Eigenvector                
                x(:) = x(:) / l

                if (dabs((l - ll) / l) < t_tol) then
                    ok = .TRUE.
                    return
                end if
            end do
            ok = .FALSE.
            return
        end function
        
        function Jacobi_eigen(A, n, L, X, tol, max_iter) result (ok)
            implicit none
            logical :: ok
            integer :: n, i, j, k, u, v, t_max_iter
            integer, optional :: max_iter
            double precision :: A(n, n), L(n, n), X(n, n), P(n, n)
            double precision :: y, z, t_tol
            double precision, optional :: tol

            if (.NOT. PRESENT(tol)) then
                t_tol = D_TOL
            else
                t_tol = tol
            end if

            if (.NOT. PRESENT(max_iter)) then
                t_max_iter = D_MAX_ITER
            else
                t_max_iter = max_iter
            end if

            X(:, :) = id_matrix(n)
            L(:, :) = A(:, :)

            do k=1, t_max_iter
                z = 0.0D0
                do i = 1, n
                    do j = 1, i - 1
                        y = DABS(L(i, j))

!                       Found new maximum absolute value                        
                        if (y > z) then
                            u = i
                            v = j
                            z = y
                        end if
                    end do
                end do

                if (z >= t_tol) then
                    P(:, :) = given_matrix(L, n, u, v)
                    L(:, :) = matmul(matmul(transpose(P), L), P)
                    X(:, :) = matmul(X, P)
                else
                    ok = .TRUE.
                    return
                end if
            end do
            ok = .FALSE.
            return
        end function

!        _      _____  _____ _______         _____  
!       | |    |_   _|/ ____|__   __|/\     |__   |
!       | |      | | | (___    | |  /  \       )  /
!       | |      | |  \___ \   | | / /\ \     |_  \
!       | |____ _| |_ ____) |  | |/ ____ \   ___)  |
!       |______|_____|_____/   |_/_/    \_\ |_____/ 
!       ==========================================

        function least_squares(x, y, s, n) result (ok)
            implicit none
            integer :: n

            logical :: ok

            double precision :: A(2,2), b(2), s(2), r(2), x(n), y(n)
            
            A(1, 1) = n
            A(1, 2) = SUM(x)
            A(2, 1) = SUM(x)
            A(2, 2) = dot_product(x, x)

            b(1) = SUM(y)
            b(2) = dot_product(x, y)

            ok = Cholesky_solve(A, s, r, b, n)
            return
        end function

!       ========== Extra Stuff ========     
        
        function Gauss_solve(A0, x, b0, n) result (ok)
            implicit none 
            integer n
            double precision, dimension(n, n), intent(in) :: A0
            double precision, dimension(n, n) :: A
            double precision, dimension(n), intent(in) :: b0
            double precision, dimension(n) :: b, x, s
            double precision :: c, pivot, store
            integer i, j, k, l

            logical :: ok
            
            ok = .TRUE.
            
            A(:, :) = A0(:, :)
            b(:) = b0(:)

            do k=1, n-1
                do i=k,n
                    s(i) = 0.0
                    do j=k,n
                        s(i) = MAX(s(i), DABS(A(i,j)))
                    end do
                end do
            
                pivot = DABS(A(k,k) / s(k))
                l = k
                do j=k+1,n
                    if(DABS(A(j,k) / s(j)) > pivot) then
                        pivot = DABS(A(j,k) / s(j))
                        l = j
                    end if
                end do
            
                if(pivot == 0.0) then
                    ok = .FALSE.
                    return
                end if
            
                if (l /= k) then
                    do j=k,n
                        store = A(k,j)
                        A(k,j) = A(l,j)
                        A(l,j) = store
                    end do
                    store = b(k)
                    b(k) = b(l)
                    b(l) = store
                end if
            
                do i=k+1,n
                    c = A(i,k) / A(k,k)
                    A(i,k) = 0.0D0
                    b(i) = b(i)- c*b(k)
                    do j=k+1,n
                        A(i,j) = A(i,j) - c * A(k,j)
                    end do
                end do
            end do
            
            x(n) = b(n) / A(n,n)
            do i=n-1,1,-1
                c = 0.0D0
                do j=i+1,n
                    c = c + A(i,j) * x(j)
                end do 
                x(i) = (b(i)- c) / A(i,i)
            end do
            
            return
        end function
            
        function solve(A, b, n, kind) result (x)
            implicit none
            integer :: n
            double precision, dimension(n), intent(in) :: b
            double precision, dimension(n) :: x, y
            double precision, dimension(n, n), intent(in) :: A
            character(len=*), optional :: kind
            character(len=:), allocatable :: t_kind

            logical :: ok = .TRUE.

            if (.NOT. PRESENT(kind)) then
                call debug("Indeed, not present.")
                t_kind = "gauss"
            else
                t_kind = kind
            end if

            call debug("Now it is: "//t_kind)
            if (t_kind == "LU") then
                ok = LU_solve(A, x, y, b, n)
            else if (t_kind == "PLU") then
                ok = PLU_solve(A, x, y, b, n)
            else if (t_kind =="cholesky") then
                ok = Cholesky_solve(A, x, y, b, n)
            else if (t_kind =="gauss") then
                ok = Gauss_solve(A, x, b, n)
            else
                ok = .FALSE.
            end if

            call debug(":: Solved via `"//t_kind//"` ::")

            if (.NOT. ok) then
                call error("Failed to solve system Ax = b.")
            end if

            return
        end function
        
    end module Matrix
