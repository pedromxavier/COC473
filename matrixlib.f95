!   Matrix Module

    module Matrix
        implicit none
        integer :: NMAX = 1000

    contains

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

        subroutine ill_cond()
!           Prompts the user with an ill-conditioning warning.
            implicit none
            call error('Matriz mal-condicionada.')
        end subroutine

        subroutine print_matrix(A, m, n)
            implicit none

            integer :: m, n
            double precision :: A(m, n)

            integer :: i, j

20          format(' |', F10.5, ' ')
21          format(F10.5, '|')
22          format(F10.5, ' ')

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

30          format(' |', F10.5, '|')

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

!       =========== Matrix Methods ============
        recursive function det(A, n) result (d)
            implicit none

            integer :: n
            double precision :: A(n, n)
            double precision :: X(n-1, n-1)

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

        function rand_vector(n) result (x)
            implicit none
            integer :: n
            double precision :: x (n)

            integer :: i

            do i = 1, n
                x(i) = 2 * ran(0) - 1
            end do
            return
        end function

        function rand_matrix(m, n) result (A)
            implicit none
            integer :: m, n
            double precision :: A(m, n)

            integer :: i

            do i = 1, m
                A(i, :) = rand_vector(n)
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

        recursive function positive_definite(A, n) result (x)
!       Checks wether a matrix is positive definite
!       according to Sylvester's criterion.
            implicit none

            integer :: n
            double precision A(n, n)

            logical :: x

            if (n == 1) then
                x = (A(1, 1) > 0)
                return
            else
                x = positive_definite(A(:n-1, :n-1), n-1) .AND. (det(A, n) > 0)
                return
            end if
        end function

        function symmetrical(A, n) result (x)
!           Check if the Matrix is symmetrical
            integer :: n

            double precision :: A(n, n)

            integer :: i, j
            logical :: x
            
            do i = 1, n
                do j = 1, i-1
                    if (A(i, j) /= A(j, i)) then
                        x = .FALSE.
                        return
                    end if
                end do
            end do
            x = .TRUE.
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

        function matrix_norm(A, n) result (s)
!           Frobenius norm
            implicit none
            integer :: n
            double precision :: A(n, n)
            double precision :: s

            s = sqrt(sum(A * A))
            return 
        end function

        function spectral_radius(A, n, k) result (r)
            implicit none

            integer :: n
            double precision :: A(n, n), M(n, n)
            double precision :: r

            integer :: i, k

            M(:, :) = A(:, :)

            do i = 1, k
                M = matmul(M, M)
            end do
            r = matrix_norm(M, n)
            do i = 1, k
                r = sqrt(r)
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
        function Cholesky_cond(A, n) result (x)
            implicit none

            integer :: n
            double precision :: A(n, n)

            logical :: x

            x = symmetrical(A, n) .AND. positive_definite(A, n)
            return

        end function

        function PLU_cond(A, n) result (x)
            implicit none

            integer :: n
            double precision A(n, n)

            integer :: i, j
            double precision :: s

            logical :: x

            do j = 1, n
                s = 0.0D0
                do i = 1, j
                    if (A(i, j) > s) then
                        s = A(i, j)
                    end if
                end do
            end do
            
            x = (s < 0.01D0)

            return
        end function

        function LU_cond(A, n) result (x)
            implicit none

            integer :: n
            double precision A(n, n)

            logical :: x

            x = positive_definite(A, n)

            return
        end function

!       ======= Matrix Factorization Methods ========
        function PLU_decomp(A, P, L, U, n) result (x)
            implicit none

            integer :: n
            double precision :: A(n,n), P(n,n), L(n,n), U(n,n)

            logical :: x

!           Permutation Matrix
            P = pivot_matrix(A, n)

!           Decomposition over Row-Swapped Matrix
            x = LU_decomp(matmul(P, A), L, U, n)
            return
        end function
            
        function LU_decomp(A, L, U, n) result (x)
            implicit none

            integer :: n
            double precision :: A(n, n), L(n, n), U(n,n), M(n, n)

            logical :: x

            integer :: i, j, k

!           Results Matrix
            M(:, :) = A(:, :)

            if (.NOT. LU_cond(A, n)) then
                call ill_cond()
                x = .FALSE.
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

            call LU_matrix(M, L, U, n)

            x = .TRUE.
            return

        end function

        function Cholesky_decomp(A, L, n) result (x)
            implicit none

            integer :: n
            double precision :: A(n, n), L(n, n)

            logical :: x

            integer :: i, j

            if (.NOT. Cholesky_cond(A, n)) then
                call ill_cond()
                x = .FALSE.
                return
            end if

            do i = 1, n
                L(i, i) = sqrt(A(i, i) - sum(L(i, :i-1) * L(i, :i-1)))
                do j = 1 + 1, n
                    L(j, i) = (A(i, j) - sum(L(i, :i-1) * L(j, :i-1))) / L(i, i)
                end do
            end do

            x = .TRUE.
            return

        end function

!       === Linear System Solving Conditions ===
        function Jacobi_cond(A, n) result (x)
            implicit none

            integer :: n

            double precision :: A(n, n)

            logical :: x

            x = (spectral_radius(A, n, 1000) < 1)

            return
        end function

!       == Linear System Solving Methods ==
        function Jacobi(A, x, b, e, n) result (ok)
            implicit none
            
            integer :: n

            double precision :: A(n, n)
            double precision :: b(n), x(n), x0(n)

            logical :: ok

            double precision :: e

            integer :: i, k

            x0 = rand_vector(n)

            ok = Jacobi_cond(A, n)

            if (.NOT. ok) then
                call ill_cond()
                return
            end if

            do k = 1, NMAX
                do i = 1, n
                    x(i) = (b(i) - dot_product(A(i, :), x0)) / A(i, i)
                end do
                x0(:) = x(:)

                e = vector_norm(x - b, n)
            end do
            return
        end function

    end module Matrix
