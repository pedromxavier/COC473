!   Matrix Module

    module Matrix
        implicit none

    contains

        function ID_MATRIX(n) result (A)
            implicit none

            integer :: n
            double precision :: A(n, n)

            integer :: j

            A(:, :) = 0.0D0

            do j = 1, n
                A(j, j) = 1.0D0
            end do

        end function

        subroutine SWAP_ROWS(A, i, j, n)
            implicit none

            integer :: n
            integer :: i, j
            double precision A(n, n)
            double precision temp(n)

            temp(:) = A(i, :)
            A(i, :) = A(j, :)
            A(j, :) = temp(:)

        end subroutine 

        function ROW_MAX(A, j, n) result(k)
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

        end function

        function COND(A, n) result (x)
            implicit none
            
            integer :: n
            double precision :: A(n, n)
            double precision :: x


        end function

        function PIVOT_MATRIX(A, n) result (M)
            implicit none

            integer :: n
            double precision :: A(n, n)
            double precision :: M(n, n)

            integer :: i, j, k

            M = ID_MATRIX(n)

            do j = 1, n
                k = ROW_MAX(A, j, n)
                if (j /= k) then 
                    call SWAP_ROWS(M, i, j, n)
                end if
            end do

            return
        end function
            
        subroutine LU_DECOMP(A, L, U, n)
            implicit none

            integer :: n
            double precision :: A(n, n), L(n, n), U(n, n)
            double precision :: s

            integer :: i, j, k

!           Lower-Triangular Matrix
            L(:, :) = 0.0D0

!           Upper-Triangular Matrix
            U(:, :) = 0.0D0

            do i = 1, n
!               Upper-Triangular
                do k = 1, n
                    s = 0.0D0
                    do j = 1, i-1
                        s = s + L(i, j) * U(j, k)
                    end do
                    U(i, k) = A(i, k) - s
                end do
!               Lower-Triangular
                do k = 1, n
                    if (i == k) then
                        L(i, i) = 1.0D0
                    else
                        s = 0.0D0
                        do j = 1, i-1
                            s = s + L(k, j) * U(j, i)
                        end do
                        L(k, i) = (A(k, i) - s) / U(i, i)
                    end if
                end do
            end do

        end subroutine

        subroutine CHOLESKY_DECOMP(A, n)
            implicit none

            integer :: n
            double precision :: A(n, n)

        end subroutine

        subroutine print_matrix(A, m, n)
            implicit none

            integer :: m, n
            double precision :: A(m, n)

            integer :: i, j

20          format('|', F6.2, ' ')
21          format(F6.2, '|')
22          format(F6.2, ' ')

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

        subroutine print_vector(x, n)
            implicit none

            integer :: n
            double precision :: x(n)
            
            integer :: i

30          format('|', F6.2, '|')

            do i = 1, n
                write(*, 30) x(i)
            end do
        end subroutine

        recursive function det(A, n) result (d)
            implicit none

            integer :: n
            double precision :: A(n, n)
            double precision :: X(n-1, n-1)

            integer :: i, j, k
            double precision :: d, s

            if (n == 2) then
                d = A(1, 1) * A(2, 2) - A(1, 2) * A(2, 1)
            else
                d = 0.0D0
                s = 1.0D0
                do i = 1, n
                    X(:,  :i-1) = A(2:,    :i-1)
                    X(:, i:   ) = A(2:, i+1:   )

                    d = s * det(X, n-1) * A(1, i) + d
                    s = -s
                end do
            end if
                
        end function

    end module Matrix

!   Tests

    program matrix_test
        use Matrix

        implicit none

        integer, parameter :: n = 3

!       Matrices
        double precision :: A(n, n), LU(n, n), L(n, n), U(n, n), I(n, n)

!       Vectors
        double precision :: b(n)

!       Define matrix
        A(1,1) = 1.0D0
        A(1,2) = 2.0D0
        A(1,3) = 4.0D0

        A(2,1) = 0.0D0
        A(2,2) = 1.0D0
        A(2,3) = 2.0D0

        A(3,1) = 1.0D0
        A(3,2) = 0.0D0
        A(3,3) = 1.0D0

!       Print Matrix Name
        write(*, *) 'A:'

!       Print Matrix
        call print_matrix(A, n, n)

!       Print its Determinant

        write(*, *) 'Matrix Det:', det(A, n)

        write(*, *) ''

!       Define Identity Matrix
        I = ID_MATRIX(n)

!       Print Matrix Name
        write(*, *) 'I:'

!       Print Matrix
        call print_matrix(I, n, n)

!       Print its Determinant

        write(*, *) 'Matrix Det:', det(I, n)

        write(*, *) ''

!       Define Vector
        b(1) = 3.0D0
        b(2) = 1.0D0
        b(3) = -5.0D0

!       Print Vector
        call print_vector(b, n)

!       LU - Decomposition
        call LU_DECOMP(A, L, U, n)

!       Print matrices
        write(*, *) 'L:'
        call print_matrix(L, n, n)

        write(*, *) ''

        write(*, *) 'U:'
        call print_matrix(U, n, n)

        write(*, *) ''

        write(*, *) 'LU:'
        LU = matmul(L, U)
        call print_matrix(LU, n, n)

    end program matrix_test
