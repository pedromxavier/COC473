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
            
        subroutine LU_DECOMP(A, b, n, L, U)
            implicit none

            integer :: n
            double precision :: A(n, n), b(n), L(n, n), U(n, n)
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

        recursive function det(A, n) result (D)
            implicit none

            integer :: n
            double precision :: A(n, n)
            double precision :: X(n-1, n-1)

            integer :: i, j, k
            double precision :: D

            if (n == 2) then
                D = A(1, 1) * A(2, 2) - A(1, 2) * A(2, 1)
            else
                D = 0.0D0
                do i = 1, n
                    X(:,  :i-1) = A(2:,    :i-1)
                    X(:, i:   ) = A(2:, i+1:   )

                    D = det(X, n-1) * A(1, i) + D
                end do
            end if
                
        end function

    end module Matrix

!   Tests

    program matrix_test
        use Matrix

        implicit none

        integer, parameter :: n = 3

        double precision :: A(n, n), L(n, n), U(n, n), I(n, n)

        double precision :: x(n)

        integer :: i, j

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
        call print_matrix(A, m, n)

!       Print its Determinant

        write(*, *) 'Matrix Det:', det(A, n)

        write(*, *) ''

!       Define Identity Matrix

!       Print Matrix Name
        write(*, *) 'I:'

!       Print Matrix
        call print_matrix(ID_MATRIX(n), n, n)

!       Print its Determinant

        write(*, *) 'Matrix Det:', det(ID_MATRIX(n), n)

        write(*, *) ''

!       Define Vector
        x(1) = 3.0D0
        x(2) = 1.0D0
        x(3) = -5.0D0

!       Print Vector
        call print_vector(x, m)

!       LU - Decomposition

        

    end program matrix_test
