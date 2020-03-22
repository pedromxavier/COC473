!   Matrix Module

    module Matrix
        implicit none

    contains

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

            do i = 1, n
                read(33,*) A(i,:)
                write(*,*) A(i,:)
            enddo

            close(33)
        end subroutine

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

        recursive function POSITIVE_DEFINITE(A, n) result (x)
!       Checks wether a matrix is positive definite
!       according to Sylvester's criterion.
            implicit none

            integer :: n
            double precision A(n, n)

            logical :: x

            if (n == 1) then
                x = (A(1, 1) > 0)
            else
                x = POSITIVE_DEFINITE(A(:n-1, :n-1), n-1) .AND. (DET(A, n) > 0)
            end if

            return

        end function

        function SYMMETRICAL(A, n) result (x)
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

        subroutine ILL_COND()
            implicit none

            print *, ''//achar(27)//'[95m Matriz mal-condicionada '//achar(27)//'[0m.'
        end subroutine

        function CHOLESKY_COND(A, n) result (x)
            implicit none

            integer :: n
            double precision :: A(n, n)

            logical :: x

            x = SYMMETRICAL(A, n) .AND. POSITIVE_DEFINITE(A, n)
            return

        end function

        function PLU_COND(A, n) result (x)
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

        function LU_COND(A, n) result (x)
            implicit none

            integer :: n
            double precision A(n, n)

            integer :: i, j
            logical :: x

            do j = 1, n
                do i = 1, j
                    if (A(i, j) < 0.01D0) then
                        x = .FALSE.
                        return
                    end if 
                end do
            end do
            x = .TRUE.
            return
        end function

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

            return 
        end function

        function PIVOT_MATRIX(A, n) result (P)
            implicit none

            integer :: n
            double precision :: A(n, n)
            double precision :: P(n, n)

            integer :: j, k

            P = ID_MATRIX(n)

            do j = 1, n
                k = ROW_MAX(A, j, n)
                if (j /= k) then 
                    call SWAP_ROWS(P, j, k, n)
                end if
            end do

            return
        end function

        function PLU_DECOMP(A, P, n) result (X)
            implicit none

            integer :: n
            double precision :: P(n, n), A(n, n), PA(n, n), X(n, n)

            if (.NOT. PLU_COND(A, n)) then
                call ILL_COND()
                return
            end if

!           Permutation Matrix
            P = PIVOT_MATRIX(A, n)

!           Row-Swapped Matrix
            PA = matmul(P, A)

            X = LU_DECOMP(PA, n)
            return

        end function
            
        function LU_DECOMP(A, n) result (X)
            implicit none

            integer :: n
            double precision :: A(n, n), X(n, n)

            integer :: i, j, k

!           Results Matrix
            X(:, :) = A(:, :)

            if (.NOT. LU_COND(A, n)) then
                call ILL_COND()
                return
            end if

            do k = 1, n-1
                do i = k+1, n
                    X(i, k) = X(i, k) / X(k, k)
                end do
                
                do j = k+1, n
                    do i = k+1, n
                        X(i, j) = X(i, j) - X(i, k) * A(k, j)
                    end do
                end do
            end do

            return

        end function

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

20          format('|', F10.5, ' ')
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

        subroutine LU_MATRIX(A, L, U, n)
            implicit none

            integer :: n
            double precision :: A(n, n), L(n, n), U(n, n)

            integer :: i

            do i = 1, n
                L(i, i) = 1.0D0
                if (i >= 1) then
                    L(i,  :i-1) = A(i,  :i-1)
                end if
                U(i, i:) = A(i, i:)
            end do
        end subroutine

        subroutine print_vector(x, n)
            implicit none

            integer :: n
            double precision :: x(n)
            
            integer :: i

30          format('|', F10.5, '|')

            do i = 1, n
                write(*, 30) x(i)
            end do
        end subroutine

        recursive function DET(A, n) result (d)
            implicit none

            integer :: n
            double precision :: A(n, n)
            double precision :: X(n-1, n-1)

            integer :: i
            double precision :: d, s

            if (n == 2) then
                d = A(1, 1) * A(2, 2) - A(1, 2) * A(2, 1)
            else
                d = 0.0D0
                s = 1.0D0
                do i = 1, n
                    X(:,  :i-1) = A(2:,    :i-1)
                    X(:, i:   ) = A(2:, i+1:   )

                    d = s * DET(X, n-1) * A(1, i) + d
                    s = -s
                end do
            end if
                
        end function

!       ================== TEST FUNCTIONS BELOW ================



    end module Matrix

!   Tests

    program matrix_test
        use Matrix

        implicit none

        integer :: m, n

!       Matrices
        double precision, allocatable :: A(:, :)

!       Vectors
!       double precision :: b(n)

!       Define matrix
        call read_matrix('matrix.txt', A, m, n)

!       Print Matrix Name
        write(*, *) 'A:'

!       Print Matrix
        call print_matrix(A, n, n)

!       Print its Determinant

        write(*, *) 'Matrix Det:', DET(A, n)
    end program matrix_test
