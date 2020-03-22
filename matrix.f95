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

            do i = 1, m
                read(33,*) A(i,:)
            end do

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

            return

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

            call warn('Aviso: Matriz mal-condicionada.')
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

            logical :: x

            x = POSITIVE_DEFINITE(A, n)

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

        function PLU_DECOMP(A, P, L, U, n) result (x)
            implicit none

            integer :: n
            double precision :: A(n,n), P(n,n), L(n,n), U(n,n)

            logical :: x

!           Permutation Matrix
            P = PIVOT_MATRIX(A, n)

!           Decomposition over Row-Swapped Matrix
            x = LU_DECOMP(matmul(P, A), L, U, n)
            return
        end function
            
        function LU_DECOMP(A, L, U, n) result (x)
            implicit none

            integer :: n
            double precision :: A(n, n), L(n, n), U(n,n), M(n, n)

            logical :: x

            integer :: i, j, k

!           Results Matrix
            M(:, :) = A(:, :)

            if (.NOT. LU_COND(A, n)) then
                call ILL_COND()
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

            call LU_MATRIX(M, L, U, n)

            x = .TRUE.
            return

        end function

        function CHOLESKY_DECOMP(A, L, n) result (x)
            implicit none

            integer :: n
            double precision :: A(n, n), L(n, n)

            logical :: x

            integer :: i, j

            if (.NOT. CHOLESKY_COND(A, n)) then
                call ILL_COND()
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

        subroutine warn(text)
            implicit none
            character(len=*) :: text
            write (*, *) ''//achar(27)//'[31m'//text//''//achar(27)//'[0m'
        end subroutine

        subroutine info(text)
            implicit none
            character(len=*) :: text
            write (*, *) ''//achar(27)//'[32m'//text//''//achar(27)//'[0m'
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

            if (n == 1) then
                d = A(1, 1)
            elseif (n == 2) then
                d = A(1, 1) * A(2, 2) - A(1, 2) * A(2, 1)
            else
                d = 0.0D0
                s = 1.0D0
                do i = 1, n
!                   Compute submatrix X
                    X(:,  :i-1) = A(2:,    :i-1)
                    X(:, i:   ) = A(2:, i+1:   )

                    d = s * DET(X, n-1) * A(1, i) + d
                    s = -s
                end do
            end if
                
        end function

    end module Matrix

!   Tests

    program matrix_test
        use Matrix

        implicit none

!       Command-line Args
        integer :: argc

        character(len=32) :: matrix_fname, vector_fname

!       Matrices
        integer :: m, n

        double precision, allocatable :: A(:, :)

!       Determinant
        double precision :: d

!       For LU decomposition & Cholesky also
        double precision, allocatable :: L(:, :), U(:, :)

!       For PLU decomposition
        double precision, allocatable :: P(:, :)

!       Vectors
!       double precision, allocatable :: b(:)

!       Get Command-Line Args
        if (argc == 1) then
            call getarg(1, matrix_fname)
        elseif(argc == 2) then
            call getarg(1, matrix_fname)
            call getarg(2, vector_fname)
        elseif (argc > 2) then
            goto 91
        else
            matrix_fname = 'matrix.txt'
        end if

        call read_matrix(matrix_fname, A, m, n)

        if (m /= n) then
            goto 90
        end if

!       Print Matrix Name
        write(*, *) 'A:'

!       Print Matrix
        call print_matrix(A, m, n)

!       Print its Determinant
        d = DET(A, n)
        write(*, *) 'Matrix Det:', d

        if (d == 0.0D0) then
            goto 92
        endif

!       Decomposition Time!
!       Allocate Result Matrices
        allocate(P(n, n))
        allocate(L(n, n))
        allocate(U(n, n))

!       Try LU Decomposition
       call info(':: Decomposição LU (sem pivoteamento) ::')
        if (.NOT. LU_DECOMP(A, L, U, n)) then
            goto 81
        end if

        write(*, *) 'A:'
        call print_matrix(A, n, n)
        
        write(*, *) 'L:'
        call print_matrix(L, n, n)

        write(*, *) 'U:'
        call print_matrix(U, n, n)

!       Try PLU Decomposition
81      call info(':: Decomposição PLU (com pivoteamento) ::')
        if (.NOT. PLU_DECOMP(A, P, L, U, n)) then
            goto 82
        end if

        write(*, *) 'A:'
        call print_matrix(A, n, n)

        write(*, *) 'P:'
        call print_matrix(P, n, n)
        
        write(*, *) 'L:'
        call print_matrix(L, n, n)

        write(*, *) 'U:'
        call print_matrix(U, n, n)

!       Try Cholesky Decomposition
82      call info(':: Decomposição de Cholesky ::')
        if (.NOT. CHOLESKY_DECOMP(A, L, n)) then
            goto 85
        end if

        write(*, *) 'A:'
        call print_matrix(A, n, n)

        write(*, *) 'L:'
        call print_matrix(L, n, n)

85      call info('Sucesso!')
        deallocate(P)
        deallocate(L)
        deallocate(U)
        goto 100

90      call warn('Essa matriz não é quadrada, logo o programa não faz sentido.')
        goto 100
91      call warn('Parâmetros em excesso: O programa espera apenas o nome do arquivo da matriz.')
        goto 100
92      call warn('O Determinante é 0! Os métodos não se aplicam nesse caso.')
        goto 100

100     stop

    end program matrix_test
