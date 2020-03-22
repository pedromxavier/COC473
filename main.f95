program main
    use Matrix

    implicit none

!   Command-line Args
    integer :: argc

    character(len=32) :: matrix_fname, vector_fname

!   Matrices
    integer :: m, n
    double precision, allocatable :: A(:, :)

!   Determinant
    double precision :: d

!   For LU decomposition & Cholesky also
    double precision, allocatable :: L(:, :), U(:, :)

!   For PLU decomposition
    double precision, allocatable :: P(:, :)

!   Vectors
    integer :: t
    double precision, allocatable :: b(:)

!   Get Command-Line Args
    matrix_fname = 'matrix.txt'
    vector_fname = 'vector.txt'

    if (argc == 1) then
        call getarg(1, matrix_fname)
    elseif(argc == 2) then
        call getarg(1, matrix_fname)
        call getarg(2, vector_fname)
    elseif (argc > 2) then
        goto 92
    end if

    call read_matrix(matrix_fname, A, m, n)
    call read_vector(vector_fname, b, t)

    if (m /= n) then
        goto 90
    elseif (m /= t) then
        goto 91
    end if



!   Print Matrix
!   Print Matrix Name
    write(*, *) 'A:'
    call print_matrix(A, m, n)

!   Print Matrix
!   Print Matrix Name
    write(*, *) 'b:'
    call print_vector(b, t)

!   Print its Determinant
    d = DET(A, n)
    write(*, *) 'Matrix Det:', d

    if (d == 0.0D0) then
        goto 93
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

85  call info('Sucesso!')
    deallocate(A)
    deallocate(P)
    deallocate(L)
    deallocate(U)
    deallocate(b)
    goto 100

90  call error('Essa matriz não é quadrada! Assim o programa não faz sentido.')
    goto 100
91  call error('A matriz `A` e o vetor `b` não possuem a mesma dimensão. Não tem como fazer Ax = b!')
    goto 100
92  call error('Parâmetros em excesso: O programa espera apenas o nome do arquivo da matriz.')
    goto 100
93  call error('O Determinante é 0! Os métodos não se aplicam nesse caso.')
    goto 100

100     stop

end program main