program main2
    use Matrix

    implicit none

!   Command-line Args
    integer :: argc

    character(len=32) :: matrix_fname, vector_fname

!   Matrices
    integer :: m, n
    double precision, allocatable :: A(:, :)
    double precision, allocatable :: LL(:, :)
    double precision, allocatable :: XX(:, :)

!   Vectors
    integer :: t
    double precision, allocatable :: b(:)
    double precision, allocatable :: x(:)

!   Eigenvalue
    double precision :: l

!   Get Command-Line Args
    matrix_fname = 'matrix2.txt'
    vector_fname = 'vector2.txt'

    argc = iargc()

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

!   Eigenvalue Time
!   Power Method
    allocate(x(n))
    call info(':: Método das Potências (Power Method) ::')
    if (.not. power_method(A, n, x, l)) then
        goto 80
    end if 

    write(*, *) 'x:'
    call print_vector(x, n)

    write(*, *) 'lambda:'
    write(*, *) l

!   Jacobi Eigenvalues/vectors time
    allocate(LL(n, n))
    allocate(XX(n, n))
    call info(':: Método de autovalores de Jacobi ::')
    if (.not. Jacobi_eigen(A, n, LL, XX)) then
        goto 81
    end if
    write(*, *) 'L:'
    call print_matrix(LL, n, n)
    write(*, *) 'X:'
    call print_matrix(XX, n, n)
    goto 85

80  call error('O Método das potências não convergiu a tempo.')
    write(*, *) 'x:'
    call print_vector(x, n)
    goto 100

81  call error('O Método de autovalores de Jacobi não convergiu a tempo.')
    write(*, *) 'L:'
    call print_matrix(LL, n, n)
    write(*, *) 'X:'
    call print_matrix(XX, n, n)
    goto 100

85  call info('Sucesso!')
    deallocate(A)
    deallocate(b)
    goto 100

90  call error('Essa matriz não é quadrada! Assim o programa não faz sentido.')
    goto 100
91  call error('A matriz `A` e o vetor `b` não possuem a mesma dimensão. Não tem como fazer Ax = b!')
    goto 100
92  call error('Parâmetros em excesso: O programa espera apenas o nome do arquivo da matriz.')
    goto 100

100     stop

end program main2