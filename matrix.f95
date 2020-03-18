!   Matrix Module

    module Matrix
        implicit none

    contains

        subroutine LU_DECOMP()
            implicit none

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

    end module Matrix

!   Tests

    program matrix_test
        use Matrix

        implicit none

        integer, parameter :: m = 3
        integer, parameter :: n = 3

        double precision :: A(m, n)

        double precision :: x(m)

        integer :: i, j

!       Define matrix
        do i = 1, m
            do j = 1, n
                if (i == j) then
                    A(i, j) = -1.0D0
                else
                    A(i, j) = 0.0D0
                end if
            end do
        end do
!       Print Matrix
        call print_matrix(A, m, n)

        write(*, *) ''

!       Define Vector
        x(1) = 3.0D0
        x(2) = 1.0D0
        x(3) = -5.0D0

!       Print Vector
        call print_vector(x, m)

    end program matrix_test