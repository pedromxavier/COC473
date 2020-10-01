!   Calc Module

    module Calc
        implicit none
        integer :: MAX_ITER = 1000
        double precision :: h = 1.0D-10
        double precision :: TOL = 1.0D-8
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

        function bissection(f, a, b) result (x)
            double precision :: f
            double precision :: a, b, x

            if (b < a) then
                x = a
                a = b
                b = x
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

        function newton(f, df, x0) result (x)
            implicit none
            integer :: i

            double precision :: f, df
            double precision :: x0, x

            do i = 1, MAX_ITER
                x = x0 - f(x0) / df(x0)
                if (DABS(x - x0) > TOL) then
                    x0 = x
                else
                    return
                end if 
            end do
            call error("O método de Newton não convergiu.")
            return
        end function

        function secant(f, x0) result (x)
            implicit none
            integer :: i

            double precision :: f
            double precision :: xk(3), yk(2)
            double precision :: x0, x

            xk(1) = x0
            xk(2) = x0 + h
            do i = 1, MAX_ITER
                yk(2) = f(xk(2))
                xk(3) = xk(1) - yk(2) * (xk(2) - xk(1)) / (yk(2) - yk(1))
                
                if (DABS(xk(3) - xk(2)) > TOL) then
                    xk(1:2) = xk(2:3)
                    yk(1) = yk(2)
                else
                    x = xk(3)
                    return
                end if 
            end do

            call error("O método da Secante não convergiu.")
            return
            
            return
        end function


      
    end module Calc
