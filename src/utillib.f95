!   Util Module
    module Util
        implicit none
        character, parameter :: ENDL = ACHAR(10)
        character, parameter :: TAB = ACHAR(9)
    contains
!       ==== Random seed Initialization ====
        subroutine init_random_seed()
            integer :: i, n, clock
            integer, allocatable :: seed(:)
        
            call RANDOM_SEED(SIZE=n)
            allocate(seed(n))
            call SYSTEM_CLOCK(COUNT=clock)
            seed = clock + 37 * (/ (i - 1, i = 1, n) /)
            call RANDOM_SEED(PUT=seed)
            deallocate(seed)
        end subroutine

        function DRAND(a, b) result (y)
            implicit none
            double precision :: a, b, x, y
            ! x in [0, 1)            
            call RANDOM_NUMBER(x)
            y = (x * (b - a)) + a
            return
        end function

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

        subroutine show(var, value)
!           Violet Text
            implicit none
            character(len=*) :: var
            character(len=24) :: val
            double precision :: value
10          format(F24.16, '')
            write (val, 10) value
            write (*, *) ''//achar(27)//'[36m'//var//' = '//val//''//achar(27)//'[0m'
        end subroutine

        subroutine show_vector(var, x, n)
            implicit none
            integer :: n
            character(len=*) :: var
            double precision :: x(n)
            integer :: i
30          format(' |', F30.12, '|')
            write (*, *) ''//achar(27)//'[36m'//var//' = '
            do i = 1, n
                write(*, 30) x(i)
            end do
            write (*, *) ''//achar(27)//'[0m'
        end subroutine

        recursive subroutine cross_quick_sort(x, y, u, v, n)
            integer :: n, i, j, u, v
            double precision :: p, aux, auy
            double precision :: x(n), y(n)

            i = u
            j = v

            p = x((u + v) / 2)

            do while (i <= j)
                do while (x(i) < p)
                    i = i + 1
                end do
                do while(x(j) > p)
                    j = j - 1
                end do
                if (i <= j) then
                    aux = x(i)
                    auy = y(i)
                    x(i) = x(j)
                    y(i) = y(j)
                    x(j) = aux
                    y(j) = auy
                    i = i + 1
                    j = j - 1
                end if
            end do

            if (u < j) then
                call cross_quick_sort(x, y, u, j, n)
            end if
            if (i < v) then
                call cross_quick_sort(x, y, i, v, n)
            end if
            return
        end subroutine

        subroutine cross_sort(x, y, n)
            implicit none
            integer :: n
            double precision :: x(n), y(n)
            
            call cross_quick_sort(x, y, 1, n, n)
        end subroutine
    end module Util
