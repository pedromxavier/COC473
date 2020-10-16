!   Util Module
    module Util
        implicit none
        character, parameter :: ENDL = ACHAR(10)
        character, parameter :: TAB = ACHAR(9)

        double precision :: DINF, DNINF, DNAN
        !DATA DINF/x'7ff0000000000000'/, DNINF/x'fff0000000000000'/, DNAN/x'7ff8000000000000'/

        double precision :: PI = 4.0D0 * DATAN(1.0D0)

        logical :: ENABLE_DEBUG = .FALSE.

        type StringArray
            character (:), allocatable :: str
        end type StringArray
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
        function STR(k) result (t)
!       "Convert an integer to string."
            integer, intent(in) :: k
            character(len=128) :: s
            character(len=:), allocatable :: t
            write(s, *) k
            t = TRIM(ADJUSTL(s))
            return
            return
        end function

        function DSTR(x) result (t)
            double precision, intent(in) :: x
            character(len=128) :: s
            character(len=:), allocatable :: t
            write(s, *) x
            t = TRIM(ADJUSTL(s))
            return
        end function

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

        subroutine debug(text)
!           Yellow Text
            implicit none
            character(len=*) :: text
            if (ENABLE_DEBUG) then
                write (*, *) 'DEBUG: '//achar(27)//'[93m'//text//''//achar(27)//'[0m'
            end if
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

        subroutine linspace(a, b, dt, n, t)
            implicit none
            integer :: k, n
            double precision, intent(in) :: a, b, dt
            double precision, dimension(:), allocatable :: t
            n = 1 + (b - a) / dt
            allocate(t(n))
            t(:) = dt * (/ (k, k=0, n-1) /)
        end subroutine
    end module Util
