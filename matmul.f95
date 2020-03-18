!   Matmul Fortran
    program matmul
        implicit none

        integer, parameter :: N = 20000, v = 1, s = 500
        integer :: k, i
        
        double precision :: mean, stddev
        double precision :: A(N,N), x(N), y(N)

        double precision :: t_c, t_f
        double precision :: t__c(v), t__f(v)

        double precision :: t, dt

        open (unit=1, file="f_result_c.txt")
        open (unit=2, file="f_result_f.txt")

10      format (i5,"|",f22.12,"$",f22.12,"|200")

        do k=1000, N, s
            do i=1, v
                call MatMulC(A, x, y, k, t_c)
                call MatMulF(A, x, y, k, t_f)
                T__c(i) = t_c
                T__f(i) = t_f
            end do

            t = mean(T__c, v)
            dt = stddev(T__c, v)
            write (1, 10) k, t, dt

            t = mean(T__f, v)
            dt = stddev(T__f, v)
            write (2, 10) k, t, dt
        end do
    end program matmul

    function mean(x, n)
        integer n
        double precision mean, x(n), s
        s = 0
        do i = 1, n
            s = s + x(i)
        end do
        mean = s/n
    end function

    function stddev(x, n)
        integer n
        double precision mean, stddev, var, avg, x(n)

        if (n < 2) then
            stddev = 0
        else
            var = 0
            avg = mean(x, n)
            do i = 1, n
                var = var + (x(i) - avg)**2
            end do
            var = var/(n-1)
            stddev = sqrt(var)
        end if
    end function

    subroutine zero_vector(y, n)
        integer n
        double precision y(n)
        do i=1, n
            y(i) = 0
        end do
    end subroutine

    subroutine random_vector(x, n)
        integer n
        double precision x(n)
        do i=1, n
            x(i) = rand()
        end do
    end subroutine

    subroutine random_matrix(A, m, n)
        integer m, n
        double precision A(m,n)
        do i = 1, m
            do j = 1, n
                A(i,j) = rand()
            end do
        end do
    end subroutine

    subroutine show_vector(y, n)
        integer n
        double precision y(n)
        do i=1, n
            print *, y(i)
        end do
    end subroutine

    subroutine MatMulC(A, x, y, n, t)
        integer n
        double precision A(n,n), x(n), y(n), before, after
        double precision t

        call random_matrix(A, n, n)
        call random_vector(x,n)
        call zero_vector(y,n)

        call cpu_time(before)
        do i = 1, n
            do j = 1, n
                y(i) = y(i) + A(i, j)*x(j)
            end do
        end do
        call cpu_time(after)
        t = after - before
    end subroutine

    subroutine MatMulF(A, x, y, n, t)
        integer n
        double precision A(n,n), x(n), y(n), before, after
        double precision t

        call random_matrix(A, n, n)
        call random_vector(x,n)
        call zero_vector(y,n)

        call cpu_time(before)
        do j = 1, n
            do i = 1, n
                y(i) = y(i) + A(i, j)*x(j)
            end do
        end do
        call cpu_time(after)
        t = after - before
    end subroutine