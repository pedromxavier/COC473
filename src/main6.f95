program main6
    use Func
    use Calc
    use Util
    use Plotlib
    implicit none

!   Command-line Args
    integer :: argc

    DEBUG_MODE = .TRUE.

!   Random seed definition
    call init_random_seed()

!   Get Command-Line Args
    argc = iargc()

    if (argc == 0) then
        goto 100
    else
        goto 11
    end if

!   ====== Begin =====================================    

!   ====== Success ===================================
10  call info(':: Sucesso ::')
    goto 1
!   ====== Errors ====================================
11  call error('Este programa não aceita parâmetros.')
    goto 1
!   ====== Finish ====================================
1   stop
!   ====== Program ===================================
100 call Q1
    goto 200

200 call Q2
    goto 300

300 call Q3
    goto 400

400 call QB
    goto 10

    contains
        subroutine Q1
            implicit none
            integer :: n, i, j
            integer :: k = 0
            double precision, dimension(:), allocatable :: t, y
            double precision, dimension(4) :: dt
            type(StringArray), dimension(:), allocatable :: legend
            type(StringArray), dimension(:), allocatable :: title

            allocate(legend(4), title(4))

            dt(1) = 5.0D-1
            dt(2) = 2.5D-1
            dt(3) = 1.0D-1
            dt(4) = 0.1D-1

            title(1)%str = "Δt = 0.50"
            title(2)%str = "Δt = 0.25"
            title(3)%str = "Δt = 0.10"
            title(4)%str = "Δt = 0.01"

            legend(1)%str = "Euler"
            legend(2)%str = "Runge-Kutta II"
            legend(3)%str = "Runge-Kutta IV"
            legend(4)%str = "y(t)"

            call blue(ENDL//"1)"//ENDL//F13_NAME)

            call info(":: Plotagem ::")

            call begin_plot(fname='L6-Q1', size_w='12in', size_h='9in')
            call subplots(2, 2)

            do i=1,2
                do j=1,2
                    k = k + 1

                    call linspace(F13_A, F13_B, dt(k), n, t)

                    allocate(y(n))

                    !call info(":: Método de Euler ::")
                    y = ode_solve(df13, F13_Y0, t, n, kind='euler')
                    call subplot(i, j, t, y, n)

                    !call info(":: Método de Runge-Kutta de 2ª ordem ::")
                    y = ode_solve(df13, F13_Y0, t, n, kind='runge-kutta2')
                    call subplot(i, j, t, y, n)

                    !call info(":: Método de Runge-Kutta de 4ª ordem ::")
                    y = ode_solve(df13, F13_Y0, t, n, kind='runge-kutta4')
                    call subplot(i, j, t, y, n)

                    y = (/ (f13(t(j)), j=1, n) /)
                    call subplot(i, j, t, y, n)

                    call subplot_config(i, j, &
                    title=title(k)%str, xlabel='t', ylabel='y(t)', grid=.TRUE., legend=legend)

                    deallocate(t, y)
                end do
            end do

            call render_plot(clean=.TRUE.)

            deallocate(legend, title)
        end subroutine

        subroutine Q2
            integer :: n, i, j
            integer :: k = 0
            double precision, dimension(:), allocatable :: t, y
            double precision, dimension(4) :: dt
            type(StringArray), dimension(:), allocatable :: legend
            type(StringArray), dimension(:), allocatable :: title

            allocate(legend(2), title(4))

            dt(1) = 5.0D-1
            dt(2) = 2.5D-1
            dt(3) = 1.0D-1
            dt(4) = 0.1D-1

            title(1)%str = "Δt = 0.50"
            title(2)%str = "Δt = 0.25"
            title(3)%str = "Δt = 0.10"
            title(4)%str = "Δt = 0.01"

            legend(1)%str = "Série de Taylor"
            legend(2)%str = "Runge-Kutta-Nystrom"

            call blue(ENDL//"2)"//ENDL//F14_NAME)

            call info(":: Plotagem ::")

            call begin_plot(fname='L6-Q2', size_w='12in', size_h='9in')
            call subplots(2, 2)

            do i=1,2
                do j=1,2
                    k = k + 1

                    call linspace(F14_A, F14_B, dt(k), n, t)

                    allocate(y(n))

                    y = ode2_solve(d2f14, F14_Y0, F14_DY0, t, n, kind='taylor')
                    call subplot(i, j, t, y, n)

                    y = ode2_solve(d2f14, F14_Y0, F14_DY0, t, n, kind='runge-kutta-nystrom')
                    call subplot(i, j, t, y, n)

                    call subplot_config(i, j, &
                    title=title(k)%str, xlabel='t', ylabel='y(t)', grid=.TRUE., legend=legend)

                    deallocate(t, y)
                end do
            end do

            call render_plot(clean=.TRUE.)

            deallocate(legend, title)
        end subroutine

        subroutine Q3
            integer :: n, i, j
            integer :: k = 0
            double precision, dimension(:), allocatable :: t, y
            double precision, dimension(4) :: dt
            type(StringArray), dimension(:), allocatable :: legend
            type(StringArray), dimension(:), allocatable :: title

            allocate(legend(2), title(4))

            dt(1) = 5.0D-1
            dt(2) = 4.0D-1
            dt(3) = 3.0D-1
            dt(4) = 2.0D-1

            title(1)%str = "Δt = 0.50"
            title(2)%str = "Δt = 0.40"
            title(3)%str = "Δt = 0.30"
            title(4)%str = "Δt = 0.20"

            legend(1)%str = "Série de Taylor"
            legend(2)%str = "Runge-Kutta-Nystrom"

            call blue(ENDL//"3)"//ENDL//F15_NAME)

            call begin_plot(fname='L6-Q3', size_w='12in', size_h='9in')
            call subplots(2, 2)

            do i=1,2
                do j=1,2
                    k = k + 1

                    call linspace(F15_A, F15_B, dt(k), n, t)

                    allocate(y(n))

                    y = ode2_solve(d2f15, F15_Y0, F15_DY0, t, n, kind='taylor')
                    call subplot(i, j, t, y, n)

                    y = ode2_solve(d2f15, F15_Y0, F15_DY0, t, n, kind='runge-kutta-nystrom')
                    call subplot(i, j, t, y, n)

                    call subplot_config(i, j, &
                    title=title(k)%str, xlabel='t', ylabel='z(t)', grid=.TRUE., legend=legend)

                    deallocate(t, y)
                end do
            end do

            call render_plot(clean=.TRUE.)

            deallocate(legend, title)
        end subroutine

        subroutine QB
            integer :: n, i, j
            integer :: k = 0
            double precision, dimension(:), allocatable :: t, y
            double precision, dimension(4) :: dt
            type(StringArray), dimension(:), allocatable :: legend
            type(StringArray), dimension(:), allocatable :: title

            allocate(legend(2), title(4))

            dt(1) = 5.0D-1
            dt(2) = 2.5D-1
            dt(3) = 1.0D-1
            dt(4) = 0.1D-1

            title(1)%str = "Δt = 0.50"
            title(2)%str = "Δt = 0.25"
            title(3)%str = "Δt = 0.10"
            title(4)%str = "Δt = 0.01"

            legend(1)%str = "Série de Taylor"
            legend(2)%str = "Runge-Kutta-Nystrom"

            call blue(ENDL//"Bônus)"//ENDL//F15_NAME)

            call info(":: Plotagem ::")

            call begin_plot(fname='L6-QB', size_w='12in', size_h='9in')
            call subplots(2, 2)

            do i=1,2
                do j=1,2
                    k = k + 1

                    call linspace(F15_A, F15_B, dt(k), n, t)

                    allocate(y(n))

                    y = ode2_solve(d2f15, F15_BY0, F15_DY0, t, n, kind='taylor')
                    call subplot(i, j, t, y, n)

                    y = ode2_solve(d2f15, F15_BY0, F15_DY0, t, n, kind='runge-kutta-nystrom')
                    call subplot(i, j, t, y, n)

                    call subplot_config(i, j, &
                    title=title(k)%str, xlabel='t', ylabel='z(t)', grid=.TRUE., legend=legend)

                    deallocate(t, y)
                end do
            end do

            call render_plot(clean=.TRUE.)

            deallocate(legend, title)
        end subroutine
end program main6