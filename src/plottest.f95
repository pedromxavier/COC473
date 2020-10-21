program plottest
    use Util
    use PlotLib
    implicit none

    integer :: m, n, k, i, j
    double precision, dimension(:), allocatable :: x, y, z

    DEBUG_MODE = .TRUE.

    m = 3
    n = 2
    k = 100

    allocate(x(k))
    allocate(y(k))
    allocate(z(k))

    x(:) = (/(i, i=1, k)/) / 100.0D0
    y(:) = (/(2 * x(j), j=1, k)/)
    z(:) = (/(x(j) * x(j), j=1, k)/)

    call begin_plot(fname='plottest')
    call subplots(m, n)

    do i=1, m
        do j=1,n
            call subplot(i, j, x, y, k)
            call subplot(i, j, x, z, k)
            call subplot_config(i, j, title='subplot ('//STR(i)//','//STR(j)//')', &
            xlabel='t', ylabel='z', grid=.TRUE.)
        end do
    end do

    call render_plot()

    deallocate(x, y, z)

    stop
end program plottest