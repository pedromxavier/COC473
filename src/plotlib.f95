module Plot
    use Util
    implicit none
    character(len=*), parameter :: DEFAULT_FNAME = 'plotfile'
    character(len=*), parameter :: PLOT_ENDL = ',\'//ENDL
    contains

    function PLOT_FNAME(fname) result (path)
        implicit none
        character(len=*) :: fname
        character(len=:), allocatable :: path
        path = 'plot/'//fname//'.plt'
        return
    end function

    function DATA_FNAME(fname) result (path)
        implicit none
        character(len=*) :: fname
        character(len=:), allocatable :: path
        path = 'plot/'//fname//'.dat'
        return
    end function

    function OUTP_FNAME(fname) result (path)
        implicit none
        character(len=*) :: fname
        character(len=:), allocatable :: path
        path = 'plot/'//fname//'.pdf'
        return
    end function

    function GNU_PLOT_CMD(fname) result (path)
        implicit none
        character(len=*) :: fname
        character(len=:), allocatable :: path
        path = 'gnuplot -p '//PLOT_FNAME(fname)
    end function

!   =========================== Plot =========================
    subroutine plot_config(fname, title, xlabel, ylabel, grid, points)
        implicit none
        integer :: file

        logical, optional :: grid
        logical :: t_grid

        character(len=:), allocatable :: set_grid

        logical, optional :: points
        logical :: t_points

        character(len=:), allocatable :: set_points

        character(len=*), optional :: fname
        character(len=:), allocatable :: t_fname

        character(len=*), optional :: title
        character(len=:), allocatable :: t_title

        character(len=*), optional :: ylabel
        character(len=:), allocatable :: t_ylabel

        character(len=*), optional :: xlabel
        character(len=:), allocatable :: t_xlabel

        character(len=:), allocatable :: config

        if (.NOT. PRESENT(fname)) then
            t_fname = DEFAULT_FNAME
        else
            t_fname = fname
        end if

        if (.NOT. PRESENT(grid)) then
            t_grid = .TRUE.
        else
            t_grid = grid
        end if

        if (.NOT. PRESENT(points)) then
            t_points = .FALSE.
        else
            t_points = points
        end if

        if (.NOT. PRESENT(title)) then
            t_title = ''
        else
            t_title = title
        end if

        if (.NOT. PRESENT(xlabel)) then
            t_xlabel = 'x'
        else
            t_xlabel = xlabel
        end if

        if (.NOT. PRESENT(ylabel)) then
            t_ylabel = 'y'
        else
            t_ylabel = ylabel
        end if

        if (t_grid) then
            set_grid = "set grid"
        else
            set_grid = ""
        end if

        if (t_points) then
            set_points = "with linespoints"
        else
            set_points = ""
        end if

        config = "#"//OUTP_FNAME(t_fname)//ENDL// &
        'set size 1,1'//ENDL// &
        'set terminal pdf'//ENDL// &
        'set output "'//OUTP_FNAME(t_fname)//'"'//ENDL// &
        'set title "'//t_title//'"'//ENDL// &
        'set nokey'//ENDL// &
        set_grid//ENDL// &
        'set xlabel "'//t_xlabel//'"'//ENDL// &
        'set ylabel "'//t_ylabel//'"'//ENDL// &
        'xydata = "'//DATA_FNAME(t_fname)//'"'//ENDL// &
        'plot xydata using 1:2 '//set_points

        open (newunit=file, action='write', file=PLOT_FNAME(t_fname), status='replace')
            write(file, *) config
        close(file)
    end subroutine 

    subroutine xy_plot(x, y, n, fname, title, xlabel, ylabel, grid, points)
        implicit none
        integer :: file, i
        integer, intent(in) :: n
        double precision, dimension(n), intent(in) :: x, y

        logical, optional :: grid
        logical :: t_grid

        logical, optional :: points
        logical :: t_points

        character(len=*), optional :: title
        character(len=:), allocatable :: t_title

        character(len=*), optional :: ylabel
        character(len=:), allocatable :: t_ylabel

        character(len=*), optional :: xlabel
        character(len=:), allocatable :: t_xlabel

        character(len=*), optional :: fname
        character(len=:), allocatable :: t_fname

        if (.NOT. PRESENT(fname)) then
            t_fname = DEFAULT_FNAME
        else
            t_fname = fname
        end if

        if (.NOT. PRESENT(grid)) then
            t_grid = .TRUE.
        else
            t_grid = grid
        end if

        if (.NOT. PRESENT(points)) then
            t_points = .FALSE.
        else
            t_points = points
        end if

        if (.NOT. PRESENT(title)) then
            t_title = ''
        else
            t_title = title
        end if

        if (.NOT. PRESENT(xlabel)) then
            t_xlabel = 'x'
        else
            t_xlabel = xlabel
        end if

        if (.NOT. PRESENT(ylabel)) then
            t_ylabel = 'y'
        else
            t_ylabel = ylabel
        end if

!       ===================== Write to data file ========================
        open (newunit=file, action='write', file=DATA_FNAME(t_fname), status='replace')
        do i=1,n
            write(file, *) x(i), y(i)
        end do
        close(file)
!       =================================================================
        call plot_config(fname=t_fname, title=t_title, xlabel=t_xlabel, ylabel=t_ylabel, grid=t_grid)
        call EXECUTE_COMMAND_LINE(GNU_PLOT_CMD(t_fname))
    end subroutine

!   ======================= MultiPlot ===================================
    subroutine multiplot_config(m, fname, title, xlabel, ylabel, legend, grid, with)
        implicit none
        integer :: file, i, m

        logical, optional :: grid
        logical :: t_grid

        character(len=:), allocatable :: set_grid

        character(len=*), optional :: with
        character(len=:), allocatable :: t_with

        character(len=*), optional :: fname
        character(len=:), allocatable :: t_fname

        character(len=*), optional :: title
        character(len=:), allocatable :: t_title

        character(len=*), optional :: ylabel
        character(len=:), allocatable :: t_ylabel

        character(len=*), optional :: xlabel
        character(len=:), allocatable :: t_xlabel

        character(len=:), allocatable :: config

        type(StringArray), dimension(:), optional :: legend
        type(StringArray), dimension(m) :: t_legend

        if (.NOT. PRESENT(legend)) then
            do i=1,m
                t_legend(i)%str = STR(i)
            end do
        else
            t_legend(:) = legend(:)
        end if

        if (.NOT. PRESENT(fname)) then
            t_fname = DEFAULT_FNAME
        else
            t_fname = fname
        end if

        if (.NOT. PRESENT(grid)) then
            t_grid = .TRUE.
        else
            t_grid = grid
        end if

        if (.NOT. PRESENT(with)) then
            t_with = ""
        else
            t_with = "with "//with
        end if

        if (.NOT. PRESENT(title)) then
            t_title = ''
        else
            t_title = title
        end if

        if (.NOT. PRESENT(xlabel)) then
            t_xlabel = 'x'
        else
            t_xlabel = xlabel
        end if

        if (.NOT. PRESENT(ylabel)) then
            t_ylabel = 'y'
        else
            t_ylabel = ylabel
        end if

        if (t_grid) then
            set_grid = "set grid"
        else
            set_grid = ""
        end if

        config = "#"//OUTP_FNAME(t_fname)//ENDL// &
        'set size 1,1'//ENDL// &
        'set terminal pdf'//ENDL// &
        'set output "'//OUTP_FNAME(t_fname)//'"'//ENDL// &
        'set title "'//t_title//'"'//ENDL// &
        !'set nokey'//ENDL// &
        set_grid//ENDL// &
        'set xlabel "'//t_xlabel//'"'//ENDL// &
        'set ylabel "'//t_ylabel//'"'//ENDL// &
        'xydata = "'//DATA_FNAME(t_fname)//'"'//ENDL

!       ===============================================================================
10      format(A, ' ')
        open (newunit=file, action='write', file=PLOT_FNAME(t_fname), status='replace')
            write(file, *) config
            write(file, 10, advance='no') 'plot'
            do i=1, m
                if (i == 1) then
                    write(file, 10, advance='no') 'xydata u 1:2 t "'//t_legend(1)%str//'" '//t_with//PLOT_ENDL
                else if (i == m) then
                    write(file, *) '"" u 1:'//STR(m + 1)//' t "'//t_legend(m)%str//'" '//t_with
                else
                    write(file, 10, advance='no') '"" u 1:'//STR(i + 1)//' t "'//t_legend(i)%str//'" '//t_with//PLOT_ENDL
                end if
            end do
        close(file)
    end subroutine 

    subroutine xy_multiplot(x, y, n, m, fname, title, xlabel, ylabel, grid, with, legend)
        implicit none
        integer :: file, i, j, m
        integer, intent(in) :: n
        double precision, dimension(n), intent(in) :: x
        double precision, dimension(n, m), intent(in) :: y

        logical, optional :: grid
        logical :: t_grid

        character(len=*), optional :: with
        character(len=:), allocatable :: t_with

        character(len=*), optional :: title
        character(len=:), allocatable :: t_title

        character(len=*), optional :: ylabel
        character(len=:), allocatable :: t_ylabel

        character(len=*), optional :: xlabel
        character(len=:), allocatable :: t_xlabel

        character(len=*), optional :: fname
        character(len=:), allocatable :: t_fname

        type(StringArray), dimension(:), optional :: legend
        type(StringArray), dimension(m) :: t_legend

        if (.NOT. PRESENT(legend)) then
            do i=1,m
                t_legend(i)%str = STR(i)
            end do
        else
            t_legend(:) = legend(:)
        end if

        if (.NOT. PRESENT(fname)) then
            t_fname = DEFAULT_FNAME
        else
            t_fname = fname
        end if

        if (.NOT. PRESENT(grid)) then
            t_grid = .TRUE.
        else
            t_grid = grid
        end if

        if (.NOT. PRESENT(with)) then
            t_with = "lines"
        else
            t_with = with
        end if

        if (.NOT. PRESENT(title)) then
            t_title = ''
        else
            t_title = title
        end if

        if (.NOT. PRESENT(xlabel)) then
            t_xlabel = 'x'
        else
            t_xlabel = xlabel
        end if

        if (.NOT. PRESENT(ylabel)) then
            t_ylabel = 'y'
        else
            t_ylabel = ylabel
        end if

!       ===================== Write to data file ========================
10      format(F24.12, ' ')
        open (newunit=file, action='write', file=DATA_FNAME(t_fname), status='replace')
        do i=1,n
            write(file, 10, advance='no') x(i)
            do j=1,m
                write(file, 10, advance='no') y(i, j)
            end do
            write(file, *) ''
        end do
        close(file)
!       =================================================================
        call multiplot_config(m, fname=t_fname, title=t_title, xlabel=t_xlabel, ylabel=t_ylabel, grid=t_grid, with=t_with, &
        legend=t_legend)
        call EXECUTE_COMMAND_LINE(GNU_PLOT_CMD(t_fname))
    end subroutine
end module Plot