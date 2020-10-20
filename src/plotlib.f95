module PlotLib
    use Util
    implicit none
    character(len=*), parameter :: DEFAULT_FNAME = 'plotfile'
    character(len=*), parameter :: PLOT_ENDL = ',\'//ENDL

    logical :: g_INPLOT = .FALSE.
    logical :: g_INMULTIPLOT = .FALSE.
    
    character(len=:), allocatable :: g_FNAME
    character(len=:), allocatable :: g_OUTP_FNAME
    character(len=:), allocatable :: g_PLOT_FNAME

    integer :: g_M, g_N

    type SPLOT
        integer :: i, j
        integer :: n = 0
        logical :: grid = .FALSE.
        logical :: done = .FALSE.
        character(len=:), allocatable :: title, xlabel, ylabel
        type(StringArray), dimension(:), allocatable :: legend, with
    end type

    type(SPLOT), dimension(:, :), allocatable :: g_SUBPLOTS

    contains

    function REMOVE_TEMP_FILES(fname) result (cmd)
        implicit none
        character(len=*) :: fname
        character(len=:), allocatable :: plt
        character(len=:), allocatable :: dat
        character(len=:), allocatable :: cmd

        plt = 'plot/'//fname//'*.plt'
        dat = 'plot/'//fname//'*.dat'

        cmd = 'rm '//plt//' '//dat
        return
    end function

    function PLOT_FNAME(fname) result (path)
        implicit none
        character(len=*) :: fname
        character(len=:), allocatable :: path
        path = 'plot/'//fname//'.plt'
        return
    end function

    function DATA_FNAME(fname, i, j, n) result (path)
        implicit none
        integer, optional, intent(in) :: i, j, n
        character(len=*) :: fname
        character(len=:), allocatable :: path, t_i, t_j, t_n

        if (.NOT. PRESENT(i)) then
            t_i = '1'
        else
            t_i = STR(i)
        end if 

        if (.NOT. PRESENT(j)) then
            t_j = '1'
        else
            t_j = STR(j)
        end if

        if (.NOT. PRESENT(n)) then
            t_n = '1'
        else
            t_n = STR(n)
        end if

        path = 'plot/'//fname//'_'//t_i//'_'//t_j//'_'//t_n//'.dat'
        return
    end function

    function OUTP_FNAME(fname) result (path)
        implicit none
        character(len=*) :: fname
        character(len=:), allocatable :: path
        path = 'plot/'//fname//'.pdf'
        return
    end function

    function GNU_PLOT_CMD(fname) result (cmd)
        implicit none
        character(len=*) :: fname
        character(len=:), allocatable :: cmd
        cmd = 'gnuplot -p '//PLOT_FNAME(fname)
    end function
    
    subroutine subplot_config(i, j, title, xlabel, ylabel, legend, with, grid)
        integer, intent(in) :: i, j
        integer :: k, n

        logical, optional :: grid

        character(len=*), optional :: title, ylabel, xlabel

        type(StringArray), dimension(:), optional :: legend, with

        type(SPLOT) :: subplot_ij
        
        subplot_ij = g_SUBPLOTS(i, j)

        if (subplot_ij%done) then
            call error("Duplicate configuration of subplot ("//STR(i)//", "//STR(j)//")")
            stop "ERROR"
        end if

        n = subplot_ij%n

        allocate(subplot_ij%legend(n))

        if (.NOT. PRESENT(legend)) then
            do k=1,n
                subplot_ij%legend(k)%str = 't '//quote(STR(i))
            end do
        else
            do k=1,n
                subplot_ij%legend(k)%str = 't '//quote(legend(i)%str)
            end do
        end if

        if (.NOT. PRESENT(with)) then
            do k=1,n
                subplot_ij%with(k)%str = 'w lines'
            end do
        else
            do k=1,n
                subplot_ij%with(k)%str = 'w '//with(k)%str
            end do
        end if

        if (.NOT. PRESENT(grid)) then
            subplot_ij%grid = .TRUE.
        else
            subplot_ij%grid = grid
        end if

        if (.NOT. PRESENT(title)) then
            subplot_ij%title = ''
        else
            subplot_ij%title = title
        end if

        if (.NOT. PRESENT(xlabel)) then
            subplot_ij%xlabel = 'x'
        else
            subplot_ij%xlabel = xlabel
        end if

        if (.NOT. PRESENT(ylabel)) then
            subplot_ij%ylabel = 'y'
        else
            subplot_ij%ylabel = ylabel
        end if

        subplot_ij%done = .TRUE.

    end subroutine

    subroutine subplot(i, j, x, y, n)
        integer :: file, k
        integer, intent(in) :: i, j, n
        double precision, dimension(n), intent(in) :: x
        double precision, dimension(n), intent(in) :: y

        character(len=:), allocatable :: s_data_fname

        type(SPLOT) :: subplot_ij
        
        subplot_ij = g_SUBPLOTS(i, j)

        if (subplot_ij%done) then
            call error("Plot over finished subplot ("//STR(i)//", "//STR(j)//")")
            stop "ERROR"
        else
            subplot_ij%n = subplot_ij%n + 1
        end if

        s_data_fname = DATA_FNAME(g_FNAME, i, j, subplot_ij%n)

!       ==================== Write to Plot File ==============================================
10      format(F16.8, ' ')
        open(newunit=file, file=s_data_fname, status="new", position="append", action="write")
        write(file, *)
        do k=1, n
            write(file, 10, advance='no') x(k)
            write(file, *) y(k)
        end do
        close(file)
!       ======================================================================================
    end subroutine

!   ======= Pipeline ============
    subroutine begin_plot(fname)
        integer :: file
        character(len=*), optional :: fname; character(len=:), allocatable :: t_fname;

        if (.NOT. PRESENT(fname)) then
            t_fname = DEFAULT_FNAME
        else
            t_fname = fname
        end if

        g_PLOT_FNAME = g_FNAME
        g_OUTP_FNAME = OUTP_FNAME(g_FNAME)

        open(newunit=file, file=g_PLOT_FNAME, status="new", action="write")
        write(file, *) 'set terminal pdf'
        write(file, *) 'set output '//quote(g_OUTP_FNAME)
        close(file)

        g_INPLOT = .TRUE.

    end subroutine

    subroutine subplots(m, n)
        integer, optional, intent(in) :: m, n
        integer :: t_m, t_n

        if ((.NOT. PRESENT(m)) .OR. (m <= 0)) then
            t_m = 1
        else
            t_m = m
        end if

        if ((.NOT. PRESENT(n)) .OR. (n <= 0)) then
            t_n = 1
        else
            t_n = n
        end if

        if (.NOT. g_INPLOT) then
            call begin_plot()
        end if

!       ===== Allocate Variables =====
        allocate(g_SUBPLOTS(t_m, t_n))
        g_M = t_m
        g_N = t_n
        g_INMULTIPLOT = .TRUE.
!       ==============================
    end subroutine

    subroutine render_plot()
        integer :: file, i, j, k, m, n
        type(SPLOT) :: subplot_ij

!       === Check Plot =========
        if (.NOT. g_INPLOT) then
            call error("No active plot to render.")
            stop "ERROR"
        end if
!       ========================

        m = g_M
        n = g_N

!       ==================== Write to Plot File ==============================================      
        open(newunit=file, file=g_PLOT_FNAME, status="old", position="append", action="write")
        write(file, *) 'set size 1,1;'
        write(file, *) 'set origin 0,0;'

        if (g_INMULTIPLOT) then
            write(file, *) 'set multiplot layout '//STR(m)//','//STR(n)//' rowsfirst;'        
        end if

10      format(A, ' ')
!       =========== Plot data ============
        do i= 1, m
            do j = 1, n
                subplot_ij = g_SUBPLOTS(i, j)

                write(file, *) 'set title '//quote(subplot_ij%title)//';'
                write(file, *) 'set xlabel '//quote(subplot_ij%xlabel)//';'
                write(file, *) 'set ylabel '//quote(subplot_ij%ylabel)//';'

                if (subplot_ij%grid) then
                    write(file, *) 'set grid;'
                else
                    write(file, *) 'unset grid;'
                end if

                write(file, 10, advance='no') 'plot'

                do k = 1, subplot_ij%n
                    write(file, 10, advance='no') quote(DATA_FNAME(g_FNAME, i, j, k))
                    write(file, 10, advance='no') 'u 1:2'
                    write(file, 10, advance='no') subplot_ij%legend(k)%str
                    write(file, 10, advance='no') subplot_ij%with(k)%str
                end do

                if (k == (subplot_ij%n)) then
                    write(file, *) ';'
                else
                    write(file, *) ',\'
                end if
            end do
        end do
!       ===================================

!       == Finish Multiplot ===        
        if (g_INMULTIPLOT) then
            write(file, *) 'unset multiplot'
            g_INMULTIPLOT = .FALSE.
        end if
!       =======================        
        close(file)
!       ======================================================================================

!       ===== Call GNUPLOT and remove temporary files =======
!        call EXECUTE_COMMAND_LINE(GNU_PLOT_CMD(g_PLOT_FNAME))

!        call EXECUTE_COMMAND_LINE(REMOVE_TEMP_FILES(g_FNAME))
!       =====================================================

!       ========= Free Variables ======================
        deallocate(g_FNAME, g_OUTP_FNAME, g_PLOT_FNAME)
        deallocate(g_SUBPLOTS)
        g_INPLOT = .FALSE.
!       ===============================================
    end subroutine
end module PlotLib