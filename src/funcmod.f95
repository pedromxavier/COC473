!   Func Module

    module Func
        implicit none

        double precision, parameter :: g = 9.80600D0
        double precision, parameter :: k = 0.00341D0

        double precision, parameter :: t1 = 0.0D0
        double precision, parameter :: t2 = 0.0D0

        character (len = *), parameter :: F1_NAME = "f(x) = log(cosh(x * sqrt(g * j))) - 50"
        character (len = *), parameter :: F2_NAME = "f(x) = 4 * cos(x) - exp(2 * x)"
        character (len = *), parameter :: F3_NAME = "f(x, y, z) := \n\t"// &
                                                    "16x^4 + 16y^4 + z^4 = 16\n\"// &
                                                    "tx^2 + y^2 + x^2 = 3\n\t"// &
                                                    "x^3 - y + z = 1"

    contains

    function f1(x) result (y)
        implicit none
        double precision :: x, y
        y = DLOG(DCOSH(x * DSQRT(g * k))) - 50.0D0
        return
    end function

    function df1(x) result (y)
        implicit none
        double precision :: x, y
        y = (DSINH(x * DSQRT(g * k)) * DSQRT(g * k)) / DCOSH(x * DSQRT(g * k))
        return
    end function
    
    function f2(x) result (y)
        implicit none
        double precision :: x, y
        y = 4 * DCOS(x) - DEXP(2 * x)
        return
    end function

    function df2(x) result (y)
        implicit none
        double precision :: x, y
        y = - 4 * DSIN(x) - 2 * DEXP(2 * x)
        return
    end function

    function f3(x) result (y)
        implicit none
        double precision :: x, y
        y = x * x
        return
    end function

    function df3(x) result (y)
        implicit none
        double precision :: x, y
        y = 2 * x
        return
    end function
    
!   ======= R^n -> R^n functions =======
    function ff1_1(x) result (y)
!       R^3 -> R
        implicit none
        double precision :: x(3)
        double precision :: y
        
!       16 x^4 + 16 y^4 + z^4 = 16
        y = (16 * x(1) ** 4 + 16 * x(2) ** 4 + x(3) ** 4) - 16.0D0
        return
    end function

    function ff1_2(x) result (y)
!       R^3 -> R
        implicit none
        double precision :: x(3)
        double precision :: y
        
!       x^2 + y^2 + z^2 = 3
        y = x(1) ** 2 + x(2) ** 2 +x(3) ** 2 - 3.0D0
    end function

    function ff1_3(x) result (y)
!       R^3 -> R
        implicit none
        double precision :: x(3)
        double precision :: y
        
!       x^3 - y + z = 1
        y = x(1) ** 3 - x(2) + x(3) - 1.0D0
        return
    end function

    function ff1(x, n) result (y)
!       R^3 -> R^3 (n == 3)
        implicit none
        integer :: n
        double precision :: x(n), y(n)

        y(:) = (/ ff1_1(x), ff1_2(x), ff1_3(x) /)
        return
    end function

!   ========== Derivatives ===========
    function dff1_1_1(x) result (y)
!       R^3 -> R
        implicit none
        double precision :: x(3)
        double precision :: y
        
!       16 x^4 + 16 y^4 + z^4 = 16
        y = 64 * x(1) ** 3
        return
    end function

    function dff1_2_1(x) result (y)
!       R^3 -> R
        implicit none
        double precision :: x(3)
        double precision :: y
        
!       x^2 + y^2 + z^2 = 3
        y = 2 * x(1)
    end function

    function dff1_3_1(x) result (y)
!       R^3 -> R
        implicit none
        double precision :: x(3)
        double precision :: y
        
!       x^3 - y + z = 1
        y = 3 * x(1) ** 2
        return
    end function

    function dff1_1_2(x) result (y)
!       R^3 -> R
        implicit none
        double precision :: x(3)
        double precision :: y
        
!       16 x^4 + 16 y^4 + z^4 = 16
        y=  64 * x(2) ** 3
        return
    end function

    function dff1_2_2(x) result (y)
!       R^3 -> R
        implicit none
        double precision :: x(3)
        double precision :: y
        
!       x^2 + y^2 + z^2 = 3
        y = 2 * x(2)
    end function

    function dff1_3_2(x) result (y)
!       R^3 -> R
        implicit none
        double precision :: x(3)
        double precision :: y
        
!       x^3 - y + z = 1
        y = - 1.0D0
        return
    end function

    function dff1_1_3(x) result (y)
!       R^3 -> R
        implicit none
        double precision :: x(3)
        double precision :: y
        
!       16 x^4 + 16 y^4 + z^4 = 16
        y= 4 * x(3) ** 3
        return
    end function

    function dff1_2_3(x) result (y)
!       R^3 -> R
        implicit none
        double precision :: x(3)
        double precision :: y
        
!       x^2 + y^2 + z^2 = 3
        y = 2 * x(3)
    end function

    function dff1_3_3(x) result (y)
!       R^3 -> R
        implicit none
        double precision :: x(3)
        double precision :: y
        
!       x^3 - y + z = 1
        y = 1.0D0
        return
    end function

    function dff1(x, n) result (J)
!       R^3 -> R^3x3 (n == 3)
        implicit none
        integer :: n
        double precision :: x(n), J(n, n)

        J(1, :) = (/dff1_1_1(x), dff1_1_2(x), dff1_1_3(x)/)
        J(2, :) = (/dff1_2_1(x), dff1_2_2(x), dff1_2_3(x)/)
        J(3, :) = (/dff1_3_1(x), dff1_3_2(x), dff1_3_3(x)/)
        return
    end function

!   ================== Another function =====================
    function ff2_1(x) result (y)
        implicit none
        double precision, dimension(3) :: x
        double precision :: y
        y = x(1)**2 + 2*x(2)**2 + 6*x(3)**2 - 1.0D0
        return
    end function

    function ff2_2(x) result (y)
        implicit none
        double precision, dimension(3) :: x
        double precision :: y
        y = 2*x(2)*( &
            3*x(1)**2 + &
            4*x(2)**2 + &
            18*x(1)*x(3) + &
            54*x(3)**4 &
            ) - t1
        return
    end function

    function ff2_3(x) result (y)
        implicit none
        double precision, dimension(3) :: x
        double precision :: y
        y = 3*( &
            x(1) + &
            20*x(1)**2*x(2)**2 + &
            20*x(2)**4 + &
            8*x(1)*(x(1)**2 + 24*x(2)**2)*x(3) + &
            12*(7*x(1)**2 + &
            62*x(2)**2)*x(3)**2 + &
            432*x(1)*x(3)**3 + &
            1116*x(3)**4 &
            ) - t2
        return
    end function

    function ff2(x, n) result (y)
        implicit none
        integer, intent(in) :: n
        double precision, dimension(n) :: x, y
        y = (/ff2_1(x), ff2_2(x), ff2_3(x)/)
        return
    end function

!   =============== Derivatives ==================

!   ========== Derivatives ===========
    function dff2_1_1(x) result (y)
!       R^3 -> R
        implicit none
        double precision :: x(3)
        double precision :: y
        y = 2*x(1)
        return
    end function

    function dff2_2_1(x) result (y)
!       R^3 -> R
        implicit none
        double precision :: x(3)
        double precision :: y
        
!       x^2 + y^2 + z^2 = 3
        y = 12*x(1)*x(2) + &
            36*x(2)*x(3)
    end function

    function dff2_3_1(x) result (y)
!       R^3 -> R
        implicit none
        double precision :: x(3)
        double precision :: y
        y = 3 + &
            120*x(1)*x(2)**2 + &
            72*x(1)**2*x(3) + &
            576*x(2)**2*x(3) + &
            504*x(1)*x(3)**2 + &
            1296*x(3)**3
        return
    end function

    function dff2_1_2(x) result (y)
!       R^3 -> R
        implicit none
        double precision :: x(3)
        double precision :: y
        y = 4*x(2)
        return
    end function

    function dff2_2_2(x) result (y)
!       R^3 -> R
        implicit none
        double precision :: x(3)
        double precision :: y
        y = 6*x(1)**2 + &
            24*x(2)**2 + &
            36*x(1)*x(3) + &
            108*x(3)**4
    end function

    function dff2_3_2(x) result (y)
!       R^3 -> R
        implicit none
        double precision :: x(3)
        double precision :: y
        y = 120*x(1)**2*x(2) + &
            240*x(2)**3 + &
            1152*x(1)*x(2)*x(3) + &
            4464*x(2)*x(3)**2
        return
    end function

    function dff2_1_3(x) result (y)
!       R^3 -> R
        implicit none
        double precision :: x(3)
        double precision :: y
        y = 12*x(3)
        return
    end function

    function dff2_2_3(x) result (y)
!       R^3 -> R
        implicit none
        double precision :: x(3)
        double precision :: y
        y = 36*x(1)*x(2) + &
            432*x(2)*x(3)**3
    end function

    function dff2_3_3(x) result (y)
!       R^3 -> R
        implicit none
        double precision :: x(3)
        double precision :: y
        y = 24*x(1)**3 + &
            576*x(1)*x(2)**2 + &
            504*x(1)**2*x(3) + &
            4464*x(2)**2*x(3) + &
            3888*x(1)*x(3)**2 + &
            13392*x(3)**3
        return
    end function
    
    function dff2(x, n) result (J)
!       R^3 -> R^3x3 (n == 3)
        implicit none
        integer :: n
        double precision :: x(n), J(n, n)

        J(1, :) = (/dff2_1_1(x), dff2_1_2(x), dff2_1_3(x)/)
        J(2, :) = (/dff2_2_1(x), dff2_2_2(x), dff2_2_3(x)/)
        J(3, :) = (/dff2_3_1(x), dff2_3_2(x), dff2_3_3(x)/)
        return
    end function


    end module Func
