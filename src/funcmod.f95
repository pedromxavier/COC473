!   Func Module

    module Func
        use Util
        implicit none
        double precision :: F1_G = 9.80600D0
        double precision :: F1_K = 0.00341D0

!       >> F13 <<        
        double precision :: F13_A = 0.0D0
        double precision :: F13_B = 10.0D0
        double precision :: F13_DT = 1.0D-1

!       >> F14 <<
        character (len = *), parameter :: F14_NAME = "m y''(t) = k y(t) "
        double precision :: F14_M = 1.0D0
        double precision :: F14_C = 0.2D0
        double precision :: F14_K = 1.0D0
        double precision :: F14_W = 0.5D0
        double precision :: F14_Y0 = 0.0D0
        double precision :: F14_DY0 = 0.0D0
        double precision :: F14_A = 0.0D0
        double precision :: F14_B = 100.0D0
        double precision :: F14_DT = 1.0D-1

        double precision :: F15_G = 9.80600D0
        double precision :: F15_KD = 1.0D0
        double precision :: F15_Y0 = 0.0D0
        double precision :: F15_DY0 = 0.0D0
        double precision :: F15_DT = 1.0D-1

        double precision :: t1 = 0.0D0
        double precision :: t2 = 0.0D0

        double precision :: wn = 1.00D0
        double precision :: xi = 0.05D0
        double precision :: Hs = 3.0D0
        double precision :: Tz = 5.0D0

        character (len = *), parameter :: F1_NAME = "f(x) = log(cosh(x * √(g * j))) - 50"
        character (len = *), parameter :: F2_NAME = "f(x) = 4 * cos(x) - exp(2 * x)"
        character (len = *), parameter :: F3_NAME = "f(x, y, z) :="//ENDL//TAB// &
        "16x⁴ + 16y⁴ + z⁴ = 16"//ENDL//TAB// &
        "x² + y² + x² = 3"//ENDL//TAB// &
        "x³ - y + z = 1"
        character (len = *), parameter :: F4_NAME = "f(c2, c3, c4) :="//ENDL//TAB// &
        "c2² + 2 c3² + 6 c4² = 1"//ENDL//TAB// &
        "8 c3³ + 6 c3 c2² + 36 c3 c2 c4 + 108 c3 c4⁴ = t1"//ENDL//TAB// &
        "60 * c3⁴ + 60 * c3² * c2² + 576 * c3² * c2 * c4 + "// & 
        "2232 * c3²  * c4² + 252 * c4² * c2² + "// &
        "1296 * c4³ c2 + 3348 c4⁴ + 24 c2³ c4 + 3 c2 = t2"
        character (len = *), parameter :: F5_NAME = "f(x) = b1 + b2 x^b3"

        character (len = *), parameter :: F6_NAME = "f(x) = exp(-x²/2) / √(2 π)"

        character (len = *), parameter :: F78_NAME = "Sσ(ω) = RAO(ω)²  Sη(ω)"//ENDL//TAB// &
        "RAO(ω) = 1 / √((1 - (ω/ωn)²)² + (2ξω/ωn)²)"

        character (len = *), parameter :: F7a_NAME = "f(ω) = Sσ(ω)"//ENDL//TAB// &
        "Sη(ω) = 2"
        character (len = *), parameter :: F7b_NAME = "f(ω) = ω² Sσ(ω)"//ENDL//TAB// &
        "Sη(ω) = 2"
        character (len = *), parameter :: F8a_NAME = "f(ω) = Sσ(ω)"//ENDL//TAB// &
        "Sη(ω) = ((4 π³ Hs²) / (ω⁵ Tz⁴)) exp(-(16 π³) / (ω⁴ Tz⁴))"
        character (len = *), parameter :: F8b_NAME = "f(ω) = ω² Sσ(ω)"//ENDL//TAB// &
        "Sη(ω) = ((4 π³ Hs²) / (ω⁵ Tz⁴)) exp(-(16 π³) / (ω⁴ Tz⁴))"

        character (len = *), parameter :: F9_NAME = "f(x) = 2 + 2x - x² + 3x³"

        character (len = *), parameter :: F10_NAME = "f(x) = 1 / (1 + x²)"

        character (len = *), parameter :: F11_NAME = "f(x) = exp(- x²/2) / √(2 π)"
        character (len = *), parameter :: F12_NAME = "f(x) = x² exp(- x²/2) / √(2 π)"

        character (len = *), parameter :: F13_NAME = "y'(t) = -2 t y(t)²"//ENDL//"y(0) = 1"
    contains

    function f1(x) result (y)
        implicit none
        double precision :: x, y
        y = DLOG(DCOSH(x * DSQRT(F1_G * F1_K))) - 50.0D0
        return
    end function

    function df1(x) result (y)
        implicit none
        double precision :: x, y
        y = (DSINH(x * DSQRT(F1_G * F1_K)) * DSQRT(F1_G * F1_K)) / DCOSH(x * DSQRT(F1_G * F1_K))
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
!       R³ -> R
        implicit none
        double precision :: x(3)
        double precision :: y
        
!       16 x⁴ + 16 y⁴ + z⁴ = 16
        y = (16 * x(1) ** 4 + 16 * x(2) ** 4 + x(3) ** 4) - 16.0D0
        return
    end function

    function ff1_2(x) result (y)
!       R³ -> R
        implicit none
        double precision :: x(3)
        double precision :: y
        
!       x² + y² + z² = 3
        y = x(1) ** 2 + x(2) ** 2 +x(3) ** 2 - 3.0D0
    end function

    function ff1_3(x) result (y)
!       R³ -> R
        implicit none
        double precision :: x(3)
        double precision :: y
        
!       x³ - y + z = 1
        y = x(1) ** 3 - x(2) + x(3) - 1.0D0
        return
    end function

    function ff1(x, n) result (y)
!       R³ -> R³ (n == 3)
        implicit none
        integer :: n
        double precision :: x(n), y(n)

        y(:) = (/ ff1_1(x), ff1_2(x), ff1_3(x) /)
        return
    end function

!   ========== Derivatives ===========
    function dff1_1_1(x) result (y)
!       R³ -> R
        implicit none
        double precision :: x(3)
        double precision :: y
        
!       16 x⁴ + 16 y⁴ + z⁴ = 16
        y = 64 * x(1) ** 3
        return
    end function

    function dff1_2_1(x) result (y)
!       R³ -> R
        implicit none
        double precision :: x(3)
        double precision :: y
        
!       x² + y² + z² = 3
        y = 2 * x(1)
    end function

    function dff1_3_1(x) result (y)
!       R³ -> R
        implicit none
        double precision :: x(3)
        double precision :: y
        
!       x³ - y + z = 1
        y = 3 * x(1) ** 2
        return
    end function

    function dff1_1_2(x) result (y)
!       R³ -> R
        implicit none
        double precision :: x(3)
        double precision :: y
        
!       16 x⁴ + 16 y⁴ + z⁴ = 16
        y=  64 * x(2) ** 3
        return
    end function

    function dff1_2_2(x) result (y)
!       R³ -> R
        implicit none
        double precision :: x(3)
        double precision :: y
        
!       x² + y² + z² = 3
        y = 2 * x(2)
    end function

    function dff1_3_2(x) result (y)
!       R³ -> R
        implicit none
        double precision :: x(3)
        double precision :: y
!       x³ - y + z = 1
        y = - 1.0D0
        return
    end function

    function dff1_1_3(x) result (y)
!       R³ -> R
        implicit none
        double precision :: x(3)
        double precision :: y
        
!       16 x⁴ + 16 y⁴ + z⁴ = 16
        y= 4 * x(3) ** 3
        return
    end function

    function dff1_2_3(x) result (y)
!       R³ -> R
        implicit none
        double precision :: x(3)
        double precision :: y
        
!       x² + y² + z² = 3
        y = 2 * x(3)
    end function

    function dff1_3_3(x) result (y)
!       R³ -> R
        implicit none
        double precision :: x(3)
        double precision :: y
        
!       x³ - y + z = 1
        y = 1.0D0
        return
    end function

    function dff1(x, n) result (J)
!       R³ -> R³x3 (n == 3)
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
        integer :: n
        double precision, dimension(n) :: x, y
        y = (/ff2_1(x), ff2_2(x), ff2_3(x)/)
        return
    end function

!   =============== Derivatives ==================

!   ========== Derivatives ===========
    function dff2_1_1(x) result (y)
!       R³ -> R
        implicit none
        double precision :: x(3)
        double precision :: y
        y = 2*x(1)
        return
    end function

    function dff2_2_1(x) result (y)
!       R³ -> R
        implicit none
        double precision :: x(3)
        double precision :: y
        
!       x² + y² + z² = 3
        y = 12*x(1)*x(2) + &
            36*x(2)*x(3)
    end function

    function dff2_3_1(x) result (y)
!       R³ -> R
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
!       R³ -> R
        implicit none
        double precision :: x(3)
        double precision :: y
        y = 4*x(2)
        return
    end function

    function dff2_2_2(x) result (y)
!       R³ -> R
        implicit none
        double precision :: x(3)
        double precision :: y
        y = 6*x(1)**2 + &
            24*x(2)**2 + &
            36*x(1)*x(3) + &
            108*x(3)**4
    end function

    function dff2_3_2(x) result (y)
!       R³ -> R
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
!       R³ -> R
        implicit none
        double precision :: x(3)
        double precision :: y
        y = 12*x(3)
        return
    end function

    function dff2_2_3(x) result (y)
!       R³ -> R
        implicit none
        double precision :: x(3)
        double precision :: y
        y = 36*x(1)*x(2) + &
            432*x(2)*x(3)**3
    end function

    function dff2_3_3(x) result (y)
!       R³ -> R
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
!       R³ -> R³x3 (n == 3)
        implicit none
        integer :: n
        double precision :: x(n), J(n, n)

        J(1, :) = (/dff2_1_1(x), dff2_1_2(x), dff2_1_3(x)/)
        J(2, :) = (/dff2_2_1(x), dff2_2_2(x), dff2_2_3(x)/)
        J(3, :) = (/dff2_3_1(x), dff2_3_2(x), dff2_3_3(x)/)
        return
    end function

!   ============ One more function =============
    function ff5(x, b, m, n) result (z)
        implicit none
        integer :: m, n
        double precision, dimension(m), intent(in) :: b
        double precision, dimension(n), intent(in) :: x
        double precision, dimension(n) :: z
        integer :: i
!       m == 3
        do i=1, n
            z(i) = b(1) + (b(2) * (x(i) ** b(3)))
        end do
        return
    end function

!   ========= Derivatives ==========
    function dff5(x, b, m, n) result (J)
        implicit none
        integer :: m, n
        double precision, dimension(m), intent(in) :: b
        double precision, dimension(n), intent(in) :: x
        double precision, dimension(n, m) :: J
        integer :: i
!       m == 3
        do i=1, n
            J(i, :) = (/ &
!               df_i / db_1
                1.0D0, &
!               df_i / db_2
                x(i) ** b(3), &
!               df_i / db_3
                b(2) * b(3) * (x(i) ** (b(3) - 1.0D0)) &
                /)
        end do
        return
    end function

!   ======== Function 6 ============
    function f6(x) result (y)
        implicit none
        double precision :: x, y

        y = DEXP(-(x*x)/2) / DSQRT(2 * PI)
        return
    end function

!   ======== Functions 7 & 8 ============
    function RAO(w) result (z)
        implicit none
        double precision :: w, z

        z = 1.0 / DSQRT((1.0D0 - (w/wn) ** 2) ** 2 + (2 * xi * (w/wn)) ** 2)
    end function

    function Sn1(w) result (z)
        implicit none
        double precision :: w, z

        z = 2.0D0
        return
    end function

    function Sn2(w) result (z)
        implicit none
        double precision :: w, z

        z = (4 * (PI ** 3) * (Hs ** 2)) * DEXP(- (16.0D0 * (PI ** 3)) / ((w ** 4) * (Tz ** 4))) / (w ** 5) * (Tz ** 4)
        return
    end function

    function Ss(w, Sn) result (z)
        implicit none
        double precision :: w, z
        interface
            function Sn(w) result (z)
                implicit none
                double precision :: w, z
            end function
        end interface
        z = (RAO(w) ** 2) * Sn(w)
        return
    end function

    function f7a(w) result (z)
        implicit none
        double precision :: w, z
        z = Ss(w, Sn1)
        return
    end function

    function f7b(w) result (z)
        implicit none
        double precision :: w, z
        z = (w ** 2) * Ss(w, Sn1)
        return
    end function

    function f8a(w) result (z)
        implicit none
        double precision :: w, z
        z = Ss(w, Sn1)
        return
    end function

    function f8b(w) result (z)
        implicit none
        double precision :: w, z
        z = (w ** 2) * Ss(w, Sn1)
        return
    end function

!   ========== Function 9 ==============
    function f9(x) result (y)
        implicit none
        double precision :: x, y
        y = 2.0D0 + 2.0D0 * x - x ** 2 + 3.0D0 * x ** 3
        return
    end function

!   ========== Function 10 ==============
    function f10(x) result (y)
        implicit none
        double precision :: x, y
        y = 1.0D0 / (1.0D0 + x ** 2)
        return
    end function

!   ========== Function 11 ==============
    function f11a(x) result (y)
        implicit none
        double precision :: x, y
        y = DEXP((x ** 2) / 2) / DSQRT(8 * PI)
        return
    end function

    function f11b(x) result (y)
        implicit none
        double precision :: x, y
        y = DEXP(-(x ** 2) / 2) / DSQRT(8 * PI)
        return
    end function

!   ========== Function 12 ==============
    function f12(x) result (y)
        implicit none
        double precision :: x, y
        y = (x ** 2) * DEXP((x ** 2) / 2) / DSQRT(2 * PI)
        return
    end function


!   ========== Function 13 ==============
    function df13(t, y) result (u)
        implicit none
        double precision :: t, y, u
        u = - 2 * t * (y ** 2)
        return
    end function

    function f13(t) result (y)
        implicit none
        double precision :: t, y
        y = 1 / (1 + (t**2))
        return
    end function

!   ========== Function 14 ===============
    function F14_F(t) result (y)
        implicit none
        double precision :: t, y
        y = 2 * DSIN(F14_W * t) + DSIN(2 * F14_W * t) + DCOS(3 * F14_W * t)
        return
    end function

    function d2f14(t, y, dy) result (u)
        implicit none
        double precision :: t, y, dy, u
        u = (F14_F(t) - F14_K * y - F14_C * dy) / F14_M
        return
    end function

!   ========== Function 15 ===============
    function d2f15(t, y, dy) result (u)
        implicit none
        double precision :: t, y, dy, u
        u = - F15_G - F15_KD * dy * DABS(dy)
        return
    end function

    end module Func