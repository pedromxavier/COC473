!   Func Module

    module Func
        use Util
        implicit none

!       >> F1 <<
        character (len = *), parameter :: F1_NAME = "f(x) = log(cosh(x * √(g * j))) - 50"
        double precision :: F1_G = 9.80600D0
        double precision :: F1_K = 0.00341D0

!       >> F2 <<
        character (len = *), parameter :: F2_NAME = "f(x) = 4 * cos(x) - exp(2 * x)"

!       >> F3 <<
        character (len = *), parameter :: F3_NAME = "f(x, y, z) :="//ENDL// &
        "16x⁴ + 16y⁴ + z⁴ = 16"//ENDL// &
        "x² + y² + x² = 3"//ENDL// &
        "x³ - y + z = 1"
        integer :: F3_N = 3

!       >> F4 <<
        character (len = *), parameter :: F4_NAME = "f(c2, c3, c4) :="//ENDL// &
        "c2² + 2 c3² + 6 c4² = 1"//ENDL// &
        "8 c3³ + 6 c3 c2² + 36 c3 c2 c4 + 108 c3 c4⁴ = θ1"//ENDL// &
        "60 * c3⁴ + 60 * c3² * c2² + 576 * c3² * c2 * c4 + "// & 
        "2232 * c3² * c4² + 252 * c4² * c2² + "// &
        "1296 * c4³ c2 + 3348 c4⁴ + 24 c2³ c4 + 3 c2 = θ2"
        double precision :: F4_TT1(3) = (/ 0.0D0, 0.75D0,  0.000D0 /)
        double precision :: F4_TT2(3) = (/ 3.0D0, 6.50D0, 11.667D0 /)
        double precision :: F4_T1 = 0.0D0
        double precision :: F4_T2 = 0.0D0
        integer :: F4_N = 3

!       >> F5 <<
        character (len = *), parameter :: F5_NAME = "f(x) = b1 + b2 x^b3"
        integer :: F5_N = 3
        
!       >> F6 <<
        character (len = *), parameter :: F6_NAME = "f(x) = exp(-x²/2) / √(2 π)"

!       >> F7 <<
        character (len = *), parameter :: F7_NAME = "Sσ(ω) = RAO(ω)²  Sη(ω)"//ENDL//TAB// &
        "RAO(ω) = 1 / √((1 - (ω/ωn)²)² + (2ξω/ωn)²)"

        character (len = *), parameter :: F7a_NAME = "Sη(ω) = 2"
        character (len = *), parameter :: F7b_NAME = "Sη(ω) = 2"

!       >> F8 <<
        character (len = *), parameter :: F8_NAME = "Sσ(ω) = RAO(ω)²  Sη(ω)"//ENDL//TAB// &
        "RAO(ω) = 1 / √((1 - (ω/ωn)²)² + (2ξω/ωn)²)"

        character (len = *), parameter :: F8a_NAME = "Sη(ω) = ((4 π³ Hs²) / (ω⁵ Tz⁴)) exp(-(16 π³) / (ω⁴ Tz⁴))"
        character (len = *), parameter :: F8b_NAME = "Sη(ω) = ((4 π³ Hs²) / (ω⁵ Tz⁴)) exp(-(16 π³) / (ω⁴ Tz⁴))"


!       >> F13 <<
        character (len = *), parameter :: F13_NAME = "y'(t) = -2 t y(t)²"//ENDL//"y(0) = 1"
        double precision :: F13_A = 0.0D0
        double precision :: F13_B = 10.0D0
        double precision :: F13_Y0 = 1.0D0

!       >> F14 <<
        character (len = *), parameter :: F14_NAME = "m y''(t) + c y'(t) + k y(t) = F(t)"//ENDL// &
        "m = 1; c = 0.2; k = 1;"//ENDL// &
        "F(t) = 2 sin(w t) + sin(2 w t) + cos(3 w t)"//ENDL// &
        "w = 0.5;"//ENDL// &
        "y'(0) = 0; y(0) = 0;"
        double precision :: F14_M = 1.0D0
        double precision :: F14_C = 0.2D0
        double precision :: F14_K = 1.0D0
        double precision :: F14_W = 0.5D0
        double precision :: F14_Y0 = 0.0D0
        double precision :: F14_DY0 = 0.0D0
        double precision :: F14_A = 0.0D0
        double precision :: F14_B = 100.0D0

!       >> F15 <<
        character (len = *), parameter :: F15_NAME = "z''(t) = -g -k z'(t) |z'(t)|"//ENDL// &
        "z'(0) = 0; z(0) = 0;"//ENDL// &
        "g = 9.806; k = 1;"
        double precision :: F15_G = 9.80600D0
        double precision :: F15_KD = 1.0D0
        double precision :: F15_BY0 = 100.0D0
        double precision :: F15_Y0 = 0.0D0
        double precision :: F15_DY0 = 0.0D0
        double precision :: F15_A = 0.0D0
        double precision :: F15_B = 20.0D0



        double precision :: t1 = 0.0D0
        double precision :: t2 = 0.0D0

        double precision :: wn = 1.00D0
        double precision :: xi = 0.05D0
        double precision :: Hs = 3.0D0
        double precision :: Tz = 5.0D0

        character (len = *), parameter :: F9_NAME = "f(x) = 2 + 2x - x² + 3x³"

        character (len = *), parameter :: F10_NAME = "f(x) = 1 / (1 + x²)"

        character (len = *), parameter :: F11_NAME = "f(x) = exp(- x²/2) / √(2 π)"
        character (len = *), parameter :: F12_NAME = "f(x) = x² exp(- x²/2) / √(2 π)"

!       >> L5-QE <<
        character (len = *), parameter :: FL5_QE1_NAME = 'f(x) = x³ + exp(-x)'
        character (len = *), parameter :: DFL5_QE1_NAME = "f'(x) = 3 x² - exp(-x)"
        character (len = *), parameter :: FL5_QE2_NAME = 'f(x) = ³√x + log(x)'
        character (len = *), parameter :: DFL5_QE2_NAME = "f'(x) = 1 / (3 ³√x²) + (1 / x)"
        character (len = *), parameter :: FL5_QE3_NAME = 'f(x) = 1 - exp(-x² / 25)'
        character (len = *), parameter :: DFL5_QE3_NAME = "f'(x) = (2 x / 25) exp(-x² / 25)"

        
    contains

    function FL5_QE1(x) result (y)
        implicit none
        double precision :: x, y
        y = x ** 3 + DEXP(-x)
        return
    end function

    function DFL5_QE1(x) result (y)
        implicit none
        double precision :: x, y
        y = 3 * x ** 2 - DEXP(-x)
        return
    end function

    function FL5_QE2(x) result (y)
        implicit none
        double precision :: x, y
        y = x ** (1.0D0/3.0D0) + DLOG(x)
        return
    end function

    function DFL5_QE2(x) result (y)
        implicit none
        double precision :: x, y
        y = 1 / (3 * x ** (2.0D0/3.0D0)) + (1 / x)
        return
    end function

    function FL5_QE3(x) result (y)
        implicit none
        double precision :: x, y
        y = 1 - DEXP(-(x ** 2) / 25)
        return
    end function

    function DFL5_QE3(x) result (y)
        implicit none
        double precision :: x, y
        y = (2 * X) / (25 * DEXP((x ** 2) / 25))
        return
    end function


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
    
!   ======= R^n -> R^n functions =======
    function f3(x, n) result (y)
!       R³ -> R³ (n == 3)
        implicit none
        integer :: n
        double precision, dimension(n) :: x, y

        y = (/ &
            (16 * x(1) ** 4 + 16 * x(2) ** 4 + x(3) ** 4) - 16.0D0, &
            x(1) ** 2 + x(2) ** 2 +x(3) ** 2 - 3.0D0, &
            x(1) ** 3 - x(2) + x(3) - 1.0D0 &
            /)
        return
    end function

!   ========== Derivative ===========
    function df3(x, n) result (J)
        implicit none
        integer :: n
        double precision, dimension(n) :: x
        double precision, dimension(n, n) :: J

        J(1, :) = (/ 64 * x(1) ** 3, 64 * x(2) ** 3, 4 * x(3) ** 3 /)
        J(2, :) = (/  2 * x(1)     ,  2 * x(2)     , 2 * x(3)      /)
        J(3, :) = (/  3 * x(1) ** 2,         -1.0D0,         1.0D0 /)
        return
    end function

!   ================== Another function =====================
    function f4(x, n) result (y)
        implicit none
        integer :: n
        double precision, dimension(n) :: x, y
        y = (/ &
            x(1)**2+2*x(2)**2+6*x(3)**2, &
            2*x(2)*(3*x(1)**2+4*x(2)**2+18*x(1)*x(3)+54*x(3)**4), &
            3*(x(1)+20*x(1)**2*x(2)**2+20*x(2)**4+8*x(1)*(x(1)**2+24*x(2)**2)*x(3)+&
            12*(7*x(1)**2+62*x(2)**2)*x(3)**2+432*x(1)*x(3)**3+1116*x(3)**4)&
            /) - (/ 1.0D0, F4_T1, F4_T2 /)
        return
    end function

!   ========== Derivatives ===========    
    function df4(x, n) result (J)
!       R³ -> R³x3 (n == 3)
        implicit none
        integer :: n
        double precision :: x(n), J(n, n)

        J(1, :) = (/ &
            2*x(1), &
            4*x(2), &
            12*x(3) &
            /)
        J(2, :) = (/ &
            12*x(1)*x(2)+36*x(2)*x(3),                     &
            6*x(1)**2+24*x(2)**2+36*x(1)*x(3)+108*x(3)**4, &
            36*x(1)*x(2)+432*x(2)*x(3)**3                  &
            /)
        J(3, :) = (/ &
            3+120*x(1)*x(2)**2+72*x(1)**2*x(3)+576*x(2)**2*x(3)+504*x(1)*x(3)**2+1296*x(3)**3,             &
            120*x(1)**2*x(2)+240*x(2)**3+1152*x(1)*x(2)*x(3)+4464*x(2)*x(3)**2,                            &
            24*x(1)**3+576*x(1)*x(2)**2+504*x(1)**2*x(3)+4464*x(2)**2*x(3)+3888*x(1)*x(3)**2+13392*x(3)**3 &
            /)
        return
    end function

!   ============ One more function =============
    function f5(x, b, m, n) result (z)
        implicit none
        integer :: m, n
        double precision, dimension(m), intent(in) :: b
        double precision, dimension(n), intent(in) :: x
        double precision, dimension(n) :: z

        z = b(1) + (b(2) * (x ** b(3)))
        return
    end function

!   ========= Derivatives ==========
    function df5(x, b, m, n) result (J)
        implicit none
        integer :: m, n
        double precision, dimension(m), intent(in) :: b
        double precision, dimension(n), intent(in) :: x
        double precision, dimension(n, m) :: J
!       m == 3
        J(:, 1) = 1.0D0
        J(:, 2) = x ** b(3)
        J(:, 3) = b(2) * DLOG(x) * (x ** b(3))
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
        z = (4 * (Hs**2) * (PI**3)) / ( DEXP( (16 * (PI**3))/((Tz*w)**4) ) * (Tz**4) * (w**5) )
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
        z = Ss(w, Sn2)
        return
    end function

    function f8b(w) result (z)
        implicit none
        double precision :: w, z
        z = (w ** 2) * Ss(w, Sn2)
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
        if (y >= 0) then
            u = - F15_G
        else
            u = - F15_G - F15_KD * dy * DABS(dy)
        end if
        return
    end function

    end module Func