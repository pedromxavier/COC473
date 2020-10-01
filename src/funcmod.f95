!   Func Module

    module Func
        implicit none

        double precision :: g = 9.80600D0
        double precision :: k = 0.00341D0

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
    
    end module Func
