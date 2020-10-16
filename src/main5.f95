program main5
    use Func
    use Calc
    use Util
    implicit none

    integer :: i

    integer :: n = 10
    double precision :: a, b, r, s

    character(len=32) :: n_chr

!   Command-line Args
    integer :: argc

!    ENABLE_DEBUG = .TRUE.

!   Random seed definition
    call init_random_seed()

!   Get Command-Line Args
    argc = iargc()

    if (argc == 0) then
        goto 099
    else if (argc == 1) then
        call getarg(1, n_chr)
        read (unit=n_chr, fmt=*) n
        goto 099
    else
        goto 101
    end if

!   ====== Begin =====================================    
099 call info("Nº de pontos de integração: "//STR(n))
    goto 200

!   ====== Success ===================================
100 call info(':: Sucesso ::')
    goto 1
!   ====== Errors ====================================
101 call error('Este programa não aceita mais do que um parâmetro (número de pontos de integração).')
    goto 1
!   ====== Finish ====================================
1   stop
!   ==================================================

!   ===============================
200 call info(ENDL//"2) "//F6_NAME)
    i = 0
    goto 201

201 i = i + 1
    go to(210, 220, 300), i
!   ===============================
300 call info(ENDL//"3) "//F78_NAME)
    i = 0
    goto 301

301 i = i + 1
    go to(310, 400), i
!   ===============================
400 call info(ENDL//"4) "//F78_NAME)
    i = 0
    goto 401

401 i = i + 1
    go to(410, 500), i
!   ===============================   
500 call info(ENDL//"5) "//F9_NAME)
    i = 0
    goto 501

501 i = i + 1
    go to(510, 600), i
!   ===============================  
600 call info(ENDL//"6) "//F10_NAME)
    i = 0
    goto 601

601 i = i + 1
    go to(610, 700), i
!   ===============================
700 call info(ENDL//"7)")
    i = 0
    goto 701

701 i = i + 1
    go to(710, 720, 800), i
!   =============================== 
800 goto 100

210 a = 0.0D0
    b = 1.0D0
    call info("[a, b] = ["//DSTR(a)//", "//DSTR(b)//"]")
    call info(":: Integração Polinomial ::")
    s = num_int(f6, a, b, n, kind="polynomial")
    call show("∫f(x) dx", s)
    call info(":: Quadratura de Gauss-Legendre ::")
    s = num_int(f6, a, b, n, kind="gauss-legendre")
    call show("∫f(x) dx", s)
    goto 201

220 a = 0.0D0
    b = 5.0D0
    call info("[a, b] = ["//DSTR(a)//", "//DSTR(b)//"]")
    call info(":: Integração Polinomial ::")
    s = num_int(f6, a, b, n, kind="polynomial")
    call show("∫f(x) dx", s)
    call info(":: Quadratura de Gauss-Legendre ::")
    s = num_int(f6, a, b, n, kind="gauss-legendre")
    call show("∫f(x) dx", s)
    goto 201

310 call info("m0) "//ENDL//TAB//F7a_NAME)
    a = 0.00D0
    b = 10.0D0
    call info("[a, b] = ["//DSTR(a)//", "//DSTR(b)//"]")
    call info(":: Integração Polinomial ::")
    s = num_int(f7a, a, b, n, kind="polynomial")
    call show("∫f(ω) dω", s)
    call info(":: Quadratura de Gauss-Legendre ::")
    s = num_int(f7a, a, b, n, kind="gauss-legendre")
    call show("∫f(ω) dω", s)

    call info("m2) "//ENDL//TAB//F7b_NAME)
    call info(":: Integração Polinomial ::")
    s = num_int(f7b, a, b, n, kind="polynomial")
    call show("∫f(ω) dω", s)
    call info(":: Quadratura de Gauss-Legendre ::")
    s = num_int(f7b, a, b, n, kind="gauss-legendre")
    call show("∫f(ω) dω", s)
    goto 301

410 call info("m0) "//ENDL//TAB//F8a_NAME)
    a = 0.00D0
    b = 10.0D0
    call info("[a, b] = ["//DSTR(a)//", "//DSTR(b)//"]")
    call info(":: Integração Polinomial ::")
    s = num_int(f8a, a, b, n, kind="polynomial")
    call show("∫f(ω) dω", s)
    call info(":: Quadratura de Gauss-Legendre ::")
    s = num_int(f8a, a, b, n, kind="gauss-legendre")
    call show("∫f(ω) dω", s)

    call info("m2) "//ENDL//TAB//F8b_NAME)
    call info(":: Integração Polinomial ::")
    s = num_int(f8b, a, b, n, kind="polynomial")
    call show("∫f(ω) dω", s)
    call info(":: Quadratura de Gauss-Legendre ::")
    s = num_int(f8b, a, b, n, kind="gauss-legendre")
    call show("∫f(ω) dω", s)
    goto 401

510 a = 0.0D0
    b = 4.0D0
    call info("[a, b] = ["//DSTR(a)//", "//DSTR(b)//"]")
    call info(":: Integração Polinomial ::")
    call show('n', 4.0D0)
    s = num_int(f9, a, b, 4, kind="polynomial")
    call show("∫f(x) dx", s)
    call info(":: Quadratura de Gauss-Legendre ::")
    call show('n', 2.0D0)
    s = num_int(f9, a, b, 2, kind="gauss-legendre")
    call show("∫f(x) dx", s)
    goto 501

610 a = 0.0D0
    b = 3.0D0
    call info("[a, b] = ["//DSTR(a)//", "//DSTR(b)//"]")
    call info(":: Integração Polinomial ::")
    s = num_int(f10, a, b, n, kind="polynomial")
    call show("∫f(x) dx", s)
    call info(":: Quadratura de Gauss-Legendre ::")
    s = num_int(f10, a, b, n, kind="gauss-legendre")
    call show("∫f(x) dx", s)
    goto 601

710 call info("A1) "//F11_NAME)
    a = DNINF
    b = 1.0D0
    call info("[a, b] = [-∞, "//DSTR(b)//"]")
    call info(":: Quadratura de Gauss-Hermite e de Gauss-Legendre ::")
    r = num_int(f11a, a, -a, n, kind="gauss-hermite")
    s = num_int(f11b, -b, b, n, kind="gauss-legendre")
    call show("∫f(x) dx", r + s)
    goto 701

720 call info("A2) "//F12_NAME)
    a = DNINF
    b = DINF
    call info("[a, b] = [-∞, ∞]")
    call info(":: Quadratura de Gauss-Hermite ::")
    s = num_int(f12, a, b, n, kind="gauss-hermite")
    call show("∫f(x) dx", s)
    goto 701
end program main5