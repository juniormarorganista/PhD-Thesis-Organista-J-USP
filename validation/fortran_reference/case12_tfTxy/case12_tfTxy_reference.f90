program case12_tftxy_reference
    use iso_fortran_env, only: real64
    implicit none

    integer :: nx, ny, i, j
    character(len=512) :: outfile, arg
    real(real64) :: x, y, t, dx, dy, value

    nx = 3
    ny = 3
    t = 0.0_real64
    outfile = "case12_tfTxy_reference.csv"

    call get_command_argument(1, outfile)

    call get_command_argument(2, arg)
    if (len_trim(arg) > 0) read(arg, *) nx

    call get_command_argument(3, arg)
    if (len_trim(arg) > 0) read(arg, *) ny

    call get_command_argument(4, arg)
    if (len_trim(arg) > 0) read(arg, *) t

    if (nx < 1 .or. ny < 1) error stop "nx and ny must be positive."

    if (nx > 1) then
        dx = 1.0_real64 / real(nx - 1, real64)
    else
        dx = 0.0_real64
    end if

    if (ny > 1) then
        dy = 1.0_real64 / real(ny - 1, real64)
    else
        dy = 0.0_real64
    end if

    open(unit=10, file=trim(outfile), status="replace", action="write")
    write(10, '(A)') "x,y,tfTxy"

    do j = 1, ny
        y = real(j - 1, real64) * dy
        do i = 1, nx
            x = real(i - 1, real64) * dx
            value = tftxy_case12_default(x, y, t)
            write(10, '(ES26.17,",",ES26.17,",",ES26.17)') x, y, value
        end do
    end do

    close(10)

contains

    pure function tftxy_case12_default(x, y, t) result(value)
        real(real64), intent(in) :: x, y, t
        real(real64) :: value

        real(real64), parameter :: pii = acos(-1.0_real64)
        real(real64), parameter :: Rey = 1.0_real64
        real(real64), parameter :: Wi = 1.0_real64
        real(real64), parameter :: betann = 0.5_real64
        real(real64), parameter :: epsylon = 0.0_real64
        real(real64), parameter :: xi = 0.0_real64
        real(real64), parameter :: alphaG = 0.0_real64
        real(real64), parameter :: at = 0.0_real64

        real(real64), parameter :: a1xx = 1.0_real64, a2xx = 1.0_real64, a3xx = 1.0_real64
        real(real64), parameter :: b1xx = 1.0_real64, b2xx = 1.0_real64, b3xx = 1.0_real64
        real(real64), parameter :: a1xy = 1.0_real64, a2xy = 1.0_real64, a3xy = 1.0_real64
        real(real64), parameter :: b1xy = 1.0_real64, b2xy = 1.0_real64, b3xy = 1.0_real64
        real(real64), parameter :: a1yy = 1.0_real64, a2yy = 1.0_real64, a3yy = 1.0_real64
        real(real64), parameter :: b1yy = 1.0_real64, b2yy = 1.0_real64, b3yy = 1.0_real64

        real(real64) :: exp_at, exp_2at
        real(real64) :: cos_px, sin_px, cos_py, sin_py
        real(real64) :: x_poly
        real(real64) :: Axx, Bxx, Axy, Bxy, Ayy, Byy
        real(real64) :: dpsi_dy_like, dudx_shape
        real(real64) :: mixed_xy_shape, uxx_shape
        real(real64) :: x_derivative_xy, y_derivative_xy
        real(real64) :: shear_shape, tensor_sum
        real(real64) :: first, giesekus, relaxation, wi_inner

        exp_at = exp(at * t)
        exp_2at = exp(2.0_real64 * at * t)

        cos_px = cos(pii * x)
        sin_px = sin(pii * x)
        cos_py = cos(pii * y)
        sin_py = sin(pii * y)
        x_poly = (-1.0_real64 + x) * x

        Axx = a1xx + a3xx * cos_px + a2xx * sin_px
        Bxx = b1xx + b3xx * cos_py + b2xx * sin_py
        Axy = a1xy + a3xy * cos_px + a2xy * sin_px
        Bxy = b1xy + b3xy * cos_py + b2xy * sin_py
        Ayy = a1yy + a3yy * cos_px + a2yy * sin_px
        Byy = b1yy + b3yy * cos_py + b2yy * sin_py

        dpsi_dy_like = pii * y * cos_py + sin_py
        dudx_shape = pii * x_poly * cos_px + (-1.0_real64 + 2.0_real64 * x) * sin_px
        mixed_xy_shape = -2.0_real64 * cos_py + pii * y * sin_py
        uxx_shape = 2.0_real64 * pii * (1.0_real64 - 2.0_real64 * x) * cos_px &
            + (-2.0_real64 + pii**2 * x_poly) * sin_px

        x_derivative_xy = -(a2xy * cos_px) + a3xy * sin_px
        y_derivative_xy = b2xy * cos_py - b3xy * sin_py

        shear_shape = pii * x_poly * cos_py * sin_px &
            - y * (pii * (-1.0_real64 + 2.0_real64 * x) * cos_px + sin_px) * sin_py

        tensor_sum = &
            a1xx * b1xx &
            + a1yy * b1yy &
            + cos_py * (a1xx * b3xx + a1yy * b3yy + (a2xx * b3xx + a2yy * b3yy) * sin_px) &
            + (a1xx * b2xx + a1yy * b2yy) * sin_py &
            + sin_px * (a2xx * b1xx + a2yy * b1yy + (a2xx * b2xx + a2yy * b2yy) * sin_py) &
            + cos_px * (a3xx * b1xx + a3yy * b1yy &
                + (a3xx * b3xx + a3yy * b3yy) * cos_py &
                + (a3xx * b2xx + a3yy * b2yy) * sin_py)

        first = 2.0_real64 * (-1.0_real64 + betann) * exp_at * shear_shape / Rey

        giesekus = alphaG * (1.0_real64 - betann) * Rey * Wi * Axy * Bxy * (Axx * Bxx + Ayy * Byy)

        relaxation = (1.0_real64 - betann) * exp_at * Axy * Bxy &
            * (1.0_real64 + (epsylon * Rey * Wi * Axx * Bxx) / exp_at &
            + (epsylon * Rey * Wi * Ayy * Byy) / exp_at)

        wi_inner = &
            y * Axx * uxx_shape * sin_py * Bxx &
            + at * exp_at * Axy * Bxy &
            + pii * y * cos_py * Axy * dudx_shape * Bxy &
            + Axy * dudx_shape * sin_py * Bxy &
            - pii * x_poly * cos_px * Axy * dpsi_dy_like * Bxy &
            - (-1.0_real64 + x) * sin_px * Axy * dpsi_dy_like * Bxy &
            - x * sin_px * Axy * dpsi_dy_like * Bxy &
            + pii * x_poly * sin_px * x_derivative_xy * dpsi_dy_like * Bxy &
            + pii * y * Axy * dudx_shape * sin_py * y_derivative_xy &
            + pii * x_poly * sin_px * Ayy * Byy * (2.0_real64 * cos_py - pii * y * sin_py) &
            + xi * (-(pii * x_poly * cos_py * sin_px) &
                + y * (pii * (-1.0_real64 + 2.0_real64 * x) * cos_px + sin_px) * sin_py) * tensor_sum

        value = (first + giesekus + relaxation + (-1.0_real64 + betann) * Wi * wi_inner) / exp_2at
    end function tftxy_case12_default

end program case12_tftxy_reference
