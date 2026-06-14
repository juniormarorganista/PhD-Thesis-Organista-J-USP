program case12_fields_reference
    use iso_fortran_env, only: real64
    implicit none

    integer :: nx, ny, i, j
    character(len=512) :: outfile, arg
    real(real64) :: x, y, t, dx, dy

    nx = 3
    ny = 3
    t = 0.0_real64
    outfile = "case12_fields_reference.csv"

    call get_command_argument(1, outfile)

    call get_command_argument(2, arg)
    if (len_trim(arg) > 0) read(arg, *) nx

    call get_command_argument(3, arg)
    if (len_trim(arg) > 0) read(arg, *) ny

    call get_command_argument(4, arg)
    if (len_trim(arg) > 0) read(arg, *) t

    if (nx < 1 .or. ny < 1) then
        error stop "nx and ny must be positive."
    end if

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
    write(10, '(A)') "x,y,u,v,wz,psi,txx,txy,tyy"

    do j = 1, ny
        y = real(j - 1, real64) * dy
        do i = 1, nx
            x = real(i - 1, real64) * dx
            write(10, '(ES26.17,",",ES26.17,",",ES26.17,",",ES26.17,",",ES26.17,",",ES26.17,",",ES26.17,",",ES26.17,",",ES26.17)') &
                x, y, &
                u_case12_default(x, y, t), &
                v_case12_default(x, y, t), &
                wz_case12_default(x, y, t), &
                psi_case12_default(x, y, t), &
                txx_case12_default(x, y, t), &
                txy_case12_default(x, y, t), &
                tyy_case12_default(x, y, t)
        end do
    end do

    close(10)

contains

    pure function u_case12_default(x, y, t) result(value)
        real(real64), intent(in) :: x, y, t
        real(real64) :: value

        real(real64), parameter :: pii = acos(-1.0_real64)
        real(real64), parameter :: at = 0.0_real64

        value = ((-1.0_real64 + x) * x * sin(pii * x) &
            * (pii * y * cos(pii * y) + sin(pii * y))) / exp(at * t)
    end function u_case12_default

    pure function v_case12_default(x, y, t) result(value)
        real(real64), intent(in) :: x, y, t
        real(real64) :: value

        real(real64), parameter :: pii = acos(-1.0_real64)
        real(real64), parameter :: at = 0.0_real64

        value = -((y * (pii * (-1.0_real64 + x) * x * cos(pii * x) &
            + (-1.0_real64 + 2.0_real64 * x) * sin(pii * x)) * sin(pii * y)) / exp(at * t))
    end function v_case12_default

    pure function wz_case12_default(x, y, t) result(value)
        real(real64), intent(in) :: x, y, t
        real(real64) :: value

        real(real64), parameter :: pii = acos(-1.0_real64)
        real(real64), parameter :: at = 0.0_real64

        value = (2.0_real64 * (pii * (-1.0_real64 + x) * x * cos(pii * y) * sin(pii * x) &
            + y * (pii * (-1.0_real64 + 2.0_real64 * x) * cos(pii * x) &
            + (1.0_real64 - pii**2 * (-1.0_real64 + x) * x) * sin(pii * x)) * sin(pii * y))) &
            / exp(at * t)
    end function wz_case12_default

    pure function psi_case12_default(x, y, t) result(value)
        real(real64), intent(in) :: x, y, t
        real(real64) :: value

        real(real64), parameter :: pii = acos(-1.0_real64)
        real(real64), parameter :: at = 0.0_real64

        value = ((-1.0_real64 + x) * x * y * sin(pii * x) * sin(pii * y)) / exp(at * t)
    end function psi_case12_default

    pure function txx_case12_default(x, y, t) result(value)
        real(real64), intent(in) :: x, y, t
        real(real64) :: value

        real(real64), parameter :: pii = acos(-1.0_real64)
        real(real64), parameter :: at = 0.0_real64
        real(real64), parameter :: betann = 0.5_real64

        real(real64), parameter :: a1xx = 1.0_real64
        real(real64), parameter :: a2xx = 1.0_real64
        real(real64), parameter :: a3xx = 1.0_real64
        real(real64), parameter :: b1xx = 1.0_real64
        real(real64), parameter :: b2xx = 1.0_real64
        real(real64), parameter :: b3xx = 1.0_real64

        value = ((1.0_real64 - betann) * (a1xx + a3xx * cos(pii * x) + a2xx * sin(pii * x)) &
            * (b1xx + b3xx * cos(pii * y) + b2xx * sin(pii * y))) / exp(at * t)
    end function txx_case12_default

    pure function txy_case12_default(x, y, t) result(value)
        real(real64), intent(in) :: x, y, t
        real(real64) :: value

        real(real64), parameter :: pii = acos(-1.0_real64)
        real(real64), parameter :: at = 0.0_real64
        real(real64), parameter :: betann = 0.5_real64

        real(real64), parameter :: a1xy = 1.0_real64
        real(real64), parameter :: a2xy = 1.0_real64
        real(real64), parameter :: a3xy = 1.0_real64
        real(real64), parameter :: b1xy = 1.0_real64
        real(real64), parameter :: b2xy = 1.0_real64
        real(real64), parameter :: b3xy = 1.0_real64

        value = ((1.0_real64 - betann) * (a1xy + a3xy * cos(pii * x) + a2xy * sin(pii * x)) &
            * (b1xy + b3xy * cos(pii * y) + b2xy * sin(pii * y))) / exp(at * t)
    end function txy_case12_default

    pure function tyy_case12_default(x, y, t) result(value)
        real(real64), intent(in) :: x, y, t
        real(real64) :: value

        real(real64), parameter :: pii = acos(-1.0_real64)
        real(real64), parameter :: at = 0.0_real64
        real(real64), parameter :: betann = 0.5_real64

        real(real64), parameter :: a1yy = 1.0_real64
        real(real64), parameter :: a2yy = 1.0_real64
        real(real64), parameter :: a3yy = 1.0_real64
        real(real64), parameter :: b1yy = 1.0_real64
        real(real64), parameter :: b2yy = 1.0_real64
        real(real64), parameter :: b3yy = 1.0_real64

        value = ((1.0_real64 - betann) * (a1yy + a3yy * cos(pii * x) + a2yy * sin(pii * x)) &
            * (b1yy + b3yy * cos(pii * y) + b2yy * sin(pii * y))) / exp(at * t)
    end function tyy_case12_default

end program case12_fields_reference
