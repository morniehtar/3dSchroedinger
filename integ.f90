module integ
    implicit none
    private

    abstract interface
        real(8) function fct(x, y, z, nrg)
        real(8), intent(in) :: x, y, z, nrg
        end function fct
    end interface

    public :: integL, integR
contains
    subroutine integL(fptr, gptr, nrg, x, stp, ySt, zSt, yEnd, zEnd)
        implicit none

        real(8), intent(in) :: ySt, zSt, nrg, x, stp
        procedure(fct), pointer, intent(in) :: fptr, gptr
        real(8), intent(out) :: yEnd, zEnd

        real(8), dimension(4) :: k, l
        integer :: i

        k(1) = fptr(x, ySt, zSt, nrg)
        l(1) = gptr(x, ySt, zSt, nrg)

        k(2) = fptr(x+stp/2, ySt+k(1)*stp/2, zSt+k(1)*stp/2, nrg)
        l(2) = gptr(x+stp/2, ySt+l(1)*stp/2, zSt+l(1)*stp/2, nrg)

        k(3) = fptr(x+stp/2, ySt+k(2)*stp/2, zSt+k(2)*stp/2, nrg)
        l(3) = gptr(x+stp/2, ySt+l(2)*stp/2, zSt+l(2)*stp/2, nrg)

        k(4) = fptr(x+stp, ySt+k(3)*stp, zSt+k(3)*stp, nrg)
        l(4) = gptr(x+stp, ySt+l(3)*stp, zSt+l(3)*stp, nrg)

        yEnd = ySt + stp*(k(1)+2*k(2)+2*k(3)+k(4))/6
        zEnd = zSt + stp*(l(1)+2*l(2)+2*l(3)+l(4))/6

    end subroutine integL

    subroutine integR(fptr, gptr, nrg, x, stp, ySt, zSt, yEnd, zEnd)
        implicit none

        real(8), intent(in) :: ySt, zSt, nrg, x, stp
        procedure(fct), pointer, intent(in) :: fptr, gptr
        real(8), intent(out) :: yEnd, zEnd

        real(8), dimension(4) :: k, l
        integer :: i

        k(1) = fptr(x, ySt, zSt, nrg)
        l(1) = gptr(x, ySt, zSt, nrg)

        k(2) = fptr(x-stp/2, ySt-k(1)*stp/2, zSt-k(1)*stp/2, nrg)
        l(2) = gptr(x-stp/2, ySt-l(1)*stp/2, zSt-l(1)*stp/2, nrg)

        k(3) = fptr(x-stp/2, ySt-k(2)*stp/2, zSt-k(2)*stp/2, nrg)
        l(3) = gptr(x-stp/2, ySt-l(2)*stp/2, zSt-l(2)*stp/2, nrg)

        k(4) = fptr(x-stp, ySt-k(3)*stp, zSt-k(3)*stp, nrg)
        l(4) = gptr(x-stp, ySt-l(3)*stp, zSt-l(3)*stp, nrg)

        yEnd = ySt - stp*(k(1)+2*k(2)+2*k(3)+k(4))/6
        zEnd = zSt - stp*(l(1)+2*l(2)+2*l(3)+l(4))/6

    end subroutine integR

end module integ
