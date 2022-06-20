module condition
    use equation
    use integ
    implicit none
    private

    abstract interface
        real(8) function cnd(nrg)
        real(8), intent(in) :: nrg
        end function cnd
    end interface

! Solution search method
    procedure(cnd), pointer :: cndPtr => cndLog

! Solution search prescision (affects stationary levels)
    integer, parameter :: cdots = 2000

    public :: cndWron, cndLog, getCND, cndPtr
contains

    real(8) function cndWron(nrg)
        implicit none
        real(8), intent(in) :: nrg

        abstract interface
            real(8) function fct(x, y, z, nrg)
            real(8), intent(in) :: x, y, z, nrg
            end function fct
        end interface
        procedure(fct), pointer :: fptr, gptr

        real(8) :: yl, zl, yr, zr, stp, ystore, zstore, x
        integer :: i

        fptr => func
        gptr => gunc

        yl = yleft
        zl = zleft
        stp = (xcrs - xnou) / cdots
        do i = 0, cdots-1
            x = xnou + i*stp
            call integL(fptr, gptr, nrg, x, stp, yl, zl, ystore, zstore)
            yl = ystore
            zl = zstore
        end do

        yr = yright(nrg)
        zr = zright(nrg)
        stp = (xinf - xcrs) / cdots
        do i = 0, cdots-1
            x = xinf - i*stp
            call integR(fptr, gptr, nrg, x, stp, yr, zr, ystore, zstore)
            yr = ystore
            zr = zstore
        end do

        cndWron = yl*zr - zl*yr

    end function cndWron

    real(8) function cndLog(nrg)
        implicit none
        real(8), intent(in) :: nrg

        abstract interface
            real(8) function fct(x, y, z, nrg)
            real(8), intent(in) :: x, y, z, nrg
            end function fct
        end interface
        procedure(fct), pointer :: fptr, gptr

        real(8) :: yl, zl, yr, zr, stp, ystore, zstore, x
        integer :: i

        fptr => func
        gptr => gunc

        yl = yleft
        zl = zleft
        stp = (xcrs - xnou) / cdots
        do i = 0, cdots-1
            x = xnou + i*stp
            call integL(fptr, gptr, nrg, x, stp, yl, zl, ystore, zstore)
            yl = ystore
            zl = zstore
        end do

        yr = yright(nrg)
        zr = zright(nrg)
        stp = (xinf - xcrs) / cdots
        do i = 0, cdots-1
            x = xinf - i*stp
            call integR(fptr, gptr, nrg, x, stp, yr, zr, ystore, zstore)
            yr = ystore
            zr = zstore
        end do

        cndLog = zr/yr - zl/yl

    end function cndLog

    subroutine getCND()
        implicit none
        real(8) :: erg, stp
        integer :: i

        stp = (einf-enou)/cdots
        open(unit = 1, file = "graphics/cndPlot.dat")
        do i = 1, cdots+1
            erg = -einf + (i-1)*stp
            write(1, *) erg, cndPtr(erg)
        end do
        close(1)

        print *, "Cnd(erg) is drawn"
    end subroutine getCND

end module condition
