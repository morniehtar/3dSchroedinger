module drawWF
    use equation
    use integ
    implicit none
    private

! WF smoothness
    integer, parameter :: dots = 10000 ! 10000 works

    public :: getWF, dots
contains

    function getWF(nrg)
        implicit none
        real(8), intent(in) :: nrg
        real(8), dimension(0:2*dots,2) :: getWF

        abstract interface
            real(8) function fct(x, y, z, nrg)
            real(8), intent(in) :: x, y, z, nrg
            end function fct
        end interface

        procedure(fct), pointer :: fptr, gptr
        real(8), dimension(4) :: k, l
        real(8) :: stp, zst, ycoeff, ystore, zstore, C
        integer :: i

        fptr => func
        gptr => gunc

    ! Calculate left-side WF
        getWF(0,2) = yleft
        zst = zleft

        stp = abs(xcrs - xnou) / dots
        do i = 0, dots
            getWF(i,1) = xnou + i*stp
            call integL(fptr, gptr, nrg, getWF(i,1), stp, getWF(i,2), zst, ystore, zstore)
            if (i.ne.dots) getWF(i+1,2) = ystore
            zst = zstore
        end do

        ycoeff = getWF(dots,2)

    ! Calculate right-side WF
        getWF(2*dots,2) = yright(nrg)
        zst = zright(nrg)

        stp = abs(xinf - xcrs) / dots
        do i = 2*dots, dots, -1
            getWF(i,1) = xinf - (2*dots-i)*stp
            call integR(fptr, gptr, nrg, getWF(i,1), stp, getWF(i,2), zst, ystore, zstore)
            if (i.ne.dots) getWF(i-1,2) = ystore
            zst = zstore
        end do


    ! Raw data continuity fixing
        C = getWF(dots,2)/ycoeff
        do i = 0, dots-1
            getWF(i,2) = getWF(i,2)*C
        end do

    ! Raw data normalization fixing
        C = discrItgr(wfSq(getWF))
        do i = 0, 2*dots
            getWF(i,2) = getWF(i,2)/dsqrt(C)
        end do

        print*, "WF is retrieved"

    end function getWF

    pure function wfSq(farr)
        real(8), dimension(0:2*dots,2), intent(in) :: farr
        real(8), dimension(0:2*dots,2) :: wfSq
        integer :: i
        do i=0, 2*dots
            wfSq(i, 1) = farr(i, 1)
            wfSq(i, 2) = farr(i, 2)**2
        end do
    end function wfSq

    real(8) function discrItgr(farr)
        implicit none
        real(8), dimension(0:2*dots,2), intent(in) :: farr

        real(8) :: sres, er, next, stp
        integer :: i

        discrItgr = 0
        er = 0
        if (mod(dots, 3).eq.0) then

            stp = farr(1,1) - farr(0,1)
            do i = 0, dots-3, 3
                next = (3*stp/8)*(farr(i,2)+3*farr(i+1,2)+3*farr(i+2,2)+farr(i+3,2)) - er
                sres = discrItgr + next
                er = (sres - discrItgr) - next
                discrItgr = sres
            end do

            stp = farr(dots+1,1) - farr(dots,1)
            do i = dots, 2*dots-3, 3
                next = (3*stp/8)*(farr(i,2)+3*farr(i+1,2)+3*farr(i+2,2)+farr(i+3,2)) - er
                sres = discrItgr + next
                er = (sres - discrItgr) - next
                discrItgr = sres
            end do

        elseif (mod(dots, 3).eq.1) then

            stp = farr(1,1) - farr(0,1)
            do i = 0, dots-7, 3
                next = (3*stp/8)*(farr(i,2)+3*farr(i+1,2)+3*farr(i+2,2)+farr(i+3,2)) - er
                sres = discrItgr + next
                er = (sres - discrItgr) - next
                discrItgr = sres
            end do
            do i = dots-4, dots-2, 2
                next = (stp/3)*(farr(i,2)+4*farr(i+1,2)+farr(i+2,2)) - er
                sres = discrItgr + next
                er = (sres - discrItgr) - next
                discrItgr = sres
            end do

            stp = farr(dots+1,1) - farr(dots,1)
            do i = dots, 2*dots-7, 3
                next = (3*stp/8)*(farr(i,2)+3*farr(i+1,2)+3*farr(i+2,2)+farr(i+3,2)) - er
                sres = discrItgr + next
                er = (sres - discrItgr) - next
                discrItgr = sres
            end do
            do i = 2*dots-4, 2*dots-2, 2
                next = (stp/3)*(farr(i,2)+4*farr(i+1,2)+farr(i+2,2)) - er
                sres = discrItgr + next
                er = (sres - discrItgr) - next
                discrItgr = sres
            end do

        elseif (mod(dots, 3).eq.2) then

            stp = farr(1,1) - farr(0,1)
            do i = 0, dots-5, 3
                next = (3*stp/8)*(farr(i,2)+3*farr(i+1,2)+3*farr(i+2,2)+farr(i+3,2)) - er
                sres = discrItgr + next
                er = (sres - discrItgr) - next
                discrItgr = sres
            end do
            next = (stp/3)*(farr(dots-2,2)+4*farr(dots-1,2)+farr(dots,2)) - er
            sres = discrItgr + next
            er = (sres - discrItgr) - next
            discrItgr = sres

            stp = farr(dots+1,1) - farr(dots,1)
            do i = dots, 2*dots-5, 3
                next = (3*stp/8)*(farr(i,2)+3*farr(i+1,2)+3*farr(i+2,2)+farr(i+3,2)) - er
                sres = discrItgr + next
                er = (sres - discrItgr) - next
                discrItgr = sres
            end do
            next = (stp/3)*(farr(2*dots-2,2)+4*farr(2*dots-1,2)+farr(2*dots,2)) - er
            sres = discrItgr + next
            er = (sres - discrItgr) - next
            discrItgr = sres

        end if
    end function discrItgr

end module drawWF
