module config
    implicit none
    public

! Integration precision (also wave function smoothness)
    integer, parameter :: prec = 10000 ! 10000 works
! Cnd(erg) smoothness
    integer, parameter :: dots = 1000
! Threads count
    !integer, parameter :: thr = 4
!---------------------------------------------------------------
! Effective zero
    real(8), parameter :: xnou = 7d-3
! Effective infinity
    real(8), parameter :: xinf = 7d0
! Cross-linking point
    real(8), parameter :: xcrs = 0.434d0
!---------------------------------------------------------------
! Energy effective nought
    real(8), parameter :: enou = 4d-3 ! >0; 4d-3 for U1; 5d-2 generally
! Energy effective infinity
    real(8), parameter :: einf = 10d0 ! >0
! Energy search precision
    integer, parameter :: eprec = 600
!---------------------------------------------------------------
! Machine epsilon (root search limiter)
    real(8), parameter :: eps = 1d-9 ! 1d-16 is machine epsilon
!---------------------------------------------------------------
end module config

module potent
    implicit none
    private
    procedure(Yl), pointer, public :: uptr => Yl

    abstract interface
        pure real(8) function ufct(x)
        real(8), intent(in) :: x
        end function ufct
    end interface

contains

    pure real(8) function Yl(r)
        implicit none
        real(8), intent(in) :: r
        real(8), parameter :: g = 2d0
        integer, parameter :: l = 0
        Yl = -g*exp(-r)/r + 0.5d0*l*(l+1)/r**2
    end function Yl

end module potent

program hello
    use config
    use potent
    use search
    implicit none


    real(8) :: bt, tp, estep
    real(8), dimension(:), allocatable :: stErg
    real(8), dimension(0:2*prec, 2) :: WF

    integer :: i, j

    abstract interface
        real(8) function cnd(nrg)
        real(8), intent(in) :: nrg
        end function cnd
    end interface

    procedure(cnd), pointer :: cndPtr
    character(len=13) :: name = "k??wfData.dat"

    !vvvvvvvvvvvvvvv
    cndPtr => cndWron
    !^^^^^^^^^^^^^^^

! Drawing cnd(erg)
    call drawCND()

! Searching stationary levels
    allocate(stErg(0))
    estep = (einf-enou)/eprec
    do i = 1, eprec
        bt = -einf + (i-1)*estep
        tp = -einf + i*estep
        call getRoot(cndPtr, bt, tp, stErg)
    end do

! WF output
    do j = 1, size(stErg)
        WF=getWF(stErg(j))
        write(unit=name(2:3), fmt="(i2.2)") j
        open(unit = 1, file = name)
        do i = 0, 2*prec
            write(1,*) WF(i, 1), WF(i, 2)
        end do
        close(1)
    end do

    deallocate(stErg)
contains

    real(8) function cndWron(nrg)
        implicit none

        real(8), intent(in) :: nrg
        real(8) :: yleft, zleft, yright, zright

        call integL(nrg, yleft, zleft)
        call integR(nrg, yright, zright)

        cndWron = yleft*zright - zleft*yright

    end function cndWron

    !y'=func(x, y, z); y = psi
    !z'=gunc(x, y, z); zdx = dpsi

    pure real(8) function func(x, y, z, nrg)
        implicit none
        real(8), intent(in) :: x, y, z, nrg
        func = z
    end function func

    pure real(8) function gunc(x, y, z, nrg)
        implicit none
        real(8), intent(in) :: x, y, z, nrg
        !ħ=m=1
        gunc = 2.d0*(uptr(x)-nrg)*y
    end function gunc

    subroutine integL(nrg, yleft, zleft)
    ! tested
        use config
        implicit none

        real(8), intent(inout) :: yleft, zleft
        real(8), intent(in) :: nrg

        abstract interface
            real(8) function fct(x, y, z, nrg)
            real(8), intent(in) :: x, y, z, nrg
            end function fct
        end interface

        procedure(fct), pointer :: fptr, gptr
        real(8), dimension(4) :: k, l
        real(8) :: x, step
        integer :: i

        fptr => func
        gptr => gunc

        yleft = 0
        zleft = 1

        step = (xcrs - xnou) / prec

        do i = 0, prec
            x = xnou + i*step

            k(1) = fptr(x, yleft, zleft, nrg)
            l(1) = gptr(x, yleft, zleft, nrg)

            k(2) = fptr(x+step/2, yleft+k(1)*step/2, zleft+k(1)*step/2, nrg)
            l(2) = gptr(x+step/2, yleft+l(1)*step/2, zleft+l(1)*step/2, nrg)

            k(3) = fptr(x+step/2, yleft+k(2)*step/2, zleft+k(2)*step/2, nrg)
            l(3) = gptr(x+step/2, yleft+l(2)*step/2, zleft+l(2)*step/2, nrg)

            k(4) = fptr(x+step, yleft+k(3)*step, zleft+k(3)*step, nrg)
            l(4) = gptr(x+step, yleft+l(3)*step, zleft+l(3)*step, nrg)

            yleft = yleft + step*(k(1)+2*k(2)+2*k(3)+k(4))/6
            zleft = zleft + step*(l(1)+2*l(2)+2*l(3)+l(4))/6

        end do

    end subroutine integL

    subroutine integR(nrg, yright, zright)
    ! tested
        use config
        implicit none

        real(8), intent(inout) :: yright, zright
        real(8), intent(in) :: nrg

        abstract interface
            real(8) function fct(x, y, z, nrg)
            real(8), intent(in) :: x, y, z, nrg
            end function fct
        end interface

        procedure(fct), pointer :: fptr, gptr
        real(8), dimension(4) :: k, l
        real(8) :: x, step
        integer :: i

        fptr => func
        gptr => gunc

        yright = exp(-sqrt(-2d0*nrg)*xinf)
        zright = -sqrt(-2d0*nrg)*exp(-sqrt(-2d0*nrg)*xinf)

        step = (xinf - xcrs) / prec

        do i = 0, prec
            x = xinf - i*step

            k(1) = fptr(x, yright, zright, nrg)
            l(1) = gptr(x, yright, zright, nrg)

            k(2) = fptr(x-step/2, yright-k(1)*step/2, zright-k(1)*step/2, nrg)
            l(2) = gptr(x-step/2, yright-l(1)*step/2, zright-l(1)*step/2, nrg)

            k(3) = fptr(x-step/2, yright-k(2)*step/2, zright-k(2)*step/2, nrg)
            l(3) = gptr(x-step/2, yright-l(2)*step/2, zright-l(2)*step/2, nrg)

            k(4) = fptr(x-step, yright-k(3)*step, zright-k(3)*step, nrg)
            l(4) = gptr(x-step, yright-l(3)*step, zright-l(3)*step, nrg)

            yright = yright - step*(k(1)+2*k(2)+2*k(3)+k(4))/6
            zright = zright - step*(l(1)+2*l(2)+2*l(3)+l(4))/6

        end do

    end subroutine integR

    subroutine drawCND()
    ! tested
        implicit none
        real(8) :: erg, stp

        stp = (einf-enou)/dots
        open(unit = 1, file = "cndPlot.dat")
        do i = 1, dots+1
            erg = -einf + (i-1)*stp
            write(1, *) erg, cndPtr(erg)
        end do
        close(1)

        print *, "Cnd(erg) is drawn"
    end subroutine drawCND

    function getWF(nrg)
        implicit none
        real(8), intent(in) :: nrg
        real(8), dimension(0:2*prec,2) :: getWF

        abstract interface
            real(8) function fct(x, y, z, nrg)
            real(8), intent(in) :: x, y, z, nrg
            end function fct
        end interface

        procedure(fct), pointer :: fptr, gptr
        real(8), dimension(4) :: k, l
        real(8) :: stp, zstore, ystore, C
        integer :: i

        fptr => func
        gptr => gunc

    ! Calculate left-side WF
        getWF(0,2) = 0
        zstore = 1

        stp = abs(xcrs - xnou) / prec
        do i = 0, prec
            getWF(i,1) = xnou + i*stp

            k(1) = fptr(getWF(i,1), getWF(i,2), zstore, nrg)
            l(1) = gptr(getWF(i,1), getWF(i,2), zstore, nrg)

            k(2) = fptr(getWF(i,1)+stp/2, getWF(i,2)+k(1)*stp/2, zstore+k(1)*stp/2, nrg)
            l(2) = gptr(getWF(i,1)+stp/2, getWF(i,2)+l(1)*stp/2, zstore+l(1)*stp/2, nrg)

            k(3) = fptr(getWF(i,1)+stp/2, getWF(i,2)+k(2)*stp/2, zstore+k(2)*stp/2, nrg)
            l(3) = gptr(getWF(i,1)+stp/2, getWF(i,2)+l(2)*stp/2, zstore+l(2)*stp/2, nrg)

            k(4) = fptr(getWF(i,1)+stp, getWF(i,2)+k(3)*stp, zstore+k(3)*stp, nrg)
            l(4) = gptr(getWF(i,1)+stp, getWF(i,2)+l(3)*stp, zstore+l(3)*stp, nrg)

            if (i.ne.prec) getWF(i+1,2) = getWF(i,2) + stp*(k(1)+2*k(2)+2*k(3)+k(4))/6d0
            zstore = zstore + stp*(l(1)+2*l(2)+2*l(3)+l(4))/6d0

        end do

        ystore = getWF(prec,2)

    ! Calculate right-side WF
        getWF(2*prec,2) = exp(-sqrt(-2d0*nrg)*xinf)
        zstore = -sqrt(-2d0*nrg)*exp(-sqrt(-2d0*nrg)*xinf)

        stp = abs(xinf - xcrs) / prec
        do i = 2*prec, prec, -1
            getWF(i,1) = xinf - (2*prec-i)*stp

            k(1) = fptr(getWF(i,1), getWF(i,2), zstore, nrg)
            l(1) = gptr(getWF(i,1), getWF(i,2), zstore, nrg)

            k(2) = fptr(getWF(i,1)-stp/2, getWF(i,2)-k(1)*stp/2, zstore-k(1)*stp/2, nrg)
            l(2) = gptr(getWF(i,1)-stp/2, getWF(i,2)-l(1)*stp/2, zstore-l(1)*stp/2, nrg)

            k(3) = fptr(getWF(i,1)-stp/2, getWF(i,2)-k(2)*stp/2, zstore-k(2)*stp/2, nrg)
            l(3) = gptr(getWF(i,1)-stp/2, getWF(i,2)-l(2)*stp/2, zstore-l(2)*stp/2, nrg)

            k(4) = fptr(getWF(i,1)-stp, getWF(i,2)-k(3)*stp, zstore-k(3)*stp, nrg)
            l(4) = gptr(getWF(i,1)-stp, getWF(i,2)-l(3)*stp, zstore-l(3)*stp, nrg)

            if (i.ne.prec) getWF(i-1,2) = getWF(i,2) - stp*(k(1)+2*k(2)+2*k(3)+k(4))/6
            zstore = zstore - stp*(l(1)+2*l(2)+2*l(3)+l(4))/6

        end do

    ! Raw data continuity fixing
        C = getWF(prec,2)/ystore
        do i = 0, prec-1
            getWF(i,2) = getWF(i,2)*C
        end do

    ! Raw data normalization fixing
        C = discrItgr(getWF**2)
        do i = 0, 2*prec
            getWF(i,2) = getWF(i,2)/C
        end do

        print*, "WF is retrieved"
    end function getWF

    real(8) function discrItgr(farr)
    ! tested
        implicit none
        real(8), dimension(0:2*prec,2), intent(in) :: farr

        real(8) :: sres, er, next, stp
        integer :: i

        discrItgr = 0
        er = 0
        if (mod(prec, 3).eq.0) then

            stp = farr(1,1) - farr(0,1)
            do i = 0, prec-3, 3
                next = (3*stp/8)*(farr(i,2)+3*farr(i+1,2)+3*farr(i+2,2)+farr(i+3,2)) - er
                sres = discrItgr + next
                er = (sres - discrItgr) - next
                discrItgr = sres
            end do

            stp = farr(prec+1,1) - farr(prec,1)
            do i = prec, 2*prec-3, 3
                next = (3*stp/8)*(farr(i,2)+3*farr(i+1,2)+3*farr(i+2,2)+farr(i+3,2)) - er
                sres = discrItgr + next
                er = (sres - discrItgr) - next
                discrItgr = sres
            end do

        elseif (mod(prec, 3).eq.1) then

            stp = farr(1,1) - farr(0,1)
            do i = 0, prec-7, 3
                next = (3*stp/8)*(farr(i,2)+3*farr(i+1,2)+3*farr(i+2,2)+farr(i+3,2)) - er
                sres = discrItgr + next
                er = (sres - discrItgr) - next
                discrItgr = sres
            end do
            do i = prec-4, prec-2, 2
                next = (stp/3)*(farr(i,2)+4*farr(i+1,2)+farr(i+2,2)) - er
                sres = discrItgr + next
                er = (sres - discrItgr) - next
                discrItgr = sres
            end do

            stp = farr(prec+1,1) - farr(prec,1)
            do i = prec, 2*prec-7, 3
                next = (3*stp/8)*(farr(i,2)+3*farr(i+1,2)+3*farr(i+2,2)+farr(i+3,2)) - er
                sres = discrItgr + next
                er = (sres - discrItgr) - next
                discrItgr = sres
            end do
            do i = 2*prec-4, 2*prec-2, 2
                next = (stp/3)*(farr(i,2)+4*farr(i+1,2)+farr(i+2,2)) - er
                sres = discrItgr + next
                er = (sres - discrItgr) - next
                discrItgr = sres
            end do

        elseif (mod(prec, 3).eq.2) then

            stp = farr(1,1) - farr(0,1)
            do i = 0, prec-5, 3
                next = (3*stp/8)*(farr(i,2)+3*farr(i+1,2)+3*farr(i+2,2)+farr(i+3,2)) - er
                sres = discrItgr + next
                er = (sres - discrItgr) - next
                discrItgr = sres
            end do
            next = (stp/3)*(farr(prec-2,2)+4*farr(prec-1,2)+farr(prec,2)) - er
            sres = discrItgr + next
            er = (sres - discrItgr) - next
            discrItgr = sres

            stp = farr(prec+1,1) - farr(prec,1)
            do i = prec, 2*prec-5, 3
                next = (3*stp/8)*(farr(i,2)+3*farr(i+1,2)+3*farr(i+2,2)+farr(i+3,2)) - er
                sres = discrItgr + next
                er = (sres - discrItgr) - next
                discrItgr = sres
            end do
            next = (stp/3)*(farr(2*prec-2,2)+4*farr(2*prec-1,2)+farr(2*prec,2)) - er
            sres = discrItgr + next
            er = (sres - discrItgr) - next
            discrItgr = sres

        end if
    end function discrItgr

end program

module search
    implicit none
    private

    abstract interface
        real(8) function func(x)
        real(8), intent(in) :: x
        end function func
    end interface

! Root search limiter
    real(8), parameter :: eps = 1d-9 ! 1d-16 is machine epsilon

    public :: getRoot

contains

    subroutine getRoot(fptr, bot, top, arr)
        implicit none

        procedure(func), pointer, intent(in) :: fptr
        real(8), intent(in) :: bot, top
        real(8), dimension(:), intent(inout), allocatable :: arr
        real(8) :: md, bt, tp
        integer :: j

        if (fptr(bot)*fptr(top).gt.0) return
        call ext(arr)

        bt = bot
        tp = top
        md = (top + bot) / 2
        j = 0
        do while ((abs(fptr(md)).gt.eps).and.((tp-bt).gt.eps))
            md = (tp + bt) / 2
            if (fptr(bt) * fptr(md) .gt. 0) then
                bt = md
            else if (fptr(bt) * fptr(md) .lt. 0) then
                tp = md
            end if
            j = j + 1
        end do

        arr(size(arr)) = md
        print *, "(k =", size(arr), ") E =", md, j , "steps"
    end subroutine getRoot

    subroutine ext(arr)
    ! tested
        implicit none
        real(8), dimension(:), intent(inout), allocatable :: arr

        real(8), dimension(size(arr)) :: tmp
        integer :: j, k

        k = size(arr)
        do j = 1, k
            tmp(j) = arr(j)
        end do

        deallocate(arr)
        allocate(arr(k+1))

        do j = 1, k
            arr(j) = tmp(j)
        end do
    end subroutine ext

end module search

#if 0
module draw
    implicit none
    private

    public:

contains



end module
#endif
