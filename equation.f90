module equation
    implicit none
    private

    abstract interface
        pure real(8) function ufct(x)
        real(8), intent(in) :: x
        end function ufct
    end interface

! Energy effective zero
    real(8), parameter :: enou = 4d-3 ! >0
! Energy effective infinity
    real(8), parameter :: einf = 150d0 ! >0
! Energy search precision
    integer, parameter :: edots = 6000

! WF effective zero
    real(8), parameter :: xnou = 7d-10
! WF effective infinity
    real(8), parameter :: xinf = 2.5d0
! Cross-linking point (for solution-finding algorythm)
    real(8), parameter :: xcrs = 0.108d0
    ! 0.338d0 for g=2d0; 0.108d0 perfect for l=0

! Left boundary WF
    real(8), parameter :: yleft = 0d0
! Left boundary WF derivative
    real(8), parameter :: zleft = 1d0

! Using potential
    procedure(Yl), pointer, public :: uptr => Yl

    public :: func, gunc, yleft, zleft, yright, zright, enou, einf, edots, xnou, xinf, xcrs
contains

! Right boundary WF
    pure real(8) function yright(nrg)
        implicit none
        real(8), intent(in) :: nrg
        yright = dexp(-dsqrt(-2d0*nrg)*xinf)
    end function yright
! Right boundary WF derivative
    pure real(8) function zright(nrg)
        implicit none
        real(8), intent(in) :: nrg
        zright = -dsqrt(-2d0*nrg)*dexp(-dsqrt(-2d0*nrg)*xinf)
    end function zright

    pure real(8) function Yl(r)
        implicit none
        real(8), intent(in) :: r
        real(8), parameter :: g = 40d0
        integer, parameter :: l = 0
        Yl = -g*exp(-r)/r + 0.5d0*l*(l+1)/r**2
    end function Yl

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

end module equation
