program calculateStationary
    use equation
    use search
    use drawWF
    use condition
    implicit none

    real(8) :: ebt, etp, estep
    real(8), dimension(:), allocatable :: stErg
    real(8), dimension(0:2*dots, 2) :: WF
    integer :: i, j
    character(len=22) :: name = "graphics/k??wfData.dat"

! Drawing cnd(erg)
    call getCND()

! Searching stationary levels
    allocate(stErg(0))
    estep = (einf-enou)/edots
    do i = 1, edots
        ebt = -einf + (i-1)*estep
        etp = -einf + i*estep
        call getRoot(cndPtr, ebt, etp, stErg)
    end do

! WF output
    do j = 1, size(stErg)
        WF=getWF(stErg(j))
        write(unit=name(11:12), fmt="(i2.2)") j
        open(unit = 1, file = name)
        do i = 0, 2*dots
            write(1,*) WF(i, 1), WF(i, 2)
        end do
        close(1)
    end do

! Fixing memory leak
    deallocate(stErg)

end program calculateStationary

