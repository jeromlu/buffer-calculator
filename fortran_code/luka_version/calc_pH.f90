    program calc_pH_values
    ! This program calculates the grams of each species required to get a certain pH, respecting certain constraints!

    ! Defining modules used
    use read_parameters
    use calc_A_B

    ! Defining variables
    implicit none
    character (len = 11) :: file_version, current_ver
    logical :: file_exists
    integer :: loop
    ! Source: pH transitions in cation exchange chromatographic columns containing
    ! weak acid groups by Carta

    ! Checking for version
    current_ver = "RK190215_05"

    ! Seeing if a version file exists
    INQUIRE(FILE="version.dat", EXIST=file_exists)
    if (file_exists) then
        ! Checking if right version!
        open(unit=20, file='version.dat', status='old', action='read')
        read(20,*) file_version
        close(20)
        if (file_version .eq. current_ver ) then ! if the right version
            ! Calling subroutines
            call read_all
            call calc_buffs
        else
            do loop = 1,20 ! if wrong version
                write (*,*) " Please download the most current version from "
                write (*,*) "\\PHCHBS-WL9525\Tools"
                write (*,*)
            end do
            write (*,*) "This program will terminate when you press enter"
            pause
        end if
    else
        do loop = 1,20 ! if version file does not exist
                write (*,*) " Please download the most current version from "
                write (*,*) "\\PHCHBS-WL9525\Tools"
                write (*,*)
        end do
        write (*,*) "This program will terminate when you press enter"
        pause
    end if



    end program calc_pH_values

