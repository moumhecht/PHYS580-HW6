!-----------------------------------------------------------------------
!Module: read_write
!-----------------------------------------------------------------------
!! By: Malida Hecht, Zach Barvian, Dustin Wheeler
!!
!! This module contains routines that handle most of the program's I/O to 
!! and from console, as well as to file.
!!----------------------------------------------------------------------
!! Included subroutines:
!! 
!! read_input
!! write_potential
!! print_wall_clock_time
!! write_potential_parallel
!!----------------------------------------------------------------------
!! Included functions:
!!
!! 
!-----------------------------------------------------------------------
module read_write
use types
use potential, only : phi_potential, phi_potential_parallel
implicit none

private
public :: read_input, write_potential, print_wall_clock_time, write_potential_parallel

contains

!-----------------------------------------------------------------------
!! Subroutine: read_input
!-----------------------------------------------------------------------
!! By: Malida Hecht
!!
!! Takes in a namelist file as a command argument when running the program,
!! then reads said namelist file to use for calculations.
!!
!!----------------------------------------------------------------------
!! Output:
!! length		real	2 element array containg the dimensions of the box
!! rho_zero		real	coefficient describing the max strength of the source term
!! center		real	the coordinates of the center of the source term
!! width		real 	quantitative description of how far the source term stretches in either direction
!! n_max	   integer	number of fourier terms
!! n_sample	   integer	number of spatial samples per dimension
!!
!! output_file 		character 	string containing name of file containing output of the serial calculations
!! output_parallel 	character 	string containing name of file containing output of the parallel calculations
!-----------------------------------------------------------------------
subroutine read_input(length, rho_zero, center, width, n_max, output_file, n_sample, output_parallel)
    implicit none
    real(dp), intent(out) :: length(1:2)
    real(dp), intent(out) :: rho_zero
    real(dp), intent(out) :: center(1:2), width(1:2)
    integer, intent(out) :: n_max, n_sample
    character(len = 200) :: namelist_file
    character(len =*), intent(out) :: output_file, output_parallel
    integer :: unit1, ios, n_arguments, ierror, file_unit
    logical :: file_exists

    namelist /box/ length
    namelist /charge_distribution/ rho_zero, center, width
    namelist /sampling/ n_max,  n_sample
    namelist /output/ output_file, output_parallel


    print*, "This program solves the poisson equation"
    !write stuff here 
    ! Initialize default values
    length = [2._dp, 2._dp]
    rho_zero = 3._dp
    center = [1._dp, 1._dp]
    width = [4._dp, 4._dp]
    n_max = 10
    n_sample = 1001
    output_file = 'results.dat'
    output_parallel = 'results_parallel.dat'

    ! get namelist file name from command line

    n_arguments = command_argument_count()
  
    IF (n_arguments == 1) THEN
        ! Verify if File exists
        CALL get_command_argument(1,namelist_file)
        INQUIRE(file = TRIM(namelist_file), exist =  file_exists)
    
        ! read namelists
        IF (file_exists) THEN
            OPEN(newunit = file_unit, file = namelist_file)
            ! Read namelist for box dimensions
            READ(file_unit, nml=box, iostat= ierror) 
            
            IF (ierror /= 0) THEN
                PRINT *, "Error Reading box namelist."
                STOP 
            END IF
            ! Read namelist for charge distribution conditions
            READ(file_unit, nml=charge_distribution, iostat = ierror)
            IF (ierror /= 0) THEN
                
                PRINT *, "Error Reading charge_distribution namelist."
                STOP 
            END IF
            !Read Namelist for sampling parameters
            READ(file_unit, nml=sampling, iostat = ierror)
            IF (ierror .NE. 0) THEN
                PRINT *, "Error Reading sampling namelist."
                STOP
            END IF
            
            !Reaad Namelist for output file
            READ(file_unit, nml=output, iostat = ierror)
            IF (ierror .NE. 0) THEN
                PRINT *, "Error Reading output namelist."
                STOP
            END IF
         
            CLOSE(file_unit)
        ELSE
            PRINT *, namelist_file, 'File not found.'
            STOP 
        END IF 
    ELSE IF (n_arguments .NE. 0) THEN
        PRINT *, 'Too many files, the program can only take up to one file only'
        PRINT *, 'See details in README.md'
        STOP 
    END IF



end subroutine read_input

!-----------------------------------------------------------------------
!! Subroutine: write_potential
!-----------------------------------------------------------------------
!! By: Malida hecht
!!
!! This subroutine writes the phi potential at a given  position (x,y)
!! using phi_potential (serial calculation)
!!----------------------------------------------------------------------
!! Input:
!!	
!! filename		character 	string containing the name of the file to be written to
!! c_fourier	real		2d array containing the nth,mth calculated fourier coefficient
!! length		real		2 element array containing the dimensions of the box
!! n_max		integer		number of fourier coefficients per dimension
!! rho_0		real		coefficient and max value of the source term
!! r_width		real		2 element array that quantitatively describes how far the exponential source term "reaches"
!! n_sample		integer		number of spatial gridpoints per dimension
!!----------------------------------------------------------------------
!! Output:
!! N/A
!-----------------------------------------------------------------------
subroutine write_potential(filename, c_fourier, length, n_max, rho_0, r_width,n_sample)
    implicit none

    REAL(dp), INTENT(IN) :: c_fourier(:,:), length(:), rho_0, r_width(:)
    INTEGER, INTENT(IN) :: n_max, n_sample
    CHARACTER(len = *), INTENT(IN) :: filename
    REAL (dp) :: x, y , dx, dy
    INTEGER :: i_x, i_y , unit
    REAL(dp) :: phi    

    OPEN (newunit = unit, file = trim(filename))

    dx = length(1) / n_sample
    dy = length(2) / n_sample


    WRITE(unit, *) "L_x = " , length(1) , "L_y = " , length(2) , "R_x = " , r_width(1), 'R_y = ' , r_width(2), "Rho_0 = ", rho_0
     WRITE(unit, *) "n_max =  " , n_max , "n_sample " , n_sample
    WRITE(unit, '(3a28)') "x", "y", "phi(x,y)"

    do i_x = 1 , n_sample

        x = (i_x - 1) * dx

        do i_y = 1 , n_sample

            y = (i_y - 1) * dy

            phi = phi_potential(c_fourier, length, x, y)
            

            WRITE (unit, '(3e28.6)') x, y, phi

        end do
    end do

    CLOSE(unit)
end subroutine write_potential

!-----------------------------------------------------------------------
!! Subroutine: print_wall_clock_time
!-----------------------------------------------------------------------
!! By: Zach Barvian
!!
!! This subroutine computes the computational time and prints to screen.
!!
!!----------------------------------------------------------------------
!! Input:
!!
!! count_1		integer		initial time 
!! count_2		integer		final time 
!! count_rate	integer		counts per second; lets you convert count_1 and count_2 to seconds										
!-----------------------------------------------------------------------
!!----------------------------------------------------------------------
!! Output:
!! N/A
!-----------------------------------------------------------------------
subroutine print_wall_clock_time(count_1, count_2, count_rate)
    implicit none
    integer, intent(in) :: count_1, count_2, count_rate
    real(dp) :: time
    
    time = (count_2 - count_1)/real(count_rate,kind = dp)

    print *,   time 
    
end subroutine print_wall_clock_time




!-----------------------------------------------------------------------
!! Subroutine: write_potential_parallel
!-----------------------------------------------------------------------
!! By: Malida hecht
!!
!! This subroutine writes the phi potential at a given position (x,y) using
!! phi_potential_parallel
!!
!!----------------------------------------------------------------------
!! Input:
!!
!! filename		character 	string containing the name of the file to be written to
!! c_fourier	real		2d array containing the nth,mth calculated fourier coefficient
!! length		real		2 element array containing the dimensions of the box
!! n_max		integer		number of fourier coefficients per dimension
!! rho_0		real		coefficient and max value of the source term
!! r_width		real		2 element array that quantitatively describes how far the exponential source term "reaches"
!! n_sample		integer		number of spatial gridpoints per dimension
!!----------------------------------------------------------------------
!! Output:
!! N/A
!-----------------------------------------------------------------------
subroutine write_potential_parallel(filename, c_fourier, length, n_max, rho_0, r_width,n_sample)
    implicit none

    REAL(dp), INTENT(IN) :: c_fourier(:,:), length(:), rho_0, r_width(:)
    INTEGER, INTENT(IN) :: n_max, n_sample
    CHARACTER(len = *), INTENT(IN) :: filename
    REAL (dp) :: x, y , dx, dy
    INTEGER :: i_x, i_y , unit
    REAL(dp) :: phi    

    OPEN (newunit = unit, file = trim(filename))

    dx = length(1) / n_sample
    dy = length(2) / n_sample


    WRITE(unit, *) "L_x = " , length(1) , "L_y = " , length(2) , "R_x = " , r_width(1), 'R_y = ' , r_width(2), "Rho_0 = ", rho_0
     WRITE(unit, *) "n_max =  " , n_max , "n_sample " , n_sample
    WRITE(unit, '(3a28)') "x", "y", "phi(x,y)"

    do i_x = 1 , n_sample

        x = (i_x - 1) * dx

        do i_y = 1 , n_sample

            y = (i_y - 1) * dy

            phi = phi_potential_parallel(c_fourier, length, x, y)
            

            WRITE (unit, '(3e28.6)') x, y, phi

        end do
    end do

    CLOSE(unit)
end subroutine write_potential_parallel

end module read_write