!-----------------------------------------------------------------------
!! Module: potential
!-----------------------------------------------------------------------
!! By: Malida Hecht, Dustin Wheeler, Zach Barvian
!!
!! This module contains the subroutines and functions needed to compute
!! the electric potential from the poisson equation.
!!----------------------------------------------------------------------
!! Included subroutines:
!!
!! calculate_coefficients
!! parallel_calculate_coefficients
!!----------------------------------------------------------------------
!! Included functions:
!!
!! phi_potential
!! phi_potential_parallel
!! source_func_rho
!-----------------------------------------------------------------------
module potential
use types
use quadrature, only : booles_quadrature
implicit none

private
public :: calculate_coefficients, phi_potential, parallel_calculate_coefficients, phi_potential_parallel


contains

!-----------------------------------------------------------------------
!! Subroutine: calculate_coefficients
!-----------------------------------------------------------------------
!! By: Malida Hecht, Dustin Wheeler
!!
!! Calculates the fourier coefficients of the phi potential function
!!
!!----------------------------------------------------------------------
!! Input:
!!
!! box_length	real	array containing the dimensions of the box that contains the problem
!! rho_zero		real	constant coefficient of the source term
!! center		real	coordinates of the centerpoint of the source term
!! rho_width	real	parameter that quantitatively describes how far the exponential 
!!							source term 'reaches'
!! n_max	  integer	number of fourier samples per dimensions
!! n_sample	  integer 	number of spatial samples per dimension to be used in integrating over space
!-----------------------------------------------------------------------
!!----------------------------------------------------------------------
!! Output:
!!
!! c_fourier	real	fourier coefficients of the phi_potential function
!-----------------------------------------------------------------------
subroutine calculate_coefficients(box_length, rho_zero, center, rho_width, n_max, c_fourier,n_sample)
    implicit none
    real(dp), intent(in) :: box_length(1:2), rho_zero, center(1:2), rho_width(1:2)
    integer, intent(in) :: n_max, n_sample
    real(dp), allocatable, intent(out) :: c_fourier(:,:)
    integer :: m, n, i_x, i_y
    real(dp):: x, y, dx, dy, A, rho_mn
    real(dp), allocatable :: g_xy(:), f_x(:)

    ! Allocate
    allocate(c_fourier(1:n_max, 1:n_max))
    allocate(g_xy(1:n_sample))
    allocate(f_x(1:n_sample))

    ! Set up constant values
    A = 2._dp / SQRT( box_length(1) * box_length(2) )
    dx = box_length(1) / n_sample
    dy = box_length(2) / n_sample

    ! Integrate for c_mn
    do m = 1, n_max
        do n = 1, n_max
        
        !Set Up double integral integral for each m,n.
            do i_x = 1, n_sample
                x = ( i_x - 1) * dx
                
                do i_y = 1, n_sample
                    y = ( i_y - 1) * dy
                    g_xy(i_y) = source_func_rho(x,y, rho_zero, center, rho_width) * cos(m* pi* x / box_length(1)) * cos(n * pi * y &
                    
                            / box_length(2))
                end do
                
                f_x(i_x) = booles_quadrature(g_xy,dy)
            end do
            
            rho_mn = A  * booles_quadrature(f_x, dx)
            
           
            c_fourier(m, n) = (4._dp * pi) * rho_mn / ((m*pi/box_length(1))**2 + (n*pi/box_length(2))**2)
        enddo
    enddo

end subroutine calculate_coefficients

!-----------------------------------------------------------------------
!! Subroutine: parallel_calculate_coefficients
!-----------------------------------------------------------------------
!! By: Malida Hecht, Dustin Wheeler, Zach Barvian
!!
!! Calculates the fourier coefficients of the phi potential function
!! but does so in parallel.
!!
!!----------------------------------------------------------------------
!! Input:
!!
!! box_length	real	array containing the dimensions of the box that contains the problem
!! rho_zero		real	constant coefficient of the source term
!! center		real	coordinates of the centerpoint of the source term
!! rho_width	real	parameter that quantitatively describes how far the exponential 
!!							source term 'reaches'
!! n_max	  integer	number of fourier samples per dimensions
!! n_sample	  integer 	number of spatial samples per dimension to be used in integrating over space
!-----------------------------------------------------------------------
!!----------------------------------------------------------------------
!! Output:
!!
!! c_fourier	real	fourier coefficients of the phi_potential function
!-----------------------------------------------------------------------
subroutine parallel_calculate_coefficients(box_length, rho_zero, center, rho_width, n_max, c_fourier, n_sample)
    implicit none
    real(dp), intent(in) :: box_length(1:2), rho_zero, center(1:2), rho_width(1:2)
    integer, intent(in) :: n_max, n_sample
    real(dp), allocatable, intent(out) :: c_fourier(:,:)
    integer :: m, n, i_x, i_y
    real(dp):: x, y, dx, dy, A, rho_mn
    real(dp), allocatable :: g_xy(:), f_x(:)

    ! Allocate
    allocate(c_fourier(1:n_max, 1:n_max))
    allocate(g_xy(1:n_sample))
    allocate(f_x(1:n_sample))

    ! Set up constant values
    A = 2._dp / SQRT( box_length(1) * box_length(2) )
    dx = box_length(1) / n_sample
    dy = box_length(2) / n_sample

!$omp parallel default(none) private(m, n, i_x, i_y, x, y, g_xy, f_x, rho_mn), &
!$omp shared(c_fourier, rho_zero, center, box_length, A, dx, dy, rho_width, n_max, n_sample)
!$omp do schedule(static)
    
    ! Integrate for c_mn
    do m = 1, n_max
        do n = 1, n_max
        
        !Set Up double integral integral for each m,n.
            do i_x = 1, n_sample
                x = ( i_x - 1) * dx
                
                do i_y = 1, n_sample
                    y = ( i_y - 1) * dy
                    g_xy(i_y) = source_func_rho(x,y, rho_zero, center, rho_width) * cos(m* pi* x / box_length(1)) * cos(n * pi * y &
                    
                            / box_length(2))
                end do
                
                f_x(i_x) = booles_quadrature(g_xy,dy)
            end do
            
            rho_mn = A  * booles_quadrature(f_x, dx)
            
           
            c_fourier(m, n) = (4._dp * pi) * rho_mn / ((m*pi/box_length(1))**2 + (n*pi/box_length(2))**2)
        enddo
    enddo

!$omp end do
!$omp end parallel

end subroutine parallel_calculate_coefficients

!-----------------------------------------------------------------------
!! Function: phi_potential
!-----------------------------------------------------------------------
!! By: Malida Hecht
!!
!! This function computes the potential, phi as the solution to the Poisson
!! equation.
!!----------------------------------------------------------------------
!! Input:
!! c_fourier	real	fourier coefficients
!! length		real 	dimensions of the box containing the system
!! x			real	x position
!! y			real 	y position
!!----------------------------------------------------------------------
!! Output:
!! r			real	value of the potential at position (x,y)
!-----------------------------------------------------------------------
real(dp) function phi_potential(c_fourier, length, x, y) result(r)
    implicit none
    real(dp), intent(in) :: c_fourier(:,:), length(1:2), x, y
    integer :: m, n

    r = 0._dp

    do m = 1 , size(c_fourier, 1)
        do n = 1 , size(c_fourier, 2)
            r = r + c_fourier(m,n) * (2 / sqrt(length(1) * length(2))) * cos(m * pi * x / length(1)) * cos(n * pi * y / length(2))
        enddo
    enddo

end function phi_potential

!-----------------------------------------------------------------------
!! Function: phi_potential_parallel
!-----------------------------------------------------------------------
!! By: Malida Hecht, Dustin Wheeler
!!
!! This function computes the potential, phi as the solution to the Poisson
!! equation.
!!----------------------------------------------------------------------
!! Input:
!! c_fourier	real	fourier coefficients
!! length		real 	dimensions of the box containing the system
!! x			real	x position
!! y			real 	y position
!!----------------------------------------------------------------------
!! Output:
!! r			real	value of the potential at position (x,y)
!-----------------------------------------------------------------------
real(dp) function phi_potential_parallel(c_fourier, length, x, y) result(r)
    use OMP_LIB
	implicit none
    real(dp), intent(in) :: c_fourier(:,:), length(1:2), x, y
	real(dp), allocatable :: vals(:)
    integer :: m, n, p, thread_id
	
	allocate(vals(1:OMP_GET_MAX_THREADS()))
	vals = 0.0_dp
	
!$omp parallel default(none) private(m,n,thread_id), shared(vals, c_fourier, length,x,y)
	thread_id = OMP_GET_THREAD_NUM()+1
!$omp do schedule(static)
    do m = 1 , size(c_fourier, 1)
        do n = 1 , size(c_fourier, 2)
			!if(thread_id==0) print*, "hello from thread ",thread_id, "i'm on step", m,n
            vals(thread_id) = vals(thread_id) + c_fourier(m,n) * (2 / sqrt(length(1) * length(2))) * &
				cos(m * pi * x / length(1)) * cos(n * pi * y / length(2))
        enddo
    enddo
!$omp end do
!$omp end parallel

	r = 0.0_dp
	do p = 1,OMP_GET_MAX_THREADS()
		r = r + vals(p)
	enddo
end function phi_potential_parallel



!-----------------------------------------------------------------------
!! Function: source_func_rho
!-----------------------------------------------------------------------
!! By: Dustin Wheeler
!!
!! This function provides the source term rho(x,y)
!!
!!----------------------------------------------------------------------
!! Input:
!! x			real	x position
!! y			real	y position
!! rho_0		real	source term coefficient
!! center		real	position of the center of the source
!! rho_width	real	quantitative description of how far the source reaches
!!----------------------------------------------------------------------
!! Output:
!! p			real 	intensity of the source at position (x,y)
!-----------------------------------------------------------------------
real(dp) function source_func_rho(x,y, rho_0, center, rho_width) result(p)
    implicit none
    real(dp), intent(in) :: x, y, rho_0, center(2), rho_width(2)

	p = rho_0 / (pi * rho_width(1) * rho_width(2)) * &
		exp( -((x-center(1))/rho_width(1))**2 - ((y-center(2))/rho_width(2))**2 )

end function source_func_rho


end module potential