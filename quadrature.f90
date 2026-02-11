!-----------------------------------------------------------------------
!Module: quadrature
!-----------------------------------------------------------------------
!! By Rodrigo Navarro Perez and Dustin Wheeler
!!
!! This module contains methods to numerically calculate integrals. 
!!----------------------------------------------------------------------
!! Included subroutines:
!!
!! 
!! monte_carlo_quad
!!----------------------------------------------------------------------
!! Included functions:
!!
!! real_booles
!! complex_booles
!! booles_rule
!-----------------------------------------------------------------------
!! Included interfaces:
!!
!! func
!! booles_quadrature
!-----------------------------------------------------------------------
module quadrature
use types

implicit none

private
public :: booles_quadrature, monte_carlo_quad

!-----------------------------------------------------------------------
!Interface: func
!-----------------------------------------------------------------------
!! This defines a new type of procedure in order to allow callbacks
!! in the Monte Carlo quadrature subroutine of an arbitrary function that is given
!! as input and declared as a procedure
!!
!! The arbitrary function receives two rank 1 arrays of arbitrary size.
!! The first array contains an n-dimensional vector representing the
!! point sampled by the Monte Carlo method. The second is a "work array"
!! that contains parameters  necessary to calculate the function to be
!! integrated.
!!----------------------------------------------------------------------
interface
    real(dp) function func(x, data)
        use types, only : dp
        implicit none
        real(dp), intent(in) :: x(:), data(:)
        ! This is the interface we saw in class that allows callbacks
    end function func
end interface

!-----------------------------------------------------------------------
!Interface: booles_quadrature
!-----------------------------------------------------------------------
!! This interface allows one to integrate real or complex valued functions 
!! over the real x axis, instead of just real ones.
!!
!!----------------------------------------------------------------------
interface booles_quadrature
    module procedure real_booles, complex_booles
end interface booles_quadrature

 
contains

!-----------------------------------------------------------------------
!! Function: real_booles
!-----------------------------------------------------------------------
!! By: Rodrigo Navarro Perez and Dustin Wheeler
!!
!! Use's Boole's method to numerically calculate the integral of f(x) from
!! xa to xb using a given array of fn's where fn = f(xn) from fa to fb.
!! ----------------------------------------------------------------------
!! Input:
!!
!! fx           real        Array containing the evaluated function
!! delta_x      real        Distance between the evaluation points
!!----------------------------------------------------------------------
!! Output:
!!
!! s            real        Result of the Boole's quadrature
!-----------------------------------------------------------------------
real(dp) function real_booles(fx, delta_x) result(s)
    implicit none
    real(dp), intent(in) :: fx(1:), delta_x

    integer :: fx_size, i
	integer :: num_intervals
	complex(dp) :: temp_sum
	
    fx_size = size(fx)

    ! As the diagram below shows, only certain number of grid points
    ! fit the scheme of Boole's quadrature. Implement a test 
    ! to make sure that the number of evaluated points in the fx array
    ! is the correct one

    ! |--interval 1---|--interval 2---|--interval 3---|
    ! 1   2   3   4   5   6   7   8   9   10  11  12  13
    ! |---|---|---|---|---|---|---|---|---|---|---|---|
    ! x0  x1  x2  x3  x4
    !                 x0  x1  x2  x3  x4
    !                                 x0  x1  x2  x3  x4

	! bins=fx_size-1 must be a multiple of 4, so mod(fx_size-1,4)=0
    if(mod(fx_size,4) /= 1) then
		print *, 'ERROR::quadrature:booles_quadrature: fx array size in booles_quadrature ', &
			'has to be a multiple of 4 plus 1. Currenty, however, it is a multiple of 4 plus ', &
			mod(fx_size,4)
        stop
    endif

	num_intervals = (fx_size-1)/4
    ! We could implement the full integration here, however to make a cleaner,
    ! easy to read (and debug or maintain) code we will define a smaller
    ! function that returns Boole's five point rule and pass slices (1:5), (5:9),
    ! (9:13), ... of fx to such function to then add all the results. 

    temp_sum = 0._dp
    do i = 0, num_intervals-1
		temp_sum = temp_sum + booles_rule(dcmplx(fx(i*4+1: i*4+5),0.0_dp),delta_x)
    enddo
	
	s = dble(temp_sum)
end function real_booles

!-----------------------------------------------------------------------
!! Function: complex_booles
!-----------------------------------------------------------------------
!! By: Rodrigo Navarro Perez and Dustin Wheeler
!!
!! Same as real_booles, but with complex valued functions.
!! ----------------------------------------------------------------------
!! Input:
!!
!! fx           complex        Array containing the evaluated function
!! delta_x      real           Distance between the evaluation points
!!----------------------------------------------------------------------
!! Output:
!!
!! s            complex        Result of the Boole's quadrature
!-----------------------------------------------------------------------
complex(dp) function complex_booles(fx, delta_x) result(s)
    implicit none
	complex(dp), intent(in) :: fx(1:)
    real(dp), intent(in) :: delta_x

    integer :: fx_size, i
	integer :: num_intervals
	
    fx_size = size(fx)

	! bins=fx_size-1 must be a multiple of 4, so mod(fx_size-1,4)=0
    if(mod(fx_size,4) /= 1) then
		print *, 'ERROR::quadrature:booles_quadrature: fx array size in booles_quadrature ', &
			'has to be a multiple of 4 plus 1. Currenty, however, it is a multiple of 4 plus ', &
			mod(fx_size,4)
        stop
    endif

	num_intervals = (fx_size-1)/4

    s = 0._dp
    do i = 0, num_intervals-1
		s = s + booles_rule(fx(i*4+1: i*4+5),delta_x)
    enddo
	
end function complex_booles

!-----------------------------------------------------------------------
!! Function: booles_rule
!-----------------------------------------------------------------------
!! By: Rodrigo Navarro Perez and Dustin Wheeler
!!
!! Use's Boole's method to numerically calculate the integral of f(x) from
!! x_a to x_(a+4) using a given array of f_n's where f_n = f(x_n) from f_a 
!! to f_(a+4).
!! ----------------------------------------------------------------------
!! Input:
!!
!! fx           complex     Array containing the evaluated function
!! delta_x      real        Distance between the evaluation points
!!----------------------------------------------------------------------
!! Output:
!!
!! s            complex      Result of the Boole's quadrature
!-----------------------------------------------------------------------
complex(dp) function booles_rule(fx, delta_x) result(s)
    implicit none
    complex(dp), intent(in) :: fx(1:)
	real(dp), intent(in) :: delta_x

    integer :: fx_size
    complex(dp) :: fx0, fx1, fx2, fx3, fx4

    fx_size = size(fx)

    ! Let's make an additional test to make sure that the array
    ! received has 5 and only 5 points 

    if(fx_size/=5) then
		print *, 'ERROR::quadrature:booles_rule: Received array has ',fx_size, &
			'elements, should have exactly 5.'
		stop
    endif
    
    fx0 = fx(1)
    fx1 = fx(2)
    fx2 = fx(3)
	fx3 = fx(4)
    fx4 = fx(5)
    
    s = delta_x*(2.0/45)*(7*fx4 + 32*fx3 + 12*fx2 + 32*fx1 + 7*fx0)
end function booles_rule

!-----------------------------------------------------------------------
!! Subroutine: monte_carlo_quad
!-----------------------------------------------------------------------
!! By: Rodrigo Navarro Perez and Dustin Wheeler
!!
!! Generates an n_samples amount of random positions, calculates the given 
!! function f at said random positions to generate an array of values, then
!! uses the average of those values to approximate the integral of the given 
!! function f in the rectangular prism bounds bounds of the volume that 
!! bounds the generated random positions.
!! ----------------------------------------------------------------------
!! Input:
!!
!! f            procedure   function to be integrated
!! a            real        array containing the lower limits of the integral
!! b            real        array containing the upper limits of the integral
!! data         real        array containing parameters necessary to calculate the function f
!! n_samples    integer     number of sample points in the Monte Carlo integration
!!----------------------------------------------------------------------
!! Output:
!!
!! s            real        Result of the Monte Carlo integral
!! sigma_s      real        Estimate of the uncertainty in the Monte Carlo integral
!-----------------------------------------------------------------------
subroutine monte_carlo_quad(f, a, b, data, n_samples, s, sigma_s)
    implicit none
    procedure(func) :: f
    real(dp), intent(in) :: a(:), b(:), data(:)
    integer, intent(in) :: n_samples
    real(dp), intent(out) :: s, sigma_s
	
    integer :: i, vector_size
    real(dp), allocatable :: x_vector(:), fx(:)! ...you might need to declare other arrays here
	
	real(dp) :: volume, fx_sum, fx2_sum, variance
	
    vector_size = size(a)

    ! We're defining a Monte Carlo routine that works for an arbitrary number of 
    ! dimensions in the integral (Remember, that's the advantage of Monte Carlo integration,
    ! it's very efficient for high dimensional integrals)

    ! Since a and b give the lower and upper limits they need to have the same size.
    ! Make a check to see if they do have the same size

    if(vector_size /= size(b)) then
		print *, 'a and b arrays in monte_carlo_quad have to be the same size'
		stop       
    endif
	
	! Calculate the volume of the rectangular prism shaped reactor (d*w*h)
	volume = 1.0_dp
	do i = 1,vector_size
		volume = volume * abs(b(i)-a(i))
	enddo
	
    ! Here we allocate memory for the vector containing the sample points and 
    ! for a vector that contains the evaluated function
    allocate(x_vector(1:vector_size))
    allocate(fx(1:n_samples))

    do i=1,n_samples
        call random_number(x_vector) !generates an array with random numbers in the [0,1) interval
        x_vector = a + x_vector*(b-a) !rescaling to the integration volume [a,b)
        fx(i) = f(x_vector,data)
    enddo
	
	!--------------------------------------------------------------------
	! s = Integral ~= Volume*<f> = Volume/Samples * (sum of fx(n) over n)
	!--------------------------------------------------------------------
	!! fx_sum = sum of fx(n) over n
	fx_sum = 0.0_dp
	do i=1,n_samples
		fx_sum = fx_sum + fx(i)
	enddo
	!! scale fx_sum by volume/samples to get approximation for integral
	s = fx_sum*volume/n_samples
	
	!--------------------------------------------------------------------
	! sigma_s = Uncertainty/Statistical error in integral = (Volume/sqrt(Samples)) * sqrt(variance)
	!--------------------------------------------------------------------
	!! variance = <f^2> - <f>^2; need <f^2> = 1/samples * (sum of fx(n)^2 over n)
	!! fx2_sum = sum of fx(n)^2 over n
	fx2_sum = 0.0_dp
	do i=1,n_samples
		fx2_sum = fx2_sum + fx(i)**2
	enddo
	!! calc variance then full integral error
	variance = fx2_sum/n_samples - (fx_sum/n_samples)**2
	sigma_s = (volume/sqrt(dble(n_samples))) * sqrt(variance)
	
end subroutine monte_carlo_quad

end module quadrature