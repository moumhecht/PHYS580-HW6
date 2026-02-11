# Solving the Poisson Equation 

 This program numerically solves the 2-D Poisson equation for an electrostatic potential with a defined source charge density in a metal rectangular box with sides of length Lx and Ly. We impose Neumann boundary conditions, requiring the gradient of the potential be zero at the sides of the rectangular box.

 The solution to the potential is given by:

$$(\frac{\partial^2}{\partial x^2} + \frac{\partial^2}{\partial y^2}  \phi(x,y) = - 4 \pi \rho(x,y)$$
This program uses Boole's quadrature to compute the Fourier coefficients. 

file that is specified by the namelist or default `results.dat`.

## Parallelization 

Parallelization was done using OpenMP. We parallelizaed two subroutines: computing the fourier coefficients and calculating the Phi potential:


```
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

```

We also parallelized the phi potential subroutine that sums the fourier coefficients:

```
	
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
```
For this subroutine, each thread had their own value for the sum, which was stored into an array and then summed to contribute to the entire potential.
## Output

## Determination of reasonable n_max


## File Descriptions  (╯°□°)╯︵ ┻━┻
- `main.f90`: calls subroutines to read namelist files for initial conditions, computing fourier coefficients, writing the output to results files, and printing computational time to terminal.
- `read_write.f90`: contains subroutines needed 
- `mechanics.f90`: contains subroutines and functions to describe system of differential equations describing planetary motion, computing total energy and magnitude of total angular momentum
- `ode_solver.f90`: contains subroutines to solve a system of differential equations using the fourth order runge-kutta numerical method
- `plots.ipynb` : that takes namelist file and results from output file and plots orbit, total energy as a function of time, and total angular momentum as a function of time
- `types.f90` : file declares and defines variable types
- `makefile` : code to compile all fortran files

## Namelist Format (╯°□°)╯︵ ┻━┻

The namelist file must have the following format in order to be properly read:

```
&box
    length = 5., 5.
/
&charge_distribution
    rho_zero = 1.,
    center = 0.5 0.5,
    width = 3. 3.
/
&sampling
   n_max = 10,
   n_sample = 1001
/
&output
    output_file = 'results.dat',
    output_parallel = 'output_parallel.dat'
/

```

## Running the Code and Generating Plots (╯°□°)╯︵ ┻━┻

To compile, type in terminal:

`$ make`

The code will then compile and create the program 'poisson'.

To run the program, type in terminal: `$ ./poisson`

If a namelist file exists, type in terminal: `$ ./poisson filename.namelist`
Where the `'filename'` should be changed according to the namelist file the user wants to use.

The program will then write the x and y coordinates of each mass, the total energy and total angular momentum to an output .dat file specified by namelist file.

To view the plots, open 'plots.ipynb'. Be sure to provide the proper results output file.
Once set, hit run to plot.


To clear output files, type in terminal:

`$ make clean`
