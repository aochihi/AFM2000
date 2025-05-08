# 3D-BIEM - version AFM2000
3D boundary integral equation method (BIEM) for earthquake dynamic rupture simulation by Aochi et al. (Pageoph, 2000). A planar fault in a 3D infinite, homogeneous medium is assumed. The code is parallelized using MPI. 

# Example 1 (illapel-mpi-distrib.f90) = planar fault

All the input parameters are written directly in the code. You may modify the variables within the indicated parts (*1, *2, *3). In the given example, 

*1 : Fault dimension (X, Y) = (-20:70, -120: 20) km. Note that the fault slip is always in X positive. A square element size is set as ds = 2 km, so that (90/2 + 1) x (140/2 + 1) = 46 x 71 = 3266 elements. A time step is automatically set as dt = ds/(2 Vp). The code calculates up to 400 time steps. 

*2 : Output directory is defined in "dir". The rigidity and Vp of the homogeneous medium is set as "mu" and "alpha". The slip-weakening law is assumed with the standard values of peak strength, residual stress, and critical slip distance as "tp0", "tr0", "dc0". The initial shear stress is by defulat set as "t0" and initial crack size as "r0". Note that these are the default values and we can introduce any heterogeniety on the fault later. 

*3 : the fault heterogeneity is introduced. Note that the the rupture initiates suddenly at t = 0 around the hypocenter (xhypo, yhypo). Stress drops instantaneously from "t0" to "tr0" independently the friction law in a circular crack of radius "r0". 



## COMPILE
On standard computing centers, you may compile as 
> mpif90 -O illapel-mpi-distrib.f90 kernel31s_05Avril.f


## RUN
No Input file. On the same directory, make sure that the output directory defind by "dir" exists. In the given example, 
> make test1

Set number of MPI processes upto hundreds to 1000 (but should be less than the total element number=3266), and run.  

## OUTPUTS
The code writes some messages as standard out. All the results are saved under the directory "dir", written by the master CPU. 

### hoge2.txt: (input parameters saved)

Line 1: ixmn, ixmx, iymn, iymx, itmx : Grid numbers and time steps for simulations

Line 2: tp0, tr0, t0: Initial reference conditions for rupture.

Line 3: ds, dt, dc0, r0: grid size, time step, controling parameters

Line 4: mu, alpha

Line 5: xhypo, yhypo

Line 6: x0, y0: The largest heterogeneity position. 

### struct.dat: (initial condition saved)

Line 1-end: Element number, X, Y, initial shear stress, peak strength, residual stress, dc 

### output???.dat : output file named by time step

Line 1-end: X (km), Y (km), slip rate (m/s), fault slip (m), shear stress (MPa), strength (MPa)

## For dveloppers

In this code, we evalute direcly the spatio-temporal convolution of the BIEM. We use only tau_{31} component (external function ker31s), as the fault slip is supposed in X axis on a planar fault. In such configuration, the kernel is unique. We can apply FFT for the fast convolution (See IA2005 code). For non-planar fault geometry, we need to combine the other components tau_{ij} and estimate shear and normal stresses every element. (6 December 2023)

# Example 2 ( ) = non-planar vertical fault


# REFERENCE 
Please cite

Aochi, H., E. Fukuyama and M. Matsu'ura (2000). Spontaneous Rupture Propagation on a Non-planar Fault in 3D Elastic Medium, Pure appl. Geophys., 157, 2003-2027. https://doi.org/10.1007/PL00001072

## Application examples

Aochi, H. and S. Ruiz (2021). Early stage and main ruptures of the 2015 Mw8.3 Illapel, Chile, megathrust earthquake : Kinematic elliptical inversions and dynamic rupture simulations, J. Geophys. Res., 126, e2020JB021207. https://doi.org/10.1029/2020JB021207

Ide, S. and H. Aochi (2013). Historical seismicity and dynamic rupture process of the 2011 Tohoku-Oki earthquake, Tectonophys., 600, 1-13. https://doi.org/10.1016/j.tecto.2012.10.018

Aochi, H. and S. Ide (2011). Conceptual multi-scale dynamic rupture model for the 2011 Off-the-Pacific-Coast-of-Tohoku earthquake, Earth Planets Space, 63, 761-765. https://doi.org/10.5047/eps.2011.05.008
