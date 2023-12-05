# AFM2000
3D boundary integral equation method (BIEM) for earthquake dynamic rupture simulation by Aochi et al. (Pageoph, 2000)

# COMPILE
On standard computing centers, you may compile as 
> mpif90 -O illapel-mpi-omp-distrib.f90 kernel31s_05Avril.f

if necesary, precise the options
> mpif90 -O -openmp -lmpi illapel-mpi-omp-distrib.f90 kernel31s_05Avril.f

# RUN
On the same directory, make sure that the output directory exists. In the given example, 
> make test1

Set number of MPI processes upto hundreds to 1000, and run.  

# REFERENCE 
Please cite

Aochi, H., E. Fukuyama and M. Matsu'ura (2000). Spontaneous Rupture Propagation on a Non-planar Fault in 3D Elastic Medium, Pure appl. Geophys., 157, 2003-2027. https://doi.org/10.1007/PL00001072
