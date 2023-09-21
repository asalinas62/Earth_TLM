# Earth_TLM
The package **Earth_TLM** is a full 3D code for study the electromagnetic wave propagation in planetary atmospheres. 
Using this code, Schumann resonances can be obtained from the transient electromagnetic field.

The physical system is modeled as a spherical cavity formed by two very large and concentric perfectly-conducting spheres: the Earth’s or planet surface with radius $R_1$ and the lower Ionosphere with radius $R_2$. The cavity is filled with the atmosphere, an isotropic and inhomogeneous dielectric with low conductivity values. The system is driven by current sources which model lightning strokes. The sources present a Gaussian-shape time evolution, with important content at low frequencies, below 100 Hz in the Earth’s case. Cartesian coordinates will be used, and the system is included inside a cubic TLM mesh with a side of $2\cdot R_2$.

**The package Earth_TLM has the following programs**:
- ‘input_Earth_box.py’: conversion of spherical coordinates referred to the Earth’s surface to Cartesian coordinates referred to the cubic Cartesian TLM mesh. It reads the file ‘input_Earth.txt’ and writes in the file ‘input_box.txt’
- ´Earth_TLM_MoiT.f95´: obtention of the electric and magnetic field at each node and time iteration. This block is parallelized using OpenMP directives. Output information is provided in two files:
  - voltages_box: a binary file containing voltages and currents at the parallel and series nodes of each 3D TLM output nodes, respectively. These magnitudes correspond to the electric and magnetic fields in arbitrary units.
  - time_box: a formatted file containing the total CPU time spent in the calculation, that, in the examples presented in the paper and using a laptop with i7 processor, is about several hours.

The Jupyter notebook ‘Case_2E.ipynb’ contains the post-processing analysis of the results obtained from ´Earth_TLM_MoiT.f95´.

The Fortran compiler used is **gfortran** (gfortran compiler is free to download and use). **OpenMP directives** are included in the gfortran compiler. **Python3** and **Jupyter notebooks** should also be installed in the system for the pre-processing and post-processing conversion and analysis.

**Compiling and linking the fortran codes**:

`gfortran -o main Earth_TLM_main.f95`

`gfortran -o temporal -O2 -fopenmp Earth_TLM_MoiT.f95`

An example of ‘input_Earth.txt’ file, with comments, is shown below:

| <!-- -->      | <!-- -->        |
|:-------------:|:---------------:|
| 1024        | # Number of timestep calculations       | 
| 10.0E3	         | # Side length of the cubic TLM node        | 
| 6370.E3			|  # Earth’s radius    |
| 6470.E3			|  # Inner Ionosphere’s radius (outer cavity radius)|
| 10				  |  # Conductivity code number |
| 1				    |  # Number of source points |
| 6.E3 90. 0. |	# Source coordinates referred to the Earth’s surface: altitude, zenithal and azimuthal angles in degrees  |
| 1. 0. 0.    |  # Orientation of the current density at each source point, expressed in spherical components |
| 1				    |  # Number of outputs points  |
| 6.E3 90. 0.2 |	# Coordinates of the output points referred to the Earth’s surface: altitude, zenithal and azimuthal angles in degrees |
	
This example is a shortened version of the Case 10 given in the repository that contains all the cases studied in the paper. 
