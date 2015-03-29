# FEM
## Parallel FEM for linear and non-linear small deformation static problems using PETSC in Fortran 
### PETSC --- Parallel Extensible Toolkit for Scientific Computation (http://www.mcs.anl.gov/petsc/)
  * Requirements
	1. Requires PETSC installation with MUMPS, METIS and PARMETIS
  * Input Parameters
	1. Mesh input file along with Boundary Conditions need to be given 
  * Output 
	vtk files from each processor numbered by rank
	Read using Paraview
	 1. Displacement Vectors at each nodal point
	 2. Stress Tensor, averaged at each node, represented using Scalars (1..6)
  * Implemetation
	1. Elements
	  * 2D Elements 
	    * Plane strain 3-Noded Triangular Isoparametric elements with Linear Interpolation 
	    * Plane strain 4-Noded Quadilateral Isoparametric elements with Linear Interpolation
	        * Full Integration currently 
	  * 3D Elements
	    * 3D 4-noded Tetrahedral  Isoparametric Elements with Linear Interpolation
	    * 3D 8--Noded Hexahedral Isoparametric Elements with Linear Interpolation
	      * Full Integration currently 
	  * Surface Elements
	    * 2D --- Cohesive Elements --- Representing a line at the midplane between 2 surfaces
	         * Xu-Needleman Cohesive Traction Separation Law
	2. Materials
	  * Isotropic Linear Elastic Material 	
  * Solvers 
	1. Linear Solver ---
	    Variety of solvers are available 
	    Currently MUMPS Linear Direct solver is implemented
	2. Non-linear Solver 
	    Uses SNES (petsc) to do solution currently uses default MUMPS to perform solution
  * Solution Contorls
	1. Incremental Solution from time t_init to t_final with prescribed dt
	2. Direct Solution without incrementation
Finite Element Implementation

Work in progress.
