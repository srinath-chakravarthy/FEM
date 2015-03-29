# FEM
Parallel FEM for linear and non-linear small deformation static problems using PETSC 
PETSC --- Parallel Extensible Toolkit for Scientific Computation (http://www.mcs.anl.gov/petsc/)
  Requirements
	1) Requires PETSC installation with MUMPS, METIS and PARMETIS
  Input Parameters
	1) Mesh input file along with Boundary Conditions need to be given 
	  === Description of Mesh file
	TODO ---- Read Comment lines 
	TODO ---- Read specific command lines
	TODO ---- Give Element Specific names 2DQPSLF--> 2D--Quad--Plane Strain--Linear--Full Integration
	
  Output 
	vtk files from each processor numbered by rank
	Read using Paraview
	 1) Displacement Vectors at each nodal point
	 TODO ---- a) Output stress at integration points
	 TODO ---- b) Output strain at integration points
	 TODO ---- c) Output custom quantities at integration points
	 2) Stress Tensor, averaged at each node, represented using Scalars (1..6)
  Implemetation
	1) Elements
	  2D Elements 
	    a) Plane strain 3-Noded Triangular Isoparametric elements with Linear Interpolation 
	    b) Plane strain 4-Noded Quadilateral Isoparametric elements with Linear Interpolation
	      (i) Full Integration currently 
	    TODO --- Plane Stress and Axisymmetric
	    TODO ----  (ii) Reduced Integration
	    TODO ---- 6 Noded Triangle and 8 noded Serendipity Quadilateral
	      TODO ----- Reduced integration
	  3D Elements
	    a) 3D 4-noded Tetrahedral  Isoparametric Elements with Linear Interpolation
	    TODO ---- Validation
	    b) 3D 8--Noded Hexahedral Isoparametric Elements with Linear Interpolation
	    TODO --- Validation
	      (i) Full Integration currently 
	      TODO ---- (ii) Reduced Integration
	  Surface Elements
	    a) 2D --- Cohesive Elements --- Representing a line at the midplane between 2 surfaces
	      (i) Xu-Needleman Cohesive Traction Separation Law
	    TODO --- Other Traction Separtion laws
	2) Materials
	  a) Isotropic Linear Elastic Material 
	    TODO ---- (i) Anisotropic Linear Elastic Material
	      TODO ---- a) Small Strain kinematic and Isotropic J2 Plasticity 
	      TODO ---- b) Small strain incremental viscoplasticity
	      TODO ---- c) Small strain incremental non-linear elasticity
	      TODO ---- d) Stress-Gradient Plasticity (non-local ---> averaging)
  Solvers 
	1) Linear Solver ---
	    Variety of solvers are available 
	    Currently MUMPS Linear Direct solver is implemented
	    TODO ---- Give user access to different solvers
	2) Non-linear Solver 
	    Uses SNES (petsc) to do solution
	      currently uses default MUMPS to perform solution
	    TODO --- Give user access to Solver methods
  Solution Contorls
	1) Incremental Solution from time t_init to t_final with prescribed dt
	  TODO --- Since all problems are static, normalize total time to 1 and dt is a fraction of 1
	  TODO --- Variable dt for non-linear problems
	2) Direct Solution without incrementation
Finite Element Implementation

Work in progress.
