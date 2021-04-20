# PerfDoms

This code computes the homogenized matrix for the problem  
<img src="https://render.githubusercontent.com/render/math?math=-\Delta u_\epsilon = f \, \text{ in }\, \Omega_\epsilon,">  
<img src="https://render.githubusercontent.com/render/math?math=\partial_\nu u_\epsilon = 0 \, \text{ on } \, \partial\Omega_\epsilon,">  
where <img src="https://render.githubusercontent.com/render/math?math=\Omega_\epsilon"> is a smooth domain in <img src="https://render.githubusercontent.com/render/math?math=\mathbb{R}^2"> perforated periodically by circular holes whose size and distance is proportional to <img src="https://render.githubusercontent.com/render/math?math=\epsilon">.

## Theory
As <img src="https://render.githubusercontent.com/render/math?math=\epsilon\to0"> the periodic fine structure gets "averaged out" and the solutions <img src="https://render.githubusercontent.com/render/math?math=u_\epsilon"> converge to a function <img src="https://render.githubusercontent.com/render/math?math=u"> which solves  
<img src="https://render.githubusercontent.com/render/math?math=-\nabla\cdot(A_{\text{hom}}\nabla u) = f \, \text{ in }\, \Omega,">  
<img src="https://render.githubusercontent.com/render/math?math=\partial_\nu u_\epsilon = 0 \, \text{ on } \, \partial\Omega,">  
where the so-called *homogenized coefficient matrix* <img src="https://render.githubusercontent.com/render/math?math=A_{\text{hom}}"> can be computed by solving the following problem on the unit cell (with a single copy F of the holes whose size is of order 1). For <img src="https://render.githubusercontent.com/render/math?math=\xi\in\mathbb{R}^2">, solve  
<img src="https://render.githubusercontent.com/render/math?math=\Delta N = 0 \, \text{ in }\, [0,1]^2\setminus F,">  
<img src="https://render.githubusercontent.com/render/math?math=\partial_\nu N = -\nu\cdot\xi \, \text{ on } \, \partial F,">  
<img src="https://render.githubusercontent.com/render/math?math=N \,\text{ periodic}">.  
Then <img src="https://render.githubusercontent.com/render/math?math=A_{\text{hom}}"> is given by the mean value  
<img src="https://render.githubusercontent.com/render/math?math=A_{\text{hom}}\xi = \langle\xi %2B \nabla N\rangle_{[0,1]^2\setminus F}">.

Alternatively <img src="https://render.githubusercontent.com/render/math?math=A_{\text{hom}}"> can be computed via the following variational characterization:  
<img src="https://render.githubusercontent.com/render/math?math=\xi^\top A_{\text{hom}}\xi = \inf\limits_{u\in H^1_{\text{per}}([0,1]^2)}\langle |\xi %2B \nabla u|^2\rangle_{[0,1]^2\setminus F}">.

## Implementation
The repository contains Matlab scripts that compute the homogenized matrix using either of the two methods above. The hole F is taken to be circular, but can easily be adapted to other shapes. In the first case a finite element method is used to solve the cell problem. A triangulation for a hole of radius `r=0.2` is included in the repository. In order to generate new meshes, the [`distmesh`](http://persson.berkeley.edu/distmesh/) package is required.  
**Functions:**  
Cell problem:
* `generate_mesh(r,h)` - generates mesh of size h for a circular hole of radius r (uses `distmesh`)
* `cell_problem(r,h,xi,plotting)` - calles `generate_mesh(r,h)` and then solves the cell problem for the vector xi and computes <img src="https://render.githubusercontent.com/render/math?math=A_{\text{hom}}\xi">. If `plotting==true`, the solution N is plotted.
* 
Variational problem:
* `variational_problem(r,N,xi,plotting)` - solves variational problem via finite difference method and returns the minimum. 


Any comments or queries are welcome at https://frank-roesler.github.io/contact/
