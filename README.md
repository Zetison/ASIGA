# ASIGA
The ASIGA (Acoustic Scattering using IsoGeometric Analysis) toolbox provides a framework for simulating acoustic scattering problems using IGA. 


The main emphesis is acoustic scattering using plane wave (applyLoad = 'planeWave') in which the following cases are implemented
- bistatic scattering (scatteringCase = 'BI')
- monostatic scattering (scatteringCase = 'MS')
- frequency sweep (scatteringCase = 'Sweep')


The following methods has been implemented (with availabel formulations)
- IGA using the infinite element method (method = 'IE')
	- The Bubnov--Galerkin Conjugated formulation (formulation = 'BGC')
	- The Petrov--Galerkin Conjugated formulation (formulation = 'PGC')
	- The Bubnov--Galerkin Unconjugated formulation (formulation = 'BGU')
	- The Petrov--Galerkin Unconjugated formulation (formulation = 'PGU')
- IGA using the IE method after Shirron (method = 'IENSG')
	- The Bubnov--Galerkin Conjugated formulation (formulation = 'BGC')
	- The Petrov--Galerkin Conjugated formulation (formulation = 'PGC')
	- The Bubnov--Galerkin Unconjugated formulation (formulation = 'BGU')
	- The Petrov--Galerkin Unconjugated formulation (formulation = 'PGU')
- IGA using absorbing boundary conditions (method = 'ABC')
	- Bayliss-GunzBurger-Turkel-operators (formulation = 'BGT')
- IGA using best approximation (method = 'BA')
	- Best approximation in L^2(\Gamma) (formulation = 'SL2E')
	- Best approximation in L^2(\Omega_a) (formulation = 'VL2E')
- IGA using the boundary element methdo (method = 'BEM')
	- Conventional boundary integral equation using collocation method (formulation = 'CCBIE')
	- Hypersingular boundary integral equation using collocation method (formulation = 'CHBIE')
	- Buron-Miller formulation using collocation method (formulation = 'CBM')
	- Three regularized conventional boundary integral equations using collocation method (formulation = 'CCBIE1', formulation = 'CCBIE2' and formulation = 'CCBIE3')
   The Gallerkin method may be used instead of the collocation method by replacing the prepended letter 'C' with 'G' (i.e. formulation = 'GCBIE'). Moreover, by appending the letter 'C' the CHIEF method is applied in an attempt to remove fictitious eigenfrequencies.
- IGA using Kirchhoff approximations (method = 'KDT')
	- A first order multiple bounce implementation (formulation = 'MS1')
	- (Not completed!) A second order multiple bounce implementation (formulation = 'MS2')
- The method of fundamental solutions (method = 'MFS')
	- Using fundamental solutions (Phi_k) as basis functions (formulation = 'PS')
	- Using spherical harmonics as basis functions (formulation = 'SS')
- IGA using Beam tracing (method = 'RT')


Instead of IGA (coreMethod = 'IGA') the following alternatives are implemented
- Bilinear isoparametric FEM (coreMethod = 'linear_FEM')
- Subparametric FEM with bilinear representation of geometry (coremethod = 'h_FEM')
- Isoparametric FEM (coreMethod = 'hp_FEM')


Disclaimer: No efforts has as of now been made to make the program user friendly, and certain combinations of parameters/methods are not implemented or may be erroneous. Please report any bug to jonvegard89@yahoo.com
