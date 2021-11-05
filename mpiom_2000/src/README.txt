@mainpage The Max-Planck-Institute Global Ocean/Sea-Ice Model

<center>
    Patrick Wetzel, Helmuth Haak, Johann Jungclaus and Ernst Maier-Reimer<br>
    Max-Planck-Institute for Meteorology
</center>

@section sec_intro Introduction

The Hamburg Ocean Model, MPI-OM is the successor of the 
Hamburg Ocean Primitive Equation (HOPE) model and has undergone significant
development in recent years.
Most notable is the treatment of horizontal discretization
which has undergone transition from a staggered E-grid
to an orthogonal curvilinear C-grid.
The treatment of subgridscale mixing has been improved by the
inclusion of a new formulation of bottom boundary layer (BBL) slope convection,
an isopycnal diffusion scheme,
and a Gent and McWilliams style eddy-induced mixing parameterization.


@subsection sub_details Details

The  Hamburg Ocean Model, MPI-OM, is an ocean general circulation model (OGCM) 
based on the primitive equations with representation of thermodynamic 
processes. It is capable of simulating the oceanic circulation from 
small scales (oceanic eddies) to gyre scales, in response to atmospheric 
forcing fields. For an application on horizontal scales smaller than 
about 1 km the hydrostatic assumption is no longer valid and the model 
must in parts be reformulated. The use of an ocean circulation model 
requires a comprehensive understanding of the ocean physics and the 
numerical formulation. It is not recommended to use this model as a 
black box. Many physical processes in the ocean are still not very 
well understood and are therefore only crudely parameterized. 
Each new application requires a new consideration of how to specify 
model parameters, especially coefficients for eddy viscosity and 
diffusivity. The model is thought to be a framework into which 
new ideas concerning parameterizations or forcing mechanisms might 
easily be incorporated. This manual gives a description of the MPI-OM 
model, in order to help potential users to run the model and to acquire a 
understanding of the model physics and numerics. 


@section sec_apply Applicability

This manual refers to release_1_1 of MPI-OM.

@defgroup common Global settings and Variables
