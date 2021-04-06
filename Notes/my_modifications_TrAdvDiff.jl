#=
  My (Joey) modifications to the TracerAdvDiff module to hopefully get it to work with geophysical flows:
  1. Set the problem function to take in a FourieFlows.Problem (the default setting for the problem is false).
  2. Set the the tracer advection problem with the values from the input problem.
  3. Create a new mutable struct that stores the parameters for TracerAdvectionDiffusion.Problem.
     The struct is mutable so that values of u and v (the velocity field) can be updated.
  4. Write a function to set the values in the QGFlowParams struct.
  5. Added a function vel_field_update to step forward the MultiLayerQG.Problem and update the struct values u and v to the new values.
  6. Modified the Vars struct and function to set the vars so that it sets the "C" field to have the same number of layers as the
     MultiLayerQG.Problem. It defualts to one so should not change the way it currently works.
  7. Wrote new Equation and calcN_QGflow! fucntions that also make this adjustment to allow for flow advected in multiple layers at the same time.
  8. Wrote a QG version of updatevars! and set_c!. They needed to be able to deal with the multidimensional arrays that the two layer probem has.
  9. Add function to let user drop the tracer into a flow at time they choose. This is done in two ways and not sure which is best.
     Will ask the others to take a look at it.

  All my additions have a line of # on either side.

  When I test this with the test_trcopy.jl file it works and produces passive tracer advection in the top and bottom layers.

  Note: I think this is flexible enough so that you could use any amount of layers in the QG problem.
=#
