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

  All my additions have a line of # on either side.

  When I test this with the test_trcopy.jl file it works and produces passive tracer advection in the top and bottom layers.

  Note: I think this is flexible enough so that you could use any amount of layers in the QG problem.
=#

module TracerAdvDiff_QG

export
   Problem,
   set_c!,
   QGset_c!,
   updatevars!,
   QGupdatevars!,
   vel_field_update!

using
  FFTW,
  Reexport

@reexport using FourierFlows

import LinearAlgebra: mul!, ldiv!
#################################################################################################
import GeophysicalFlows.MultiLayerQG
#################################################################################################

# --
# Problems
# --

"""
    Problem(; parameters...)
Construct a constant diffusivity problem with steady or time-varying flow.
"""
noflow(args...) = 0.0 # used as defaults for u, v functions in Problem()

function Problem(;
  #################################################################################################
  prob = nothing,
  #################################################################################################
    nx = 128,
    Lx = 2π,
    ny = nx,
    Ly = Lx,
  grid = TwoDGrid(nx, Lx, ny, Ly),
   kap = 0.1,
   eta = kap,
     u = noflow,
     v = noflow,
    dt = 0.01,
  stepper = "RK4",
  steadyflow = false, 
  nlayers = 1
  )
###################################################################################################
  #Added conditional if input is a GeophysicalFlows.Problem
  if isnothing(prob) == false
    nlayers = prob.params.nlayers
    nx = prob.grid.nx
    ny = nx
    Lx = prob.grid.Lx
    Ly = Lx
    grid = prob.grid
    dt = prob.clock.dt
    stepper = "FilteredRK4" #I had to hard code this
  end
###################################################################################################

  if steadyflow;                   pr = ConstDiffSteadyFlowParams(eta, kap, u, v, grid)
  #################################################################################################
  elseif isnothing(prob) == false; pr = set_QGFlowParams(eta, kap, prob)
  #################################################################################################
  else;                            pr = ConstDiffParams(eta, kap, u, v)
  end

  vs = Vars(grid, nlayers)
  eq = Equation(pr, grid)

  return FourierFlows.Problem(eq, stepper, dt, grid, vs, pr)
end


# --
# Params
# --

abstract type AbstractTracerParams <: AbstractParams end
abstract type AbstractConstDiffParams <: AbstractParams end
abstract type AbstractSteadyFlowParams <: AbstractParams end
#################################################################################################
abstract type AbstractQGFlowParams <: AbstractParams end
#################################################################################################

"""
    ConstDiffParams(eta, kap, kaph, nkaph, u, v)
    ConstDiffParams(eta, kap, u, v)
Returns the params for constant diffusivity problem with time-varying flow.
"""
struct ConstDiffParams{T} <: AbstractConstDiffParams
  eta::T                 # Constant isotropic horizontal diffusivity
  kap::T                 # Constant isotropic vertical diffusivity
  kaph::T                # Constant isotropic hyperdiffusivity
  nkaph::Int             # Constant isotropic hyperdiffusivity order
  u::Function            # Advecting x-velocity
  v::Function            # Advecting y-velocity
end

ConstDiffParams(eta, kap, u, v) = ConstDiffParams(eta, kap, 0eta, 0, u, v)

"""
    ConstDiffSteadyFlowParams(eta, kap, kaph, nkaph, u, v, g)
    ConstDiffSteadyFlowParams(eta, kap, u, v, g)
Returns the params for constant diffusivity problem with time-steady flow.
"""
struct ConstDiffSteadyFlowParams{T,A} <: AbstractSteadyFlowParams
  eta::T       # Constant horizontal diffusivity
  kap::T       # Constant vertical diffusivity
  kaph::T      # Constant isotropic hyperdiffusivity
  nkaph::Int   # Constant isotropic hyperdiffusivity order
  u::A         # Advecting x-velocity
  v::A         # Advecting y-velocity
end

function ConstDiffSteadyFlowParams(eta, kap, kaph, nkaph, u::Function, v::Function, g)
  x, y = gridpoints(g)
  ugrid = u.(x, y)
  vgrid = v.(x, y)
  
  return ConstDiffSteadyFlowParams(eta, kap, kaph, nkaph, ugrid, vgrid)
end

ConstDiffSteadyFlowParams(eta, kap, u, v, g) = ConstDiffSteadyFlowParams(eta, kap, 0eta, 0, u, v, g)

#################################################################################################
mutable struct QGFlowParams{T,A,Trfft} <: AbstractQGFlowParams
  nlayers::Int64  # Number of layers in the QG_flow problem
  eta::T          # Constant horizontal diffusivity
  kap::T          # Constant vertical diffusivity
  kaph::T         # Constant isotropic hyperdiffusivity
  nkaph::Int      # Constant isotropic hyperdiffusivity order
  u::A            # Advecting x-velocity. These will be updated at each step of the problem.
  v::A            # Advecting y-velocity
  rfftplan :: Trfft # rfft plan for MultiLayerQG
end

function set_QGFlowParams(eta, kap, prob)
  uvel = prob.vars.u
  vvel = prob.vars.v
  nlayers = prob.params.nlayers
  rfftplan = prob.params.rfftplan
  
  return QGFlowParams(nlayers, eta, kap, 0eta, 0, uvel, vvel, rfftplan) 
end
#################################################################################################

# --
# Equations
# --

"""
    Equation(p, g)
Returns the equation for constant diffusivity problem with params p and grid g.
"""
function Equation(p::ConstDiffParams, g::AbstractGrid{T}) where T
  L = zero(g.Krsq)
  @. L = -p.eta * g.kr^2 - p.kap * g.l^2 - p.kaph * g.Krsq^p.nkaph
  
  return FourierFlows.Equation(LC, calcN!, g)
end

function Equation(p::ConstDiffSteadyFlowParams, g::AbstractGrid{T}) where T
  L = zero(g.Krsq)
  @. L = -p.eta * g.kr^2 - p.kap * g.l^2 - p.kaph * g.Krsq^p.nkaph
  
  return FourierFlows.Equation(LC, calcN_steadyflow!, g)
end
#################################################################################################
function Equation(p::QGFlowParams, g::AbstractGrid{T}) where T
  nlayers = p.nlayers
  
  L = zeros(g.nkr, g.nl, nlayers)
  for n in 1:nlayers
    @. L[:, :, n] = -p.eta * g.kr^2 - p.kap * g.l^2 - p.kaph * g.Krsq^p.nkaph
  end
  
  return FourierFlows.Equation(L, calcN_QGflow!, g)
end
#################################################################################################

# --
# Vars
# --

# Construct Vars types
 physicalvars = [:c, :cx, :cy]
transformvars = [:ch, :cxh, :cyh]

struct Vars{Aphys,Atrans} <: AbstractVars
    c   :: Aphys
    cx  :: Aphys
    cy  :: Aphys
    ch  :: Atrans
    cxh :: Atrans
    cyh :: Atrans
end
"""
    Vars(g)
Returns the vars for constant diffusivity problem on grid g.
"""
function Vars(g, nlayers; T=typeof(g.Lx))
  ###################################################################################
  @createarrays T (g.nx, g.ny, nlayers) c cx cy
  @createarrays Complex{T} (g.nkr, g.nl, nlayers) ch cxh cyh
  #Have added the nlayers argument to these arrays for when using MultiLayerQG.Problem.
  ##################################################################################
  Vars(c, cx, cy, ch, cxh, cyh)
end

# --
# Solvers
# --

"""
    calcN!(N, sol, t, cl, v, p, g)
Calculate the advective terms for a tracer equation with constant diffusivity and time-varying flow.
"""
function calcN!(N, sol, t, cl, v, p::AbstractConstDiffParams, g)
  @. v.cxh = im*g.kr*sol
  @. v.cyh = im*g.l*sol

  ldiv!(v.cx, g.rfftplan, v.cxh) # destroys v.cxh when using fftw
  ldiv!(v.cy, g.rfftplan, v.cyh) # destroys v.cyh when using fftw

  x, y = gridpoints(g)
  @. v.cx = -p.u(x, y, cl.t) * v.cx - p.v(x, y, cl.t) * v.cy # copies over v.cx so v.cx = N in physical space
  
  mul!(N, g.rfftplan, v.cx)
  
  return nothing
end


"""
    calcN_steadyflow!(N, sol, t, cl, v, p, g)
Calculate the advective terms for a tracer equation with constant diffusivity and time-constant flow.
"""
function calcN_steadyflow!(N, sol, t, cl, v, p::AbstractSteadyFlowParams, g)
 
  @. v.cxh = im * g.kr * sol
  @. v.cyh = im * g.l  * sol


  ldiv!(v.cx, g.rfftplan, v.cxh) # destroys v.cxh when using fftw
  ldiv!(v.cy, g.rfftplan, v.cyh) # destroys v.cyh when using fftw

  @. v.cx = -p.u * v.cx - p. v* v.cy # copies over v.cx so v.cx = N in physical space
  
  mul!(N, g.rfftplan, v.cx)
  
  return nothing
end
#################################################################################################
"""
    calcN_QGflow!(N, sol, t, cl, v, p, g)
Calculate the advective terms for a tracer equation with constant diffusivity and flow from a QG problem.
"""
function calcN_QGflow!(N, sol, t, cl, v, p::AbstractQGFlowParams, g)
  @. v.cxh = im * g.kr * sol
  @. v.cyh = im * g.l  * sol

  MultiLayerQG.invtransform!(v.cx, v.cxh, p)
  MultiLayerQG.invtransform!(v.cy, v.cyh, p)

  @. v.cx = -p.u * v.cx - p.v * v.cy # copies over v.cx so v.cx = N in physical space
  
  MultiLayerQG.fwdtransform!(N, v.cx, p)
  
  return nothing
end
#################################################################################################

# --
# Helper functions
# --

"""
    updatevars!(prob)
Update the vars in v on the grid g with the solution in sol.
"""
function updatevars!(prob)
  v, g, sol = prob.vars, prob.grid, prob.sol
  v.ch .= sol
  ch1 = deepcopy(v.ch)
  ldiv!(v.c, g.rfftplan, ch1)
  
  return nothing
end

#################################################################################################
function QGupdatevars!(prob)
  v, g, sol = prob.vars, prob.grid, prob.sol
  v.ch .= sol
  ch1 = deepcopy(v.ch)
  MultiLayerQG.invtransform!(v.c, ch1, prob.params)
  
  return nothing
end
#################################################################################################

"""
    set_c!(prob, c)
    set_c!(prob, c::Function)
Set the solution sol as the transform of c and update variables v
on the grid g.
"""
function set_c!(prob, c)
  sol, v, g = prob.sol, prob.vars, prob.grid

  mul!(sol, g.rfftplan, c)
  updatevars!(prob)
  
  return nothing
end

function set_c!(prob, c::Function)
  sol, v, g = prob.sol, prob.vars, prob.grid

  x, y = gridpoints(g)
  
  cgrid = c.(x, y)
  
  mul!(sol, g.rfftplan, cgrid)
  updatevars!(prob)
  
  return nothing
end

#################################################################################################
function QGset_c!(prob, c)
  sol, v, g, nlayers = prob.sol, prob.vars, prob.grid, prob.params.nlayers
  
  C = zeros(g.nx, g.ny, nlayers)
  
  for n in 1:nlayers
    C[:, :, n] = c
  end
  
  MultiLayerQG.fwdtransform!(sol, C, prob.params)
  QGupdatevars!(prob)
  
  return nothing
end
#################################################################################################

#################################################################################################
function vel_field_update!(QG_prob,AD_prob)
  MultiLayerQG.stepforward!(QG_prob)
  MultiLayerQG.updatevars!(QG_prob)
  AD_prob.params.u = QG_prob.vars.u
  AD_prob.params.v = QG_prob.vars.v
  
  return nothing
end
#################################################################################################
end # module