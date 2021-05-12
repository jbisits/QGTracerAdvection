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
Construct a constant diffusivity problem with steady, time-varying or quasi-geostrophic flow.
"""
noflow(args...) = 0.0 # used as defaults for u, v functions in Problem()

function Problem(;
  #################################################################################################
  prob = nothing,
  delay_time = 0,
  nsubs = 1,
  #################################################################################################
    nx = 128,
    Lx = 2π,
    ny = nx,
    Ly = Lx,
  grid = TwoDGrid(nx, Lx, ny, Ly),
   κ = 0.1,
   eta = κ,
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
    stepper = "FilteredRK4"
  end
###################################################################################################

  if steadyflow;                   pr = ConstDiffSteadyFlowParams(eta, κ, u, v, grid)
  #################################################################################################
  elseif isnothing(prob) == false; pr = set_QGFlowParams(eta, κ, prob, delay_time, nsubs)
  #################################################################################################
  else;                            pr = ConstDiffParams(eta, κ, u, v)
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
    ConstDiffParams(eta, κ, κh, nκh, u, v)
    ConstDiffParams(eta, κ, u, v)
Returns the params for constant diffusivity problem with time-varying flow.
"""
struct ConstDiffParams{T} <: AbstractConstDiffParams
  eta::T                 # Constant isotropic horizontal diffusivity
  κ::T                 # Constant isotropic vertical diffusivity
  κh::T                # Constant isotropic hyperdiffusivity
  nκh::Int             # Constant isotropic hyperdiffusivity order
  u::Function            # Advecting x-velocity
  v::Function            # Advecting y-velocity
end

ConstDiffParams(eta, κ, u, v) = ConstDiffParams(eta, κ, 0eta, 0, u, v)

"""
    ConstDiffSteadyFlowParams(eta, κ, κh, nκh, u, v, g)
    ConstDiffSteadyFlowParams(eta, κ, u, v, g)
Returns the params for constant diffusivity problem with time-steady flow.
"""
struct ConstDiffSteadyFlowParams{T,A} <: AbstractSteadyFlowParams
  eta::T       # Constant horizontal diffusivity
  κ::T       # Constant vertical diffusivity
  κh::T      # Constant isotropic hyperdiffusivity
  nκh::Int   # Constant isotropic hyperdiffusivity order
  u::A         # Advecting x-velocity
  v::A         # Advecting y-velocity
end

function ConstDiffSteadyFlowParams(eta, κ, κh, nκh, u::Function, v::Function, g)
  x, y = gridpoints(g)
  ugrid = u.(x, y)
  vgrid = v.(x, y)
  
  return ConstDiffSteadyFlowParams(eta, κ, κh, nκh, ugrid, vgrid)
end

ConstDiffSteadyFlowParams(eta, κ, u, v, g) = ConstDiffSteadyFlowParams(eta, κ, 0eta, 0, u, v, g)

#################################################################################################
"""
    QGFlowParams(eta, κ, κh, nκh, u, v, g)
    QGFlowParams(eta, κ, u, v, g)
Returns the params for constant diffusivity problem with quasi-geostrophic flow.
"""
struct QGFlowParams{T,A,Trfft} <: AbstractQGFlowParams
  nlayers::Int64    # Number of layers in the QG_flow problem
  eta::T            # Constant horizontal diffusivity
  κ::T            # Constant vertical diffusivity
  κh::T           # Constant isotropic hyperdiffusivity
  nκh::Int        # Constant isotropic hyperdiffusivity order
  u::A              # Advecting x-velocity. These will be updated at each step of the problem.
  v::A              # Advecting y-velocity
  rfftplan :: Trfft # rfft plan for MultiLayerQG
end

function set_QGFlowParams(eta, κ, prob, delay_time, nsubs)
  if isequal(delay_time, 0) == false
    flow_till_delay_time!(prob, delay_time, nsubs)
  end
  uvel = prob.vars.u
  vvel = prob.vars.v
  nlayers = prob.params.nlayers
  rfftplan = prob.params.rfftplan
  
  return QGFlowParams(nlayers, eta, κ, 0eta, 0, uvel, vvel, rfftplan) 
end
#################################################################################################

# --
# Equations
# --

"""
    Equation(p, g)
Returns the equation for constant diffusivity problem with params p and grid g.
"""
function Equation(p::ConstDiffParams, grid)
  L = zero(grid.Krsq)
  @. L = -p.eta * grid.kr^2 - p.κ * grid.l^2 - p.κh * grid.Krsq^p.nκh
  
  return FourierFlows.Equation(L, calcN!, grid)
end

function Equation(p::ConstDiffSteadyFlowParams, grid)
  L = zero(grid.Krsq)
  @. L = -p.eta * grid.kr^2 - p.κ * grid.l^2 - p.κh * grid.Krsq^p.nκh
  
  return FourierFlows.Equation(L, calcN_steadyflow!, grid)
end
#################################################################################################
function Equation(p::QGFlowParams, grid)
  nlayers = p.nlayers
  
  L = zeros(grid.nkr, grid.nl, nlayers)
  for n in 1:nlayers
    @. L[:, :, n] = -p.eta * grid.kr^2 - p.κ * grid.l^2 - p.κh * grid.Krsq^p.nκh
  end
  
  return FourierFlows.Equation(L, calcN_QGflow!, grid)
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
  if nlayers > 1
    @createarrays T (g.nx, g.ny, nlayers) c cx cy
    @createarrays Complex{T} (g.nkr, g.nl, nlayers) ch cxh cyh
  else
    @createarrays T (g.nx, g.ny) c cx cy
    @createarrays Complex{T} (g.nkr, g.nl) ch cxh cyh
  end
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
"""
  function QGset_c!(prob, c)
Sets the initial condition in the QG tracer advection `prob` to be the array `c`
"""
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
"""
  function vel_field_update!(AD_prob, QG_prob, nsubs)
Updates the velocity field of the advection problem `AD_prob` with the velcoty field from the `QG_prob`.
"""
function vel_field_update!(AD_prob, QG_prob, nsubs)
  MultiLayerQG.stepforward!(QG_prob, nsubs)
  MultiLayerQG.updatevars!(QG_prob)
  @. AD_prob.params.u = QG_prob.vars.u + QG_prob.params.U
  @. AD_prob.params.v = QG_prob.vars.v
  return nothing
end
#################################################################################################

#################################################################################################
"""
  function
Flows the QG problem (`QG_prob`) forward until `delay_time` so that the tracer can be dropped into 
the QG flow at a time specified by the user.
"""
function flow_till_delay_time!(QG_prob, delay_time, nsubs)

  while QG_prob.clock.t < delay_time
    MultiLayerQG.stepforward!(QG_prob, nsubs)
    MultiLayerQG.updatevars!(QG_prob)
  end

  return nothing
end
#################################################################################################
end # module
