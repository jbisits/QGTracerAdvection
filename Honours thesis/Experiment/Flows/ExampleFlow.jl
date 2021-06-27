#Example flow from GeophysicalFlows package

dev = CPU()

#Numerical and time-stepping parameters
nx = 64        # 2D resolution = nx^2

stepper = "FilteredRK4";  # timestepper
Δt = 0.01                 # timestep
nsubs  = 1                # number of time-steps for plotting (nsteps must be multiple of nsubs)
nsteps = 6000nsubs        # total number of time-steps


#Physical parameters for a two layer QG_problem
Lx = 2π        # domain size
μ = 5e-2       # bottom drag
β = 5          # the y-gradient of planetary PV

nlayers = 2     # number of layers
f₀, g = 1, 1    # Coriolis parameter and gravitational constant
H = [0.2, 0.8]  # the rest depths of each layer
ρ = [4.0, 5.0]  # the density of each layer

U = zeros(nlayers) # the imposed mean zonal flow in each layer
U[1] = 1.0
U[2] = 0.0

#Setup QG_problem and make easier to access the parts of the struct
QGProb = MultiLayerQG.Problem(nlayers, dev; nx=nx, Lx=Lx, f₀=f₀, g=g, H=H, ρ=ρ, U=U, dt=Δt, stepper=stepper, μ=μ, β=β)
QGSol, QGClock, QGParams, QGVars, QGGrid = QGProb.sol, QGProb.clock, QGProb.params, QGProb.vars, QGProb.grid
QGx, QGy = QGGrid.x, QGGrid.y

#Set initial conditions.=
seed!(1234) # reset of the random number generator for reproducibility
q₀  = 1e-2 * ArrayType(dev)(randn((QGGrid.nx, QGGrid.ny, nlayers)))
q₀h = QGProb.timestepper.filter .* rfft(q₀, (1, 2)) # apply rfft  only in dims=1, 2
q₀  = irfft(q₀h, QGGrid.nx, (1, 2))                 # apply irfft only in dims=1, 2

MultiLayerQG.set_q!(QGProb, q₀)