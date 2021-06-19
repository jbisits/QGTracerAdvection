#Set up a flow with parameters that have mid-latitude values on a domain of 1500km * 1500km

#Mid-latitude values (all in metres and seconds)
Ω = 7.29e-5     # Earth"s rotation
ϕ = π/3         # Latitude
a = 6378e3      # Earth's radius
g = 9.81        # Gravity
H = 1000        # Total depth (in metres)
ρ₁ = 1034       # Density of top layer
ρ₂ = 1035       # Density of bottom layer

β = (2*Ω*cos(ϕ))/a          # β value computed from above values
f₀ = 2*Ω*sin(ϕ)             # Coriolis force
gprime = g*((ρ₂ - ρ₁)/ρ₂)   # Reduced gravity

Ld = sqrt(gprime*H)/(f₀) #Rossby deformation radius (horizontal domain to be 50-100 times larger)

using Random: seed!

dev = CPU() 

      n = 128               # 2D resolution = n²
stepper = "FilteredRK4"     # timestepper
     dt = 3600              # timestep (one day in seconds)
 nsteps = 365*15            # total number of time-steps (some number of years)
 nsubs  = 1                 # number of time-steps for plotting (nsteps must be multiple of nsubs)



# ## Physical parameters
 L = 1.5e6                  # domain size (in metres)
 μ = 1e-4                   # bottom drag 
 ν = 1e-3                   # viscousity
 
nlayers = 2                 # number of layers 
 H = [H/2, H/2]             # the rest depths of each layer
 ρ = [ρ₁, ρ₂]               # the density of each layer
 
 U = zeros(nlayers)         # the imposed mean zonal flow in each layer to maintain a constant shear
 U[1] = 0.21
 U[2] = -U[1]

#Set the QG problem and shortcuts
QGProb = MultiLayerQG.Problem(nlayers, dev; nx=n, Lx=L, f₀=f₀, g=g, H=H, ρ=ρ, U=U, dt=dt, stepper=stepper, μ=μ, ν=ν, β=β)
QGSol, QGClock, QGParams, QGVars, QGGrid = QGProb.sol, QGProb.clock, QGProb.params, QGProb.vars, QGProb.grid
QGx, QGy = QGGrid.x, QGGrid.y

#Set initial condition of small amplitude random noise that is reproducible by setting the seed
seed!(1234)                 
q₀  = 1e-5 * ArrayType(dev)(randn((QGGrid.nx, QGGrid.ny, nlayers)))
q₀h = QGProb.timestepper.filter .* rfft(q₀, (1, 2))                     # only apply rfft in dims=1, 2
q₀  = irfft(q₀h, QGGrid.nx, (1, 2))                                     # only apply irfft in dims=1, 2
MultiLayerQG.set_q!(QGProb, q₀)