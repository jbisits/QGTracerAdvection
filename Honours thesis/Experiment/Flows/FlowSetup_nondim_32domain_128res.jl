#Flow on the non-dimensional domain of 32 (that is 32 * Ld in real space) and 128 resolution.

dev = CPU()
nx, ny = 128, 128
stepper = "FilteredRK4"
nsubs  = 1         
nsteps = 8000nsubs

#Non dimensional paramters
Lx̂ = 32
Lŷ = 32

#Timestep
Δt̂ = 0.01

nlayers = 2
ρ̂ = [0.9, 1.0]
Ĥ = [0.5, 0.5]
U₀ = 1.0
Û = [U₀, -U₀]

#Parameters
μ̂ = 0.1
ν̂ = 0.2
β̂ = 1
f̂₀ = 1
ĝ = 10

QGProb = MultiLayerQG.Problem(nlayers, dev; nx=nx, Lx=Lx̂, f₀=f̂₀, g=ĝ, H=Ĥ, ρ=ρ̂, U=Û, dt=Δt̂, stepper=stepper, μ=μ̂, β=β̂, ν=ν̂)
QGSol, QGClock, QGParams, QGVars, QGrid = QGProb.sol, QGProb.clock, QGProb.params, QGProb.vars, QGProb.grid

seed!(1234) # reset of the random number generator for reproducibility
q₀  = 1e-2 * ArrayType(dev)(randn((QGrid.nx, QGrid.ny, nlayers)))
q₀h = QGProb.timestepper.filter .* rfft(q₀, (1, 2)) # only apply rfft in dims=1, 2
q₀  = irfft(q₀h, QGrid.nx, (1, 2)) # only apply irfft in dims=1, 2

MultiLayerQG.set_q!(QGProb, q₀)