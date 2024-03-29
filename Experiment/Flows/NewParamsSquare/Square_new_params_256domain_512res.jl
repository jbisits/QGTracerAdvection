#Flow on the non-dimensional domain of 256 (that is 256 * Ld in real space) and 512 resolution.

dev = CPU()
nx, ny = 512, 512
stepper = "FilteredRK4"

#Non dimensional paramters
Lx̂ = 256
Lŷ = 256

#Timestep
Δt̂ = 0.005

nlayers = 2
ρ̂ = [0.9, 1.0]
ĥ = 150^2
Ĥ = [ĥ/2, ĥ/2]
U₀ = 1.0
Û = [U₀, -U₀]

#Parameters
μ̂ = 0.65
ν̂ = 0.005
β̂ = 1
f̂₀ = 150
ĝ = 10

QGProb = MultiLayerQG.Problem(nlayers, dev; nx=nx, Lx=Lx̂, f₀=f̂₀, g=ĝ, H=Ĥ, ρ=ρ̂, U=Û, dt=Δt̂, stepper=stepper, μ=μ̂, β=β̂, ν=ν̂)
QGSol, QGClock, QGParams, QGVars, QGrid = QGProb.sol, QGProb.clock, QGProb.params, QGProb.vars, QGProb.grid

seed!(1230) # reset of the random number generator for reproducibility
q₀  = 1e-2 * device_array(dev)(randn((QGrid.nx, QGrid.ny, nlayers)))
q₀h = QGProb.timestepper.filter .* rfft(q₀, (1, 2)) # only apply rfft in dims=1, 2
q₀  = irfft(q₀h, QGrid.nx, (1, 2)) # only apply irfft in dims=1, 2

MultiLayerQG.set_q!(QGProb, q₀)