#Flow on the non-dimensional domain of 32 (that is 32 * Ld in real space) and 256 resolution.

dev = CPU()
nx, ny = 64, 64
stepper = "FilteredRK4"

#Non dimensional paramters
Lx̂ = 32
Lŷ = 32

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

QGProbs = Array{FourierFlows.Problem}(undef, ADSims)
seed!(1234) # reset of the random number generator for reproducibility, but change if want different intial conditions
for i ∈ 1:ADSims
    QGProbs[i] = MultiLayerQG.Problem(nlayers, dev; nx=nx, Lx=Lx̂, f₀=f̂₀, g=ĝ, H=Ĥ, ρ=ρ̂, U=Û, dt=Δt̂, stepper=stepper, μ=μ̂, β=β̂, ν=ν̂)
    QGSol, QGClock, QGParams, QGVars, QGrid = QGProbs[i].sol, QGProbs[i].clock, QGProbs[i].params, QGProbs[i].vars, QGProbs[i].grid

    q₀  = 1e-2 * ArrayType(dev)(randn((QGrid.nx, QGrid.ny, nlayers)))
    q₀h = QGProbs[i].timestepper.filter .* rfft(q₀, (1, 2)) # only apply rfft in dims=1, 2
    q₀  = irfft(q₀h, QGrid.nx, (1, 2)) # only apply irfft in dims=1, 2

    MultiLayerQG.set_q!(QGProbs[i], q₀)
end