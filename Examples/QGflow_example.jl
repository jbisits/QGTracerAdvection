#Passive tracer advection using a two layer QG flow (from the geophysical flows package).

using .TracerAdvDiff_QG
using .MeasureMixing

using GeophysicalFlows.MultiLayerQG, Plots, Distributions, JLD2, Random

#Set up the MultiLayerQG.Problem to advect-diffuse the tracer in.

#Choose CPU or GPU
dev = CPU()

#Numerical and time-stepping parameters
nx = 64        # 2D resolution = nx^2
ny = nx

stepper = "FilteredRK4";  # timestepper
Δt = 0.01                 # timestep
nsubs  = 1                # number of time-steps for plotting (nsteps must be multiple of nsubs)
nsteps = 2000nsubs        # total number of time-steps


#Physical parameters for a two layer QGproblem
Lx = 2π        # domain size
μ = 5e-2       # bottom drag
β = 5          # the y-gradient of planetary PV

nlayers = 2     # number of layers
f0, g = 1, 1    # Coriolis parameter and gravitational constant
H = [0.2, 0.8]  # the rest depths of each layer
ρ = [4.0, 5.0]  # the density of each layer

U = zeros(nlayers) # the imposed mean zonal flow in each layer
U[1] = 1.0
U[2] = 0.0

#Setup QGproblem and make easier to access the parts of the struct
QGprob = MultiLayerQG.Problem(nlayers, dev; nx=nx, Lx=Lx, f₀=f0, g=g, H=H, ρ=ρ, U=U, dt=Δt, stepper=stepper, μ=μ, β=β)
QGsol, QGclok, QGparams, QGvars, QGgrid = QGprob.sol, QGprob.clock, QGprob.params, QGprob.vars, QGprob.grid

#Set small random perturbation initial conditions for the PV field.

Random.seed!(1230) # reset of the random number generator for reproducibility
q₀  = 1e-2 * ArrayType(dev)(randn((QGgrid.nx, QGgrid.ny, nlayers)))
q₀h = QGprob.timestepper.filter .* rfft(q₀, (1, 2)) # only apply rfft in dims=1, 2
q₀  = irfft(q₀h, QGgrid.nx, (1, 2))                  # only apply irfft in dims=1, 2

MultiLayerQG.set_q!(QGprob, q₀)

#Set diffusivity
κ = 0.01
#Set delay time (that is flow for t seconds, then drop tracer in)
delay_time = Δt * 3000
#Set the tracer advection probelm by passing in the QG problem 
ADprob = TracerAdvDiff_QG.Problem(;prob = QGprob, delay_time = delay_time, inc_background_flow = true, nsubs = nsubs, κ = κ)
ADsol, ADclock, ADvars, ADparams, ADgrid = ADprob.sol, ADprob.clock, ADprob.vars, ADprob.params, ADprob.grid
ADxgrid, ADygrid = gridpoints(ADgrid)
x, y = ADgrid.x, ADgrid.y

#Set the (same) initial condition in both layers.
#A Gaussian blob centred at μIC 

μIC = [0, 0]
Σ = [1 0; 0 1]
blob = MvNormal(μIC, Σ)
blob_IC(x, y) = pdf(blob, [x, y])
C₀ = @. blob_IC(ADxgrid, ADygrid)

#A Gaussian strip around centred at μIC.
#=
μIC = 0
σ² = 0.5
strip = Normal(μIC, σ²)
strip_IC(x) = pdf(strip, x)
C₀ = Array{Float64}(undef, ADgrid.nx, ADgrid.ny)
for i in 1:ADgrid.nx
    C₀[i, :] = strip_IC(ADygrid[i, :])
end
=#
#If using strip_IC use C₀' for a vertical strip

#Set tracer initial condition in both layers
TracerAdvDiff_QG.QGset_c!(ADprob, C₀)

#Function to get the concentration field from the simulation, set frequecny at which data should be saved and create a file to save output to. 
function GetConcentration(ADprob)
    Concentration = @. ADprob.vars.c
    return Concentration
end
save_freq = 25
filename = joinpath(pwd(), "AdvectionDiffusionSim.jld2")
ADOutput = Output(ADprob, filename, (:Concentration, GetConcentration))
saveproblem(ADOutput)

#Step the tracer advection problem forward and plot at the desired time step.
while ADclock.step <= nsteps

    if ADclock.step % save_freq == 0
        saveoutput(ADOutput)
        println("Step number: ", round(Int, ADclock.step), ", saved data")
    end
    stepforward!(ADprob, nsubs)
    QGupdatevars!(ADprob)
    vel_field_update!(ADprob, QGprob, nsubs)

end

#Save the number of steps in the simulation and some other info
jldopen(ADOutput.path, "a+") do path
    path["clock/nsteps"] = ADclock.step - 1
    path["save_freq"] = save_freq
    if delay_time != 0
        path["delay_time"] = delay_time
    end
end

#Load in the saved simulation data
file = joinpath(pwd(), "AdvectionDiffusionSim.jld2")
data = load(file)

#Now the functions from the measure mixing module can be used by passing the saved data in.
tracerplots = tracer_plot(data; plot_freq = 500)
plot(tracerplots[:, 1]..., size = (1200, 1200))

concentration_variance = conc_var(data)
t = time_vec(data)
plot(t, concentration_variance, label = ["Upper layer" "Lower layer"])