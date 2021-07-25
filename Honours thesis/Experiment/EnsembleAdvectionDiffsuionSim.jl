#Ensemble tracer advection diffusion experiment

#Change to the correct directory (if it was not already correct for some reason)
SimPath = joinpath("/Users/Joey/Desktop/ThesisCode/QG_tracer_advection", "Honours thesis/Experiment")
cd(SimPath)
#Load in all the required packages for the simulation
include("PackageSetup.jl")

#Import a flow that has already been set up from the Flows folder
#include("Flows/FlowSetup_nondim_32domain_64res.jl")
#include("Flows/FlowSetup_nondim_64domain_128res.jl")
include("Flows/FlowSetup_nondim_128domain_256res.jl")
#include("Flows/FlowSetup_nondim_256domain_512res.jl")

nsubs  = 1           #Set the number of steps the simulation takes at each iteration. This is also the frequency that data is saved at.         
nsteps = 5000           #Set the total amount of time steps the advection-diffusion simulation should run for

κ = 0.01
#Set delay times (that is flow for some length of time, then drop tracer in)
delay_time = Δt̂ * 3000
#Define number of tracer advection simulations
ADSims = 10
#Set the frequency at which to save data
save_freq = 50
#Initial condition
μIC = 0
σ² = 1
IC = GaussianStripIC(μIC, σ², ADGrid)
#File name for saving, FourierFlows creates a new file each time with _i appended
filename = CreateFile(ADProb, IC, save_freq, SimPath; Ensemble = true)

#This runs a non-parallel simulation where an array of advection-diffusion problems is defined then stepped forward separately with the flow reset each time
for i ∈ 1:ADSims

    if i == 1
        ADProb = TracerAdvDiff_QG.Problem(;prob = QGProb, delay_time = delay_time, nsubs = nsubs, κ = κ)
        ADSol, ADClock, ADVars, ADParams, ADGrid = ADProb.sol, ADProb.clock, ADProb.vars, ADProb.params, ADProb.grid
    else
        #Reset the QG flow
        global QGProb = MultiLayerQG.Problem(nlayers, dev; nx=nx, Lx=Lx̂, f₀=f̂₀, g=ĝ, H=Ĥ, ρ=ρ̂, U=Û, dt=Δt̂, stepper=stepper, μ=μ̂, β=β̂, ν=ν̂)
        global QGSol, QGClock, QGParams, QGVars, QGrid = QGProb.sol, QGProb.clock, QGProb.params, QGProb.vars, QGProb.grid
        seed!( parse(Int64, join([1, 2, 3, i])) ) # reset of the random number generator for reproducibility
        local q₀  = 1e-2 * ArrayType(dev)(randn((QGrid.nx, QGrid.ny, nlayers)))
        local q₀h = QGProb.timestepper.filter .* rfft(q₀, (1, 2)) # only apply rfft in dims=1, 2
        q₀  = irfft(q₀h, QGrid.nx, (1, 2)) # only apply irfft in dims=1, 2
        MultiLayerQG.set_q!(QGProb, q₀)
        ADProb = TracerAdvDiff_QG.Problem(;prob = QGProb, delay_time = delay_time, nsubs = nsubs, κ = κ)
        ADSol, ADClock, ADVars, ADParams, ADGrid = ADProb.sol, ADProb.clock, ADProb.vars, ADProb.params, ADProb.grid
    end

    QGset_c!(ADProb, IC.C₀)
    ADOutput = Output(ADProb, filename, (:Concentration, GetConcentration))
    saveproblem(ADOutput)
    
    #Simulation loop
    while ADClock.step <= nsteps

        if ADClock.step % save_freq == 0
            saveoutput(ADOutput)
            println("Step number: ", round(Int, ADClock.step), ", saved data")
        end
        stepforward!(ADProb, nsubs)
        QGupdatevars!(ADProb)
        vel_field_update!(ADProb, QGProb, nsubs)

    end

    #Save the number of steps and other info from the simulation
    jldopen(ADOutput.path, "a+") do path
        path["clock/nsteps"] = ADClock.step - 1
        path["save_freq"] = save_freq
        if delay_time != 0
            path["delay_time"] = delay_time
        end
    end

end