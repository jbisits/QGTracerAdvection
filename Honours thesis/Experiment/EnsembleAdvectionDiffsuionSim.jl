#Ensemble tracer advection diffusion experiment

#Change to the correct directory (if it was not already correct for some reason)
SimPath = joinpath("/Users/Joey/Desktop/ThesisCode/QG_tracer_advection", "Honours thesis/Experiment")
cd(SimPath)
#Load in all the required packages for the simulation
include("PackageSetup.jl")

#Import a flow that has already been set up from the Flows folder
#include("Flows/ExampleFlow.jl")
include("Flows/FlowSetup_nondim_32domain_128res.jl")

nsubs  = 200            #Set the number of steps the simulation takes at each iteration. This is also the frequency that data is saved at.         
nsteps = 5000           #Set the total amount of time steps the advection-diffusion simulation should run for

κ = 0.01
#Set delay time (that is flow for some length of time, then drop tracer in)
delay_time = Δt̂ * 3000
#Define number of tracer advection simulations
ADSims = 2

#This runs a non-parallel simulation where an array of advection-diffusion problems is defined then stepped forward separately with the flow reset each time
for i ∈ 1:ADSims

    if i == 1
        ADProb = TracerAdvDiff_QG.Problem(;prob = QGProb, delay_time = delay_time, nsubs = nsubs, κ = κ)
        ADSol, ADClock, ADVars, ADParams, ADGrid = ADProb.sol, ADProb.clock, ADProb.vars, ADProb.params, ADProb.grid
        #Set the input for Gaussian blob initial condition
        μIC = [0, 0]
        Σ = [1 0; 0 1]
        global IC = GaussianBlobIC(μIC, Σ, ADGrid)
        #File name for saving. Only needs to be done the first time in the loop.
        global filename = CreateFile(ADProb, IC, nsubs, SimPath; Ensemble = true)
    else
        #Reset the QG flow
        global QGProb = MultiLayerQG.Problem(nlayers, dev; nx=nx, Lx=Lx̂, f₀=f̂₀, g=ĝ, H=Ĥ, ρ=ρ̂, U=Û, dt=Δt̂, stepper=stepper, μ=μ̂, β=β̂, ν=ν̂)
        global QGSol, QGClock, QGParams, QGVars, QGrid = QGProb.sol, QGProb.clock, QGProb.params, QGProb.vars, QGProb.grid
        MultiLayerQG.set_q!(QGProb, q₀)
        ADProb = TracerAdvDiff_QG.Problem(;prob = QGProb, delay_time = delay_time, nsubs = nsubs, κ = κ)
        ADSol, ADClock, ADVars, ADParams, ADGrid = ADProb.sol, ADProb.clock, ADProb.vars, ADProb.params, ADProb.grid
    end

    TracerAdvDiff_QG.QGset_c!(ADProb, IC.C₀)
    ADOutput = Output(ADProb, filename, (Symbol(:ConcentrationProb, i), GetConcentration))
    
    #Simulation loop
    while ADClock.step <= nsteps

        if ADClock.step % 1000 == 0
            println("Step number: ", round(Int, ADClock.step))
        end
        saveoutput(ADOutput)
        stepforward!(ADProb, nsubs)
        TracerAdvDiff_QG.QGupdatevars!(ADProb)
        TracerAdvDiff_QG.vel_field_update!(ADProb, QGProb, nsubs)

    end

    saveproblem(ADOutput)
    #Save the number of steps in the simulation
    jldopen(ADOutput.path, "a+") do path
        path["clock/nsteps"] = ADClock.step - 1
        path["save_freq"] = nsubs
        if delay_time != 0
            path["delay_time"] = delay_time
        end
    end

end