#Tracer advcetion diffusion experiment

#Change to the correct directory (if it was not already correct for some reason)
SimPath = joinpath("/Users/Joey/Desktop/ThesisCode/QG_tracer_advection", "Honours thesis/Experiment")
cd(SimPath)
#Load in all the required packages for the simulation
include("PackageSetup.jl")

#Import a flow that has already been set up from the Flows folder
#include("Flows/FlowSetup_nondim_32domain_64res.jl")
#include("Flows/FlowSetup_nondim_64domain_128res.jl")
#include("Flows/FlowSetup_nondim_128domain_256res.jl")

nsubs  = 200            #Set the number of steps the simulation takes at each iteration. This is also the frequency that data is saved at.         
nsteps = 22000          #Set the total amount of time steps the advection-diffusion simulation should run for

κ = 6e-4
#Set delay time (that is flow for some length of time, then drop tracer in)
delay_time = Δt̂ * 3000
#Set the tracer advection probelm by passing in the QG problem 
ADProb = TracerAdvDiff_QG.Problem(;prob = QGProb, delay_time = delay_time, nsubs = nsubs, κ = κ)
ADSol, ADClock, ADVars, ADParams, ADGrid = ADProb.sol, ADProb.clock, ADProb.vars, ADProb.params, ADProb.grid

#Set the Gaussian blob initial condition
μIC = [0, 0]
Σ = [1 0; 0 1]
IC = GaussianBlobIC(μIC, Σ, ADGrid)

QGset_c!(ADProb, IC.C₀)

filename = CreateFile(ADProb, IC, nsubs, SimPath)
ADOutput = Output(ADProb, filename, (:Concentration, GetConcentration))
saveproblem(ADOutput)

#Simulation loop
while ADClock.step <= nsteps

    if ADClock.step % 1000 == 0
        println("Step number: ", round(Int, ADClock.step))
    end
    saveoutput(ADOutput)
    stepforward!(ADProb, nsubs)
    QGupdatevars!(ADProb)
    vel_field_update!(ADProb, QGProb, nsubs)

end

#Save the number of steps in the simulation
jldopen(ADOutput.path, "a+") do path
    path["clock/nsteps"] = ADClock.step - 1
    path["save_freq"] = nsubs
    if delay_time != 0
        path["delay_time"] = delay_time
    end
end