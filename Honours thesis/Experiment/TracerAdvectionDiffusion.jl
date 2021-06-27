#Tracer advcetion diffusion experiment

#Change to the correct directory (if it was not already correct for some reason)
SimPath = joinpath("/Users/Joey/Desktop/ThesisCode/QG_tracer_advection", "Honours thesis/Experiment")
cd(SimPath)
#Load in all the required packages for the simulation
include("PackageSetup.jl")

#Import a flow that has already been set up from the Flows folder
#include("Flows/ExampleFlow.jl")
#include("Flows/FlowSetup_500domain_res128.jl")
#include("Flows/FlowSetup_1000domain_res128.jl")
#include("Flows/FlowSetup_1000domain_res256.jl")
#include("Flows/FlowSetup_1500domain_res128.jl")
#include("Flows/FlowSetup_1500domain_res256.jl")
#include("Flows/FlowSetup_nondim_32domain_128res.jl")
include("Flows/FlowSetup_nondim_32domain_128res.jl")

κ = 0.01
#Set delay time (that is flow for some number of days, then drop tracer in)
delay_time = 0
#Set the tracer advection probelm by passing in the QG problem 
ADProb = TracerAdvDiff_QG.Problem(;prob = QGProb, delay_time = delay_time, nsubs = nsubs, κ = κ)
ADSol, ADClock, ADVars, ADParams, ADGrid = ADProb.sol, ADProb.clock, ADProb.vars, ADProb.params, ADProb.grid

#Set the Gaussian blob initial condition
μIC = [0, 0]
Σ = [1 0; 0 1]
IC = GaussianBlobIC(μIC, Σ, ADGrid)

TracerAdvDiff_QG.QGset_c!(ADProb, IC.C₀)

save_freq = 200 #Frequency at which to save data
filename = CreateFile(ADProb, save_freq, SimPath)
ADOutput = Output(ADProb, filename, (:Concentration, GetConcentration))
saveproblem(ADOutput)

#Simulation loop
while ADClock.step <= nsteps

    if ADClock.step % 1000 == 0
        println("Step number: ", round(Int, ADClock.step))
    end
    if ADClock.step % save_freq == 0
        saveoutput(ADOutput)
    end
    stepforward!(ADProb, nsubs)
    TracerAdvDiff_QG.QGupdatevars!(ADProb)
    TracerAdvDiff_QG.vel_field_update!(ADProb, QGProb, nsubs)

end

#Save the number of steps in the simulation
jldopen(ADOutput.path, "a+") do path
    path["clock/nsteps"] = ADClock.step - 1
    path["save_freq"] = save_freq
end