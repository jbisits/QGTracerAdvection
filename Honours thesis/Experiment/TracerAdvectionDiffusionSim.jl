#Tracer advcetion diffusion experiment

#Change to the correct directory (if it was not already correct for some reason)
SimPath = joinpath("/Users/Joey/Desktop/ThesisCode/QG_tracer_advection", "Honours thesis/Experiment")
cd(SimPath)
#Load in all the required packages for the simulation
include("PackageSetup.jl")

#Import a flow that has already been set up from the Flows folder

# Square domain flows
include("Flows/Square/FlowSetup_nondim_32domain_64res.jl")
#include("Flows/Square/FlowSetup_nondim_64domain_128res.jl")
#include("Flows/Square/FlowSetup_nondim_128domain_256res.jl")
#include("Flows/Square/FlowSetup_nondim_256domain_512res.jl")

# Rectangular domain flows (longer in meridional direction)
#include("Flows/Rectangle/FlowSetup_nondim_32_64domain.jl")
#include("Flows/Rectangle/FlowSetup_nondim_32_128domain.jl")
#include("Flows/Rectangle/FlowSetup_nondim_32_256domain.jl")
#include("Flows/Rectangle/FlowSetup_nondim_64_128domain.jl")
#include("Flows/Rectangle/FlowSetup_nondim_64_256domain.jl")

nsubs  = 1            #Set the number of steps the simulation takes at each iteration.         
nsteps = 6000          #Set the total amount of time steps the advection-diffusion simulation should run for

κ = 0.01
#Set delay time (that is flow for some length of time, then drop tracer in)
#delay_time = Δt̂ * 3000
delay_time = 0
#Set the tracer advection probelm by passing in the QG problem 
ADProb = TracerAdvDiff_QG.Problem(;prob = QGProb, delay_time = delay_time, nsubs = nsubs, κ = κ)
ADSol, ADClock, ADVars, ADParams, ADGrid = ADProb.sol, ADProb.clock, ADProb.vars, ADProb.params, ADProb.grid

#Set the Gaussian blob initial condition
#μIC = [0, 0]
#Σ = [1 0; 0 1]
#IC = GaussianBlobIC(μIC, Σ, ADGrid)
#Set the Gaussian band initial condition
#μIC = 0
#σ² = 1
#IC = GaussianStripIC(μIC, σ², ADGrid)
#IC = PointSourceIC([64, 64], 1, ADGrid)
#Set the QGPV initial condition
IC = QGPVIC(QGProb)

QGset_c!(ADProb, IC.C₀)
save_freq = 50
filename = CreateFile(ADProb, IC, save_freq, SimPath)
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

#Save the number of steps in the simulation
jldopen(ADOutput.path, "a+") do path
    path["clock/nsteps"] = ADClock.step - 1
    path["save_freq"] = save_freq
    if delay_time != 0
        path["delay_time"] = delay_time
    end
end