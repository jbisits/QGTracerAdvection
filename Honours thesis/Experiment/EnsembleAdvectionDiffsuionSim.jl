#Ensemble tracer advection diffusion experiment

## Load all packages
#Change to the correct directory (if it was not already correct for some reason)
SimPath = joinpath("/Users/Joey/Desktop/ThesisCode/QG_tracer_advection", "Honours thesis/Experiment")
cd(SimPath)
#Load in all the required packages for the simulation
include("PackageSetup.jl")

## Run a simulation
#Define number of tracer advection simulations
ADSims = 10

#Import a an ensemble of flows on a square domain
include("Flows/EnsembleSquare/EnsembleFlow_32domain_64res.jl")
#include("Flows/EnsembleSquare/EnsembleFlow_64domain_128res.jl")
#include("Flows/EnsembleSquare/EnsembleFlow_128domain_256res.jl")
#include("Flows/EnsembleSquare/EnsembleFlow_256domain_512res.jl")

#Import an ensemble of flows on a rectanglular domain
#include("Flows/EnsembleRectangle/EnsembleFlow_32_64_domain.jl")
#include("Flows/EnsembleRectangle/EnsembleFlow_32_128_domain.jl")
#include("Flows/EnsembleRectangle/EnsembleFlow_32_256_domain.jl")
#include("Flows/EnsembleRectangle/EnsembleFlow_64_128_domain.jl")
#include("Flows/EnsembleRectangle/EnsembleFlow_64_256_domain.jl")

nsubs  = 1           #Set the number of steps the simulation takes at each iteration. This is also the frequency that data is saved at.         
nsteps = 3000           #Set the total amount of time steps the advection-diffusion simulation should run for

κ = 0.01
#Set delay times (that is flow for some length of time, then drop tracer in)
delay_time = Δt̂ * 4500
#delay_time = 0
#Set the frequency at which to save data
save_freq = 50

#This runs a non-parallel simulation where an array of QG problems is defined then used to advect the tracers
for i ∈ 1:ADSims

    ADProb = TracerAdvDiff_QG.Problem(;prob = QGProbs[i], delay_time = delay_time, nsubs = nsubs, κ = κ)
    ADSol, ADClock, ADVars, ADParams, ADGrid = ADProb.sol, ADProb.clock, ADProb.vars, ADProb.params, ADProb.grid
    #Set the Gaussian blob initial condition
    μIC = [0, 0]
    Σ = [1 0; 0 1]
    IC = GaussianBlobIC(μIC, Σ, ADGrid)

    #Set the Gaussian strip initial condition
    #μIC = 0
    #σ² = 1
    #IC = GaussianStripIC(μIC, σ², ADGrid)

    #Set QGPV as initial condition
    #IC = QGPVIC(QGProbs[i])

    #File name for saving, FourierFlows creates a new file each time with _i appended
    filename = CreateFile(ADProb, IC, save_freq, SimPath; Ensemble = true)

    QGset_c!(ADProb, IC.C₀)
    ADOutput = Output(ADProb, filename, (:Concentration, GetConcentration))
    saveproblem(ADOutput)
    
    #Simulation loop
    while ADClock.step <= nsteps

        if ADClock.step % save_freq == 0
            saveoutput(ADOutput)
            println("Advection problem: ", i, ", step number: ", round(Int, ADClock.step), ", saved data")
        end
        stepforward!(ADProb, nsubs)
        QGupdatevars!(ADProb)
        vel_field_update!(ADProb, QGProbs[i], nsubs)

    end

    #Save the number of steps and other info from the simulation
    jldopen(ADOutput.path, "a+") do path
        path["clock/nsteps"] = ADClock.step - 1
        path["save_freq"] = save_freq
        path["no_of_sims"] = ADSims
        if delay_time != 0
            path["delay_time"] = delay_time
        end
    end

end