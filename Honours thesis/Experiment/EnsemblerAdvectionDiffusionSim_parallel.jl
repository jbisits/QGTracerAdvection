#Ensemble tracer advection diffusion experiment, trying to get one to work in parallel.
using Distributed
addprocs(2)
#Change to the correct directory (if it was not already correct for some reason)
begin @everywhere SimPath = joinpath("/Users/Joey/Desktop/ThesisCode/QG_tracer_advection", "Honours thesis/Experiment")
    cd(SimPath)
    #Load in all the required packages for the simulation
    include("PackageSetup.jl")

    #Define number of tracer advection simulations
    ADSims = 2
    #Import a flow that has already been set up from the Flows folder. For an ensemble use an array of flows
    include("Flows/EnsembleSquare/EnsembleFlow_32domain_64res.jl")

    nsubs  = 1           #Set the number of steps the simulation takes at each iteration. This is also the frequency that data is saved at.         
    nsteps = 1000         #Set the total amount of time steps the advection-diffusion simulation should run for

    κ = 0.01
    #Set delay times (that is flow for some length of time, then drop tracer in)
    delay_time = Δt̂ * 3000
    #Set the frequency at which to save data
    save_freq = 100
end
#This a parallel simulation where an array of QG Problems are defined then step forward different AD Problems using these flows
Threads.@threads for i ∈ 1:ADSims

    ADProb = TracerAdvDiff_QG.Problem(;prob = QGProbs[i], delay_time = delay_time, nsubs = nsubs, κ = κ)
    ADSol, ADClock, ADVars, ADParams, ADGrid = ADProb.sol, ADProb.clock, ADProb.vars, ADProb.params, ADProb.grid
    #Initial condition
    μIC = 0
    σ² = 1
    IC = GaussianStripIC(μIC, σ², ADGrid)
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