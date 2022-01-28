# A script similar to that in the notebook `gauss-diff.ipynb` from `Gaussian-diffusion` repo

SimPath = joinpath(pwd(), "Honours thesis/Experiment")
cd(SimPath)
#Load in all the required packages for the simulation
include("PackageSetup.jl")

## Diffusion problem - no need to run again as data is saved.
dev = CPU() 
nx = 128          
stepper = "RK4"         
dt = 0.002          
nsteps = 7000            
nsubs  = 50  
Lx = 16
κ = 0.25
uvel(x, y) = 1
vvel(x, y) = 0

ADProb = TracerAdvDiff_QG.Problem(; nx=nx, Lx=Lx, κ=κ,
                                #steadyflow=true, u=uvel, v=vvel,
                                dt=dt, stepper=stepper)
ADSol, ADClock, ADVars, ADParams, ADGrid = ADProb.sol, ADProb.clock, ADProb.vars, ADProb.params, ADProb.grid

#Set the Gaussian blob initial condition
μIC = [0, 0]
Σ = [1 0; 0 1]
IC = GaussianBlobIC(μIC, Σ, ADGrid)

#Gaussian blob
set_c!(ADProb, IC.C₀)

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
    TracerAdvDiff_QG.updatevars!(ADProb)

end

#Save the number of steps in the simulation
jldopen(ADOutput.path, "a+") do path
    path["clock/nsteps"] = ADClock.step - 1
    path["save_freq"] = save_freq
    path["params/nlayers"] = 1
end

## Plots of diffusion problem

file = joinpath(pwd(), "Output/Simulation: Lx̂ = Lŷ = 16, nx = 128, save_freq = 50, IC = GaussianBlob, Ensemble = false/SimulationData.jld2")
data = load(file)

# Tracer plots
tracer_plots = tracer_plot(data)
plot(tracer_plots..., size = (1200, 1200))

conc_area_plot = concarea_plot(data)
plot(conc_area_plot..., size = (1200, 1200))
plot(conc_area_plot[1], 
    xlabel = "A",
    title = "Initial ordered concentration")
plot(conc_area_plot[end], 
    xlabel = "A",
    title = "Final ordered concentration")

# Average area growth plot
t = time_vec(data)
first_moms = first_moment(data)

plot(t, first_moms,
 xlabel = "t",
 ylabel = "⟨A⟩(t)",
 label = false,
 title = "Growth of average area of tracer\npatch during diffusion experiment"
)

Δt = t[end] - t[1]
ΔA = first_moms[end] - first_moms[1]
κ_est = ΔA / (4 * π * Δt) 
