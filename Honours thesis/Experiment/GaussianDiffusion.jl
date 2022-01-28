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
########################################################################################################

# Analysis of data

## Plots of diffusion problem

file = joinpath(pwd(), "Output/Simulation: Lx̂ = Lŷ = 16, nx = 128, save_freq = 50, IC = GaussianBlob, Ensemble = false/SimulationData.jld2")
data = load(file)

# Tracer plots
tracer_plots = tracer_plot(data)
plot(tracer_plots..., size = (1200, 1200))

inital_tracer = plot(tracer_plots[1],
                    colorbar_title = " \nConcentration",
                    thickness_scaling = 1.6,
                    right_margin = 2Plots.mm,
                    title = "Initial tracer \nconcentration over grid")
final_tracer = plot(tracer_plots[end],
                    colorbar_title = "  \nConcentration",
                    thickness_scaling = 1.5,
                    right_margin = 3Plots.mm,
                    title = "Final tracer \nconcentration over grid")

# Ordered concentration plots
conc_area_plot = concarea_plot(data)
plot(conc_area_plot..., size = (1200, 1200))

initial_conc = plot(conc_area_plot[1], 
                    xlabel = "A",
                    title = "Initial ordered concentration")
final_conc = plot(conc_area_plot[end], 
                    xlabel = "A",
                    title = "Final ordered concentration")


full_tracer_plot = plot(inital_tracer, final_tracer, initial_conc, final_conc, size = (800, 800))
savefig(full_tracer_plot, "diff_expt.png")
## Average area growth plot
t = time_vec(data)
first_moms = first_moment(data)

av_area_plot = plot(t, first_moms,
                xlabel = "t",
                ylabel = "⟨A⟩(t)",
                label = false,
                title = "Growth of average area of tracer\npatch during diffusion experiment"
                )

savefig(av_area_plot, "diff_av_area.png")
Δt = t[end] - t[1]
ΔA = first_moms[end] - first_moms[1]
κ_est = ΔA / (4 * π * Δt) 
