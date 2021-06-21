#Tracer advcetion diffusion on large domain

#Experiment setup script. Depending on domain size and time step the title, xticks and yticks need to be manually done.
include("ExperimentSetup.jl")

#Import a flow that has already been set up.
#include("Flows/ExampleFlow.jl")
include("Flows/FlowSetup_500domain_res128.jl")
#include("Flows/FlowSetup_1000domain_res128.jl")
#include("Flows/FlowSetup_1000domain_res256.jl")
#include("Flows/FlowSetup_1500domain_res128.jl")
#include("Flows/FlowSetup_1500domain_res256.jl")

κ = 0.01
#Set delay time (that is flow for some number of days, then drop tracer in)
delay_time = 0
#Set the tracer advection probelm by passing in the QG problem 
ADProb = TracerAdvDiff_QG.Problem(;prob = QGProb, delay_time = delay_time, nsubs = nsubs, κ = κ)
ADSol, ADClock, ADVars, ADParams, ADGrid = ADProb.sol, ADProb.clock, ADProb.vars, ADProb.params, ADProb.grid
ADx, ADy = gridpoints(ADGrid)
x, y = ADGrid.x, ADGrid.y

#Set the Gaussian blob initial condition
μIC = [0, 0]
Σ = [1 0; 0 1]
IC = GaussianBlobIC(μIC, Σ, ADGrid)

TracerAdvDiff_QG.QGset_c!(ADProb, IC.C₀)

max_conc = [findmax(ADVars.c[:, :, 1])[1], findmax(ADVars.c[:, :, 2])[1]]

#Define frequency at which to save a plot.
#plot_time_AD is when to get the first plot, plot_time_inc is at what interval subsequent plots are created.
#Setting them the same gives plots at equal time increments. (Might be a better work around)
plot_time_AD, plot_time_inc = 0.2, 0.2
#Blank array to save the step number so can plot histogram corresponding to the tracer advection plot.
step_nums = []

#Create a file to save data to
filename = CreateFile(ADProb)

#Save parameters from the ADProb type so have the info
jldopen(filename, "a+") do file
    file["ADParameters"] = ADParams
    file["QGParameters"] = QGParams
    file["InitialCondition"] = IC
    file["MaxConcentration"] = max_conc
end


#Set arguments for the plots from the ADProb
plotargs = Set_plotargs(ADProb)
#Step the tracer advection problem forward and plot at the desired time step.
while ADClock.step <= nsteps
    if ADClock.step == 0
        tp_u = heatmap(x, y, ADVars.c[:, :, 1]'; plotargs...)
        push!(UpperLayerTracerPlots, tp_u)
        tp_l = heatmap(x, y, ADVars.c[:, :, 2]'; plotargs...)
        push!(LowerLayerTracerPlots, tp_l)
        push!(step_nums, ADClock.step)
    elseif round(Int64, ADClock.step) == round(Int64, plot_time_AD*nsteps)
        tp_u = heatmap(x, y, ADVars.c[:, :, 1]'; plotargs...)
        push!(UpperLayerTracerPlots, tp_u)
        tp_l = heatmap(x, y, ADVars.c[:, :, 2]'; plotargs...)
        push!(LowerLayerTracerPlots, tp_l)
        push!(step_nums, ADClock.step)
        global plot_time_AD += plot_time_inc
    end
    MeasureMixing.fit_hist!(filename, ADProb, number_of_bins = 30)
    stepforward!(ADProb, nsubs)
    TracerAdvDiff_QG.QGupdatevars!(ADProb)
    TracerAdvDiff_QG.vel_field_update!(ADProb, QGProb, nsubs)
end

#Save the created plots to the .jld2 file
jldopen(filename, "a+") do file
    file["TracerPlots"] = [UpperLayerTracerPlots, LowerLayerTracerPlots]
    file["SaveStepsforTracerPlots"] = step_nums
end
