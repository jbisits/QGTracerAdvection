#Tracer advcetion diffusion on large domain
using .TracerAdvDiff_QG

using GeophysicalFlows.MultiLayerQG, Plots, Distributions, StatsBase, LinearAlgebra, JLD2

#Import the flow that has already been set up.
include("flow_setup.jl")

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
blob = MvNormal(μIC, Σ)
blob_IC(x, y) = pdf(blob, [x, y])
C₀ = @. blob_IC(ADx, ADy)

TracerAdvDiff_QG.QGset_c!(ADProb, C₀)

max_conc = [findmax(ADVars.c[:, :, 1])[1], findmax(ADVars.c[:, :, 1])[1]]

#Define blank arrays in which to store the plots of tracer diffusion in each layer.
lower_layer_tracer_plots_AD = Plots.Plot{Plots.GRBackend}[]
upper_layer_tracer_plots_AD = Plots.Plot{Plots.GRBackend}[]

#Define frequency at which to save a plot.
#plot_time_AD is when to get the first plot, plot_time_inc is at what interval subsequent plots are created.
#Setting them the same gives plots at equal time increments. (Might be a better work around)
plot_time_AD, plot_time_inc = 0.2, 0.2
#Blank array to save the step number so can plot histogram corresponding to the tracer advection plot.
step_nums = []

filepath = pwd()
filename = joinpath(filepath, "Honours thesis/Experimet/hist_conc.jld2")

#Step the tracer advection problem forward and plot at the desired time step.
while ADClock.step <= nsteps
    if ADClock.step == 0
        tp_u = heatmap(x, y, ADVars.c[:, :, 1]',
                aspectratio = 1,
                c = :deep,
                xlabel = "x",
                ylabel = "y",
                colorbar = true,
                xlims = (-ADGrid.Lx/2, ADGrid.Lx/2),
                ylims = (-ADGrid.Ly/2, ADGrid.Ly/2),
                title = "C(x,y,t), day = "*string(ADClock.t/3600));
        push!(upper_layer_tracer_plots_AD, tp_u)
        tp_l = heatmap(x, y, ADVars.c[:, :, 2]',
                aspectratio = 1,
                c = :deep,
                xlabel = "x",
                ylabel = "y",
                colorbar = true,
                xlims = (-ADGrid.Lx/2, ADGrid.Lx/2),
                ylims = (-ADGrid.Ly/2, ADGrid.Ly/2),
                title = "C(x,y,t), day = "*string(ADClock.t/3600))
        push!(lower_layer_tracer_plots_AD, tp_l)
        push!(step_nums, ADClock.step)
    elseif round(Int64, ADClock.step) == round(Int64, plot_time_AD*nsteps)
        tp_u = heatmap(x, y, ADVars.c[:, :, 1]',
                aspectratio = 1,
                c = :deep,
                xlabel = "x",
                ylabel = "y",
                colorbar = true,
                xlims = (-ADGrid.Lx/2, ADGrid.Lx/2),
                ylims = (-ADGrid.Ly/2, ADGrid.Ly/2),
                title = "C(x,y,t), day = "*string(ADClock.t/3600))
        push!(upper_layer_tracer_plots_AD, tp_u)
        tp_l = heatmap(x, y, ADVars.c[:, :, 2]',
                aspectratio = 1,
                c = :deep,
                xlabel = "x",
                ylabel = "y",
                colorbar = true,
                xlims = (-ADGrid.Lx/2, ADGrid.Lx/2),
                ylims = (-ADGrid.Ly/2, ADGrid.Ly/2),
                title = "C(x,y,t), day = "*string(ADClock.t/3600))
        push!(lower_layer_tracer_plots_AD, tp_l)
        push!(step_nums, ADClock.step)
        global plot_time_AD += plot_time_inc
    end
    #MeasureMixing.fit_hist!(filename, ADProb, number_of_bins = 30)
    stepforward!(ADProb, nsubs)
    TracerAdvDiff_QG.QGupdatevars!(ADProb)
    #Updates the velocity field in advection problem to the velocity field in the MultiLayerQG.Problem at each timestep.
    TracerAdvDiff_QG.vel_field_update!(ADProb, QGProb, nsubs)
end
#Need to set this up so this does not need to be hardcoded.
#Display the tracer advection in the upper layer.
plot_top = plot(upper_layer_tracer_plots_AD[1], upper_layer_tracer_plots_AD[2], 
                upper_layer_tracer_plots_AD[3], upper_layer_tracer_plots_AD[4],
                upper_layer_tracer_plots_AD[5], upper_layer_tracer_plots_AD[6])
     
#Display the tracer advection in the lower layer.
plot_bottom = plot(lower_layer_tracer_plots_AD[1], lower_layer_tracer_plots_AD[2], 
                   lower_layer_tracer_plots_AD[3], lower_layer_tracer_plots_AD[4],
                   lower_layer_tracer_plots_AD[5], lower_layer_tracer_plots_AD[6])

#Now can load the output of the .jld2 file created to create histograms and plots of Concentration ~ normalised area.

data = load("Experiment/hist_conc.jld2")

hist_top = Plots.Plot{Plots.GRBackend}[]
hist_bottom = Plots.Plot{Plots.GRBackend}[]
for i ∈ step_nums
    push!(hist_top, plot(data["Histograms/step"*string(i)][1], 
                             label = false, 
                            xlabel = "Concentration", 
                            ylabel = "Normalised area",
                             xlims = (0, max_conc[1] + 0.01)))
    push!(hist_bottom, plot(data["Histograms/step"*string(i)][2],
                             label = false, 
                            xlabel = "Concentration", 
                            ylabel = "Normalised area",
                             xlims = (0, max_conc[2] + 0.01)))
end
hist_top = plot(hist_top[1], hist_top[2], hist_top[3], hist_top[4], hist_top[5], hist_top[6])
hist_bottom = plot(hist_bottom[1], hist_bottom[2], hist_bottom[3], hist_bottom[4], hist_bottom[5], hist_bottom[6])

plot(plot_top, hist_top, layout=(2, 1), size=(1200, 1200))
plot(plot_bottom, hist_bottom, layout=(2, 1), size=(1200, 1200))

#Concentration ~ area example plot
plot(data["ConcentrationData/step400"][1], data["Histograms/step400"][1].edges,
         label = false,
        xlabel = "Normalised area",
        ylabel = "Concentration",
        xlims = (0, max_conc[1] + 0.01)
        )

#Animation of #Concentration ~ area 
#This animation is already saved so is commented out
#=
ConcVsArea = @animate for i in 1:10:nsteps
    p1 = plot(data["ConcentrationData/step"*string(i)][1], data["Histograms/step"*string(i)][1].edges,
                 label = false,
                xlabel = "Normalised area",
                ylabel = "Concentration",
                 xlims = (0, max_conc[1] + 0.01),
                 title = "Top layer"
                )
    p2 = plot(data["ConcentrationData/step"*string(i)][2], data["Histograms/step"*string(i)][2].edges,
                 label = false,
                xlabel = "Normalised area",
                ylabel = "Concentration",
                 xlims = (0, max_conc[2] + 0.01),
                 title = "Bottom layer"
                )
    plot(p1, p2)
end

mp4(ConcVsArea, "Movies/ConcVsArea.mp4", fps=18)
=#